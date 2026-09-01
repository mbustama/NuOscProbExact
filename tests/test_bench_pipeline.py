# -*- coding: utf-8 -*-
r"""Checks the benchmark pipeline's fairness invariants.

An audit of the previous comparison found nine ways it was unfair or
unreproducible, and every one of them was a lapse rather than a decision:
setters inside a timed loop, a batched interface the paper said did not exist,
three coexisting values of one constant, a driver that was never committed.
Care did not prevent them and will not.

So each invariant is a test here.  These do not check that the numbers are
right --- that is what the references are for --- they check that the machinery
*cannot* produce an unfair number: that no driver owns a clock, that no adapter
types a physical constant, that every code is handed one parameter set, and
that every objection raised against the comparison has a test somewhere with
its name on it.
"""

import ast
import json
import pathlib
import re

import pytest


BENCH = pathlib.Path(__file__).resolve().parent / 'bench'
TESTS = pathlib.Path(__file__).resolve().parent


def manifest():
    with open(BENCH / 'manifest.json') as handle:
        return json.load(handle)


def adapters():
    directory = BENCH / 'adapters'
    if not directory.is_dir():
        return []
    return sorted(p for p in directory.iterdir()
                  if p.suffix in ('.cpp', '.c', '.py'))



def _code_only(path):
    r"""Returns `path`'s source with comments and string literals removed.

    Mentioning a constant is not using one.  An adapter that says
    ``// GLoBES's km is built on hbar c = 1.97327e-7 eV m`` is documenting where
    its factor comes from, which is exactly what should be encouraged; the first
    version of this test failed such a comment and would have pushed the
    explanation out of the file to keep the check quiet.  Only code is checked.
    """
    text = path.read_text()
    if path.suffix == '.py':
        import io
        import tokenize
        kept = []
        try:
            for token in tokenize.generate_tokens(io.StringIO(text).readline):
                if token.type not in (tokenize.COMMENT, tokenize.STRING):
                    kept.append(token.string)
        except (tokenize.TokenError, IndentationError):
            return text
        return ' '.join(kept)
    text = re.sub(r'/\*.*?\*/', ' ', text, flags=re.S)   # /* ... */
    text = re.sub(r'//[^\n]*', ' ', text)                 # // ...
    text = re.sub(r'"(?:[^"\\]|\\.)*"', ' ', text)       # "..."
    return text


# --------------------------------------------------------------- the manifest
def test_every_code_is_pinned():
    r"""No code may be identified by anything as soft as a month.

    Three of them once were --- the paper's comparison table recorded
    ``2026-08``, the month they were cloned --- which is not a version and
    cannot be rebuilt.  A pin is a commit or a tarball hash, nothing else.
    """
    for code in manifest()['codes']:
        pin = code['pin']
        assert ('commit' in pin) or ('sha256' in pin) or ('repo_commit' in pin), (
            '%s is not pinned to a commit or a hash: %r' % (code['name'], pin))
        if 'commit' in pin and pin['commit'] != '<filled at run time by run_all.py>':
            assert re.fullmatch(r'[0-9a-f]{40}', pin['commit']), (
                '%s: %r is not a full commit sha' % (code['name'], pin['commit']))
        if 'sha256' in pin:
            assert re.fullmatch(r'[0-9a-f]{64}', pin['sha256']), (
                '%s: %r is not a sha256' % (code['name'], pin['sha256']))


def test_every_code_records_how_latest_was_verified():
    r"""Pinning the version the paper used is not the same as pinning the
    latest, and the difference is how a batched interface went unnoticed for a
    release and a half.  Each code carries the date and the method.
    """
    for code in manifest()['codes']:
        check = code.get('latest_check')
        assert check, '%s has no latest_check' % code['name']
        assert check.get('verified_on') and check.get('how')


def test_precision_knobs_include_their_exact_modes():
    r"""A code's best setting has to be in the swept domain.

    The claim that no code was more accurate "at any setting it exposes" was
    falsified by one setting we never tried: NuFast-LBL reaches exact
    eigenvalues at a negative ``N_Newton``, documented in its release notes and
    absent from our sweep.  Any knob whose registry entry names a sentinel must
    have that sentinel in its domain.
    """
    for code in manifest()['codes']:
        knobs = code.get('capabilities', {}).get('precision_knobs', {})
        for name, spec in knobs.items():
            if isinstance(spec, dict) and 'sentinel' in spec:
                domain = spec.get('domain', [])
                for sentinel in spec['sentinel']:
                    assert int(sentinel) in domain, (
                        '%s.%s: sentinel %s is not in the swept domain %r'
                        % (code['name'], name, sentinel, domain))


# ---------------------------------------------------------------- the drivers
def test_no_adapter_owns_a_clock():
    r"""The harness times; a driver supplies physics.

    Enforced rather than trusted, because the previous NuFast-Earth driver
    wrote its own loop and put the engine construction, the Earth object and
    five setters inside it.  ``bench.hpp`` and ``runner.py`` are the only files
    permitted to name a clock.
    """
    forbidden = ('chrono', 'clock_gettime', 'perf_counter', 'process_time',
                 'time.time')
    for path in adapters():
        text = path.read_text()
        for token in forbidden:
            assert token not in text, (
                '%s names %r; only bench.hpp and runner.py may time'
                % (path.name, token))


def test_no_adapter_types_a_physical_constant():
    r"""Constants reach an adapter through the generated header, never by hand.

    Three different values of one NuFast-Earth density factor once coexisted:
    ``0.99209238`` applied in a driver, ``0.9920928`` implied by its own
    comment, ``0.9920938`` printed in the README.  The generator makes a fourth
    impossible.
    """
    patterns = [r'1\.52588', r'1\.5264932', r'7\.63247', r'1\.97327',
                r'5\.06773', r'0\.99209', r'0\.99249',
                r'2\.525e-3', r'7\.39e-5', r'2\.4511']
    for path in adapters():
        text = _code_only(path)
        for pattern in patterns:
            assert not re.search(pattern, text), (
                '%s carries the literal %s in code; it must come from '
                'conversions' % (path.name, pattern))


def test_capabilities_claims_match_the_registry():
    r"""What an adapter says it does must be what the manifest says it does.

    The paper stated that NuFast-LBL exposed no batched entry point while our
    own driver was calling one.  Neither the registry nor the adapter can now
    say that alone.
    """
    by_name = {c['name']: c for c in manifest()['codes']}
    for path in adapters():
        text = path.read_text()
        found = re.search(r'return\s+"([^"]+)"\s*;\s*\}', text) or \
            re.search(r'name\s*=\s*[\'"]([^\'"]+)[\'"]', text)
        if not found:
            continue
        name = found.group(1)
        if name not in by_name:
            continue
        registry = by_name[name]['capabilities']
        if registry.get('batches_energy'):
            assert re.search(r'batches_energy\s*=\s*true|'
                             r'[\'"]batches_energy[\'"]\s*:\s*True', text), (
                '%s: the manifest says this code batches over energy, so the '
                'adapter must claim it and use the batched entry point' % path.name)


# ------------------------------------------------------- one set of parameters
def test_the_shared_parameter_set_has_exactly_one_home():
    r"""Every code is handed the same physics, from one place.

    ``bench.hpp`` briefly carried its own copy of the NuFit values beside the
    Python side's, which is how "the same parameters for all codes" quietly
    stops being true.
    """
    osc = manifest()['oscillation_parameters']
    for key in ('s12sq', 's13sq', 's23sq', 'dcp_deg', 'dmsq21_ev2',
                'dmsq31_ev2'):
        assert key in osc, 'the manifest does not define %s' % key

    header = (BENCH / 'bench.hpp').read_text()
    assert 'OSC_S12SQ' in header, 'bench.hpp must read the shared parameters'
    assert '0.310' not in header and '2.525e-3' not in header, (
        'bench.hpp types a parameter value; it must use the generated macros')


def test_conversion_factors_are_derived_from_the_pinned_sources():
    r"""Each factor is ours-over-theirs, computed from the code's own source.

    Skipped when the pinned sources are absent, since a checkout without
    ``tests/bench/build.sh`` having run cannot verify them --- but never
    silently passing on a wrong value.
    """
    import sys
    sys.path.insert(0, str(BENCH))
    import conversions

    if not pathlib.Path(conversions.BUILD).is_dir():
        pytest.skip('external sources not built; run tests/bench/build.sh')

    for code in ('NuFast-LBL', 'NuFast-Earth', 'Prob3++'):
        theirs, _ = conversions.extract(code)
        expected = conversions.OURS['matter']/theirs
        assert abs(conversions.mass_defect(code) - expected) < 1e-15


# ------------------------------------------------------------- the objections
def test_every_objection_names_a_test_that_exists_or_is_declared_pending():
    r"""No objection may be quietly dropped.

    ``tests/bench/OBJECTIONS.md`` maps each point raised against the comparison
    to the test that answers it.  A name in that table must either exist in the
    suite or appear in the file's own PENDING list, so the gap between what was
    promised and what is built stays visible rather than implied.
    """
    text = (BENCH / 'OBJECTIONS.md').read_text()
    named = set(re.findall(r'`(test_[a-z0-9_]+)`', text))
    assert named, 'OBJECTIONS.md names no tests at all'

    pending = set()
    if 'PENDING' in text:
        after = text.split('PENDING', 1)[1]
        pending = set(re.findall(r'`(test_[a-z0-9_]+)`', after))

    existing = set()
    for path in TESTS.glob('test_*.py'):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.FunctionDef) and node.name.startswith('test_'):
                existing.add(node.name)

    missing = sorted(named - existing - pending)
    assert not missing, (
        'OBJECTIONS.md names tests that neither exist nor are declared '
        'pending: %s' % ', '.join(missing))


# ------------------------------------------------------------- machine state
def test_a_busy_machine_is_refused():
    r"""A run that cannot be believed must write no artifact.

    Timings on this hardware swing by a quarter between sessions, and the audit
    that prompted this pipeline found a published pair whose ratio was noise.
    The rejection rule is therefore part of the machinery, not advice.
    """
    import sys
    sys.path.insert(0, str(BENCH))
    import machine

    steady = [1.00, 1.01, 0.99, 1.02, 0.98]
    ok, why = machine.admissible(steady)
    assert ok, why

    interrupted = [1.00, 2.20, 1.02, 0.99, 1.01]
    ok, why = machine.admissible(interrupted)
    assert not ok and 'CV' in why, why

    assert not machine.admissible([1.0])[0], 'one block cannot be judged'


def test_canary_drift_invalidates_a_session():
    r"""Three readings bound the drift across an hours-long sweep."""
    import sys
    sys.path.insert(0, str(BENCH))
    import machine

    assert machine.canary_drift(1.00, 1.01, 1.02)[0]
    ok, why = machine.canary_drift(1.00, 1.03, 1.20)
    assert not ok and 'drifted' in why, why


def test_our_own_constants_are_derived_not_transcribed():
    r"""Each factor must equal the potential ratio computed from `globaldefs`.

    This is the check that caught the module's own opening bug.
    ``conversions.OURS['matter']`` began life as ``1.514423e-4``, a seven-figure
    transcription of a value from an old driver's comment, and every factor
    derived from it inherited a 2.4e-8 error --- in a module whose entire
    purpose is to stop constants being typed by hand.

    The four external codes express the potential as ``V = k * rhoYe / 2e9``, so
    the same quantity for this library is ``sqrt(2) G_F N_e * 2e9 / rhoYe``.  If
    a factor disagrees with that, something has been transcribed again.
    """
    import math
    import sys

    sys.path.insert(0, str(BENCH))
    sys.path.insert(0, str(TESTS.parent / 'src'))
    import conversions
    import globaldefs as gd

    if not pathlib.Path(conversions.BUILD).is_dir():
        pytest.skip('external sources not built; run tests/bench/build.sh')

    rho_ye = (gd.DENSITY_MATTER_CRUST_G_PER_CM3
              * gd.ELECTRON_FRACTION_EARTH_CRUST)
    ours = math.sqrt(2.0)*gd.GF*gd.NUM_DENSITY_E_EARTH_CRUST

    for code in ('NuFast-LBL', 'NuFast-Earth', 'Prob3++', 'nuCraft'):
        theirs, _ = conversions.extract(code)
        expected = ours/(theirs*rho_ye/2.0e9)
        assert abs(conversions.mass_defect(code) - expected) < 1e-14, (
            '%s: factor %.14f but the potential ratio is %.14f'
            % (code, conversions.mass_defect(code), expected))

    assert abs(conversions.OURS['km_to_inv_ev'] - gd.CONV_KM_TO_INV_EV) == 0.0, (
        'the length conversion must be globaldefs\' own value, not a copy')
