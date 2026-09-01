# -*- coding: utf-8 -*-
r"""Checks what an adapter *does*, not where its constants came from.

``test_bench_pipeline.py`` guards provenance: who owns a clock, whether a
version is pinned, whether a constant was typed.  The adversarial pass showed
that provenance is not enough.  A throwaway adapter that computed its whole
answer in ``setup()``, returned a stored scalar from ``evaluate()``, discarded
delta_CP in ``configure()``, claimed a batched entry point it never called, and
reached the matter constant as ``float('1.5' '2588e-4')`` --- which
``_code_only`` cannot see, because it strips string literals before the check
runs --- passed all twelve of those invariants.

So these three observe behaviour instead:

* the timed scan must actually move something;
* ``evaluate()`` must consume the grid it was given;
* a precision knob must reach different answers at different settings.

Each one fails on that throwaway adapter.  None of them reads a clock: a
benchmark suite that benchmarks during collection cannot be trusted afterwards,
so every check here is a correctness probe on a deliberately tiny grid.
"""

import importlib
import pathlib
import sys

import pytest


BENCH = pathlib.Path(__file__).resolve().parent / 'bench'
sys.path.insert(0, str(BENCH))
sys.path.insert(0, str(BENCH / 'adapters'))

#: Deliberately tiny.  These are correctness probes, and a probe that takes
#: long enough to matter would disturb the machine the real runs need idle.
PROBE_ENERGIES = 2
PROBE_ZENITH = 1


def _runner():
    return importlib.import_module('runner')


def _adapter(code):
    r"""Returns an adapter instance, or skips if its code is not installed.

    nuSQuIDS and nuCraft come from ``build.sh``; a checkout without it must
    skip rather than fail, but must never pass silently.
    """
    if code == 'nuoscprobexact':
        import fastkernels
        if not (fastkernels.HAVE_NUMBA and fastkernels.available()
                and getattr(fastkernels, 'USE_NUMBA', False)):
            pytest.skip(
                'the compiled kernels are not live, and the adapter refuses to '
                'run on the NumPy path so that no measurement of it can be '
                'published as this library; these probes read no clock, but '
                'they reach that guard through setup()')
    runner = _runner()
    try:
        return runner.load_adapter(code)
    except (ImportError, ModuleNotFoundError) as exc:      # noqa: PERF203
        pytest.skip('%s not installed (%s); run tests/bench/build.sh'
                    % (code, exc))


def _problem(grid, knob=0, n_e=PROBE_ENERGIES, n_z=PROBE_ZENITH):
    return _runner().build_problem(grid, knob, n_e, n_z)


def _codes():
    return sorted(_runner().CODES)


def _grid_for(code):
    r"""nuCraft propagates through its Earth model only, so it gets a chord."""
    return 'CHORD/12x1'


@pytest.mark.parametrize('code', ['nuoscprobexact', 'nusquids', 'nucraft'])
def test_evaluate_responds_to_configure(code):
    r"""The AMORTIZED scan must move something inside the timed region.

    ``configure(dcp)`` is the one thing the protocol varies, and the whole
    protocol is meaningless if an adapter ignores it: the harness would time a
    loop over a constant.  The throwaway adapter's ``configure`` was ``pass``,
    and nothing in the twelve invariants noticed.
    """
    driver = _adapter(code)
    problem = _problem(_grid_for(code))
    driver.setup(problem)

    driver.configure(problem.dcp)
    first = driver.evaluate()
    driver.configure(problem.dcp + 0.7)          # well clear of round-off
    second = driver.evaluate()

    assert first != second, (
        '%s: evaluate() returned %r at two different delta_CP values, so the '
        'timed scan moves nothing' % (code, first))


@pytest.mark.parametrize('code', ['nuoscprobexact', 'nusquids', 'nucraft'])
def test_evaluate_consumes_the_whole_grid(code):
    r"""The checksum must depend on every point the harness normalises by.

    ``us_per_point`` divides by ``problem.points()``.  An adapter that solves a
    subset --- or a stored constant --- is divided by a grid it never touched,
    and the result is a small number that looks like speed.
    """
    driver_small = _adapter(code)
    small = _problem(_grid_for(code), n_e=PROBE_ENERGIES)
    driver_small.setup(small)
    driver_small.configure(small.dcp)
    checksum_small = driver_small.evaluate()

    driver_large = _adapter(code)
    large = _problem(_grid_for(code), n_e=2*PROBE_ENERGIES)
    driver_large.setup(large)
    driver_large.configure(large.dcp)
    checksum_large = driver_large.evaluate()

    assert large.points() == 2*small.points(), 'the probe grids are not sized as intended'
    assert checksum_small != checksum_large, (
        '%s: the checksum is %r on both a %d-point and a %d-point grid, so '
        'evaluate() is not consuming the grid it is normalised by'
        % (code, checksum_small, small.points(), large.points()))


@pytest.mark.parametrize('code', ['nuoscprobexact', 'nusquids', 'nucraft'])
def test_knob_extremes_are_distinguishable(code):
    r"""A precision knob must actually reach different precisions.

    Objection LBL-3 --- the strongest point raised against the comparison ---
    was that a setting existed which was never tried.  Declaring a domain and
    sweeping it is worthless if the endpoints compute the same thing, which is
    what a knob that is read but never applied looks like from outside.
    """
    driver = _adapter(code)
    domain = [k for k in driver.capabilities()['knob_domain']]
    assert len(domain) >= 2, '%s declares a knob domain with nothing to sweep' % code

    results = {}
    for knob in (domain[0], domain[-1]):
        instance = _adapter(code)
        problem = _problem(_grid_for(code), knob=knob)
        instance.setup(problem)
        instance.configure(problem.dcp)
        results[knob] = instance.evaluate()

    lo, hi = domain[0], domain[-1]
    assert results[lo] != results[hi], (
        '%s: knob %r and knob %r give the identical checksum %r, so the knob '
        'is declared but not applied' % (code, lo, hi, results[lo]))


def test_the_compiled_adapters_expose_the_accuracy_path():
    r"""Every compiled adapter must define ``probabilities()``.

    ``bench.hpp`` calls it from ``main``, so an adapter without one does not
    link --- which is the mechanism, not a convention.  This states it as an
    invariant too, because the link error is only felt by whoever next runs
    ``build.sh``, and a missing accuracy path is how a speed number ends up
    with nothing to check it against.
    """
    missing = []
    for path in sorted((BENCH / 'adapters').glob('*.cpp')):
        if 'void probabilities(' not in path.read_text():
            missing.append(path.name)
    assert not missing, (
        'compiled adapters with no accuracy path: %s' % ', '.join(missing))
