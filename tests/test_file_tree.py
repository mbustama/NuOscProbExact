# -*- coding: utf-8 -*-
r"""Generates the documented file tree, and checks it against the repository.

The repository layout is shown twice --- in ``README.md`` and in
``docs/source/installation.rst`` --- because the two are read in
different places.  Both are now *generated* from `TREE` below rather
than maintained by hand, so there is one place to edit and no way for
the two copies to disagree.

`TREE` is the ordered list of paths, each with the one-line comment that
appears beside it.  The order is the curated reading order, not
alphabetical, and the comments are prose: neither can be derived from
the filesystem, which is why the table exists at all.  What *is* derived
is membership --- `test_tree_matches_git` requires the file entries here
to be exactly ``git ls-files``, so adding or removing a tracked file
fails the suite until the table is updated.

To update both documents after changing `TREE`, or after adding a file::

    python tests/test_file_tree.py --write

That rewrites the block in each document in place, and prints what it
changed.  Running it with no arguments reports whether they are current
without touching anything.
"""

import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
README = os.path.join(ROOT, 'README.md')
INSTALL_RST = os.path.join(ROOT, 'docs', 'source', 'installation.rst')

ROOT_LABEL = 'NuOscProbExact/'

# Column at which a trailing comment starts.  An entry whose box-drawing
# prefix and name already reach it gets two spaces instead, which is what
# the handwritten tree did for the four longest notebook names.
COMMENT_COLUMN = 37

# The tree, in reading order.  A trailing slash marks a directory; a
# `None` comment means the entry is listed without one.  Keep parents
# immediately before their children --- the renderer derives depth and
# the box-drawing connectors from the paths alone.
TREE = [
    ('.github/', 'Continuous integration (GitHub Actions)'),
    ('.github/workflows/', None),
    ('.github/workflows/tests.yml',
     'The suite: five Pythons, all three backends'),
    ('.github/workflows/lint.yml', 'ruff, and the docs build under -W'),
    ('.github/workflows/pages.yml',
     'Builds and deploys the docs to GitHub Pages'),
    ('.github/workflows/publish.yml', 'Publishes to PyPI on a GitHub Release'),
    ('.gitignore', 'Build, cache, and generated-output artefacts'),
    ('CHANGELOG.md', 'Notable changes, rendered as a docs page'),
    ('RELEASE_NOTES.md', 'What changed since the published version'),
    ('LICENSE', 'MIT license'),
    # Not "the file that you are reading": this tree is rendered into
    # `installation.rst` as well, where that was untrue.
    ('README.md', 'Project overview and worked examples'),
    ('pyproject.toml', 'Packaging metadata and pytest configuration'),
    ('examples/', 'Runnable scripts, one per scenario, linked from README.md'),
    ('examples/example_2nu_trivial.py', 'Two-flavor, arbitrary Hamiltonian'),
    ('examples/example_2nu_vacuum.py', 'Two-flavor, oscillations in vacuum'),
    ('examples/example_2nu_vacuum_coeffs.py',
     'Two-flavor, expansion coefficients'),
    ('examples/example_3nu_trivial.py', 'Three-flavor, arbitrary Hamiltonian'),
    ('examples/example_3nu_vacuum.py', 'Three-flavor, oscillations in vacuum'),
    ('examples/example_3nu_vacuum_coeffs.py',
     'Three-flavor, expansion coefficients'),
    ('examples/example_3nu_matter.py', 'Three-flavor, oscillations in matter'),
    ('examples/example_3nu_nsi.py', 'Three-flavor, matter with NSI'),
    ('examples/example_3nu_liv.py', 'Three-flavor, LIV background'),
    ('docs/', 'Sphinx documentation'),
    ('docs/Makefile', '`make html` on Linux and macOS'),
    ('docs/make.bat', '`make html` on Windows'),
    ('docs/requirements.txt', 'Documentation-only dependencies'),
    ('docs/source/', None),
    ('docs/source/conf.py', 'Sphinx configuration'),
    ('docs/source/index.rst', 'Landing page'),
    ('docs/source/installation.rst', 'Requirements, installation, file tree'),
    ('docs/source/quickstart.rst', 'Shortest path to a probability'),
    ('docs/source/recipes.rst',
     'Numerical recipes, with pre-generated figures'),
    ('docs/source/methodology.rst', 'The SU(2), SU(3) and SU(4) expansions'),
    ('docs/source/functions.rst', 'API reference, from the docstrings'),
    ('docs/source/references.rst', 'Bibliography'),
    ('docs/source/refs.bib', 'BibTeX entries for the bibliography'),
    ('docs/source/changelog.rst', 'Includes the root CHANGELOG.md'),
    ('docs/source/_static/', None),
    ('docs/source/_static/nuoscprobexact_logo.png', None),
    ('docs/source/_static/slabs_composition.svg',
     'How slabs compose, drawn for quickstart.rst'),
    ('img/', 'Figures from earlier versions of README.md'),
    ('img/prob_3nu_vacuum_vs_baseline_ee_em_et.png', None),
    ('img/prob_3nu_vacuum_vs_energy_ee_em_et.png', None),
    ('img/gallery/', 'Figures lifted from the notebooks, shown in README.md'),
    ('img/gallery/gallery_anim_cp.png', None),
    ('img/gallery/gallery_anim_earth.png', None),
    ('img/gallery/gallery_anim_slabs.png', None),
    ('img/gallery/gallery_anim_sterile.png', None),
    ('img/gallery/gallery_biprobability.png', None),
    ('img/gallery/gallery_earth.png', None),
    ('img/gallery/gallery_matter.png', None),
    ('img/gallery/gallery_ordering.png', None),
    ('img/gallery/gallery_oscillogram.png', None),
    ('img/gallery/gallery_prem.png', None),
    ('img/gallery/gallery_profiles.png', None),
    ('img/gallery/gallery_sterile.png', None),
    ('img/gallery/gallery_sterile_earth.png', None),
    ('img/gallery/gallery_vacuum.png', None),
    ('notebooks/', 'Worked examples, with their figures stored inline'),
    ('notebooks/01_basics.ipynb', 'Units, one probability, and broadcasting'),
    ('notebooks/02_vacuum_oscillations.ipynb',
     'Against baseline and against energy'),
    ('notebooks/03_matter_nsi_liv.ipynb',
     'Constant-density matter, NSI, and LIV'),
    ('notebooks/04_oscillogram.ipynb', 'Energy-baseline maps in one call'),
    ('notebooks/05_biprobability.ipynb',
     'CP ellipses, in vacuum and in matter'),
    ('notebooks/06_earth_and_prem.ipynb', 'PREM, chord geometry, and slabs'),
    ('notebooks/07_earth_probabilities.ipynb',
     'Through the Earth, and between sites'),
    ('notebooks/08_unusual_density_profiles.ipynb',
     'Castle-wall and other hand-built profiles'),
    ('notebooks/09_performance.ipynb',
     'Looping vs broadcasting, and the backend'),
    ('notebooks/10_paper_figures.ipynb',
     'The two figures from arXiv:1904.12391'),
    ('notebooks/11_exact_vs_approximations.ipynb',
     'Where the textbook formulas break down'),
    ('notebooks/12_ordering_and_octant.ipynb',
     'Normal vs inverted, and the 23 octant'),
    ('notebooks/13_antineutrinos.ipynb',
     'Conjugate and flip, and two ways to slip'),
    ('notebooks/14_solar_and_adiabatic_msw.ipynb',
     'The MSW resonance, and the cost wall'),
    ('notebooks/15_numerical_edge_cases.ipynb',
     'Degeneracies, and what does not go NaN'),
    ('notebooks/16_four_neutrinos.ipynb',
     'A 3+1 sterile state, through the SU(4) expansion'),
    ('notebooks/17_cross_checks.ipynb',
     'Corroboration from nuSQuIDS and Zaglauer-Schwarzer'),
    ('notebooks/18_evolution_operator.ipynb',
     'The operator, and the SU(n) coefficients'),
    ('notebooks/19_animations.ipynb',
     'Four animated scenes, as stills; the reel they came from'),
    ('notebooks/make_notebooks.py', 'Generates and executes all of the above'),
    ('notebooks/coastlines_crude.json',
     'Coastlines for the Earth cutaway, from GSHHS'),
    ('src/', 'The library'),
    ('src/oscprob2nu.py', 'Two-flavor probabilities, SU(2) expansion'),
    ('src/oscprob3nu.py', 'Three-flavor probabilities, SU(3) expansion'),
    ('src/oscprob4nu.py', 'Four-flavor probabilities, SU(4) expansion'),
    ('src/hamiltonians2nu.py', 'Example two-flavor Hamiltonians'),
    ('src/hamiltonians3nu.py', 'Example three-flavor Hamiltonians'),
    ('src/hamiltonians4nu.py', 'Example four-flavor (3+1) Hamiltonians'),
    ('src/globaldefs.py', 'Physical constants and unit conversions'),
    ('src/fastkernels.py', 'Optional Numba kernels, with a NumPy fallback'),
    ('src/slabs.py', 'Propagation across adjacent slabs'),
    ('src/earth.py', 'PREM, chord geometry, and Earth crossings'),
    ('tests/', 'Regression suite, run with pytest'),
    ('tests/conftest.py', 'Shared fixtures and path setup'),
    ('tests/test_su3_algebra.py', 'd tensor, star product, SU(3) invariants'),
    ('tests/test_oscprob4nu.py', 'SU(4) algebra, quartic roots, 3+1 physics'),
    ('tests/test_evolution_operator.py',
     'U against an independent matrix exponential'),
    ('tests/test_probabilities.py', 'Normalization, positivity, P = |U|^2'),
    ('tests/test_hamiltonians.py', 'Sample Hamiltonians and sign conventions'),
    ('tests/test_reference_formulas.py',
     'Exact result against the standard formulas'),
    ('tests/test_matter_eigenvalues.py',
     'Matter spectrum, against Zaglauer-Schwarzer'),
    ('tests/test_edge_cases.py',
     'Degenerate and near-degenerate Hamiltonians'),
    ('tests/test_docstrings.py',
     'Runs the examples embedded in the docstrings'),
    ('tests/test_vectorized.py', 'The batched path, against the scalar one'),
    ('tests/test_vectorized_hamiltonians.py',
     'Hamiltonians built for an array of energies'),
    ('tests/test_annotations.py',
     'Annotations, and their agreement with the docs'),
    ('tests/test_fastkernels.py', 'Both backends, against each other'),
    ('tests/test_physical_scales.py',
     'Both backends at the scales actually used'),
    ('tests/test_slabs.py', 'Slab composition, against expm'),
    ('tests/test_earth.py', 'PREM, geometry, and Earth probabilities'),
    ('tests/test_batching_and_tolerance.py',
     'Energy batching, and the rtol/atol refinement'),
    ('tests/bit_capture.py',
     'Exact-bit capture, for refactors meant to be invisible'),
    ('tests/bit_compare.py', 'Compares two captures, in ulps'),
    ('tests/test_documented_figures.py',
     'Keeps the quoted performance figures agreeing'),
    ('tests/test_version_consistency.py',
     'Keeps the version agreeing wherever it is implied'),
    ('tests/test_nusquids_comparison.py',
     'Against nuSQuIDS, an independent external code'),
    ('tests/nusquids_reference.py',
     'Regenerates the frozen nuSQuIDS reference data'),
    ('tests/nusquids_reference.json',
     'Those reference values, with their provenance'),
    ('tests/nusquids_scan.py',
     'Regenerates the frozen nuSQuIDS energy scan'),
    ('tests/nusquids_scan.json', 'That scan, for the paper\'s figures'),
    ('tests/nufast_scan.json', 'NuFast-LBL at two Newton settings'),
    ('tests/speed_accuracy.json',
     'The six-code constant-density speed-accuracy plane'),
    ('tests/timing_other_codes.json', 'Timings behind the performance figure'),
    ('tests/prem_scan.py',
     'Regenerates the two Earth speed-accuracy planes'),
    ('tests/prem_speed_accuracy.json',
     'Those planes, at three flavors and at 3+1'),
    ('tests/external_drivers/',
     'Drivers for the codes that cannot be called from Python'),
    ('tests/external_drivers/README.md',
     'Every convention each one had to be told, and why'),
    ('tests/external_drivers/gen_prem_header.py',
     'Emits this library\'s PREM as a C header'),
    ('tests/external_drivers/nufast_drv.cpp', 'NuFast-Earth on the PREM chord'),
    ('tests/external_drivers/nufast_earth_prem.txt', None),
    ('tests/external_drivers/globes_drv.c', 'GLoBES on the PREM chord'),
    ('tests/external_drivers/globes_prem.txt', None),
    ('tests/external_drivers/prob3_drv.cpp', 'Prob3++ on the PREM chord'),
    ('tests/external_drivers/prob3_prem.txt', None),
    ('tests/test_file_tree.py', 'Keeps this tree in step with the repository'),

]


def render_tree():
    r"""Returns the file tree as a list of lines, root label included."""
    paths = [path for path, _ in TREE]
    comments = dict(TREE)

    def parent(path):
        trimmed = path.rstrip('/')
        return trimmed.rsplit('/', 1)[0] + '/' if '/' in trimmed else ''

    # An entry is the last of its siblings when no later entry shares its
    # parent.  The connectors and the vertical bars both follow from that.
    has_later_sibling = {}
    for index, path in enumerate(paths):
        has_later_sibling[path] = any(
            parent(other) == parent(path) for other in paths[index+1:])

    lines = [ROOT_LABEL]
    for path in paths:
        ancestors = []
        walk = parent(path)
        while walk:
            ancestors.append(walk)
            walk = parent(walk)
        prefix = ''.join('\u2502   ' if has_later_sibling[a] else '    '
                         for a in reversed(ancestors))
        connector = '\u2514\u2500\u2500 ' if not has_later_sibling[path] \
            else '\u251c\u2500\u2500 '
        line = prefix + connector + path.rstrip('/').rsplit('/', 1)[-1]
        if path.endswith('/'):
            line += '/'
        comment = comments[path]
        if comment is not None:
            if len(line) >= COMMENT_COLUMN:
                line += '  # ' + comment
            else:
                line = line.ljust(COMMENT_COLUMN) + '# ' + comment
        lines.append(line)
    return lines


def _readme_block(text):
    r"""Returns (start, end) line indices of the tree fence in README.md."""
    lines = text.split('\n')
    start = next(i for i, line in enumerate(lines) if line == '```text')
    end = next(i for i in range(start+1, len(lines)) if lines[i] == '```')
    return lines, start+1, end


def _install_block(text):
    r"""Returns (start, end) line indices of the tree block in the rst."""
    lines = text.split('\n')
    start = next(i for i, line in enumerate(lines)
                 if line.strip() == ROOT_LABEL)
    end = start
    while end < len(lines) and lines[end].strip():
        end += 1
    return lines, start, end


def current_readme_tree():
    r"""Returns the tree as it currently appears in ``README.md``."""
    with open(README) as handle:
        lines, start, end = _readme_block(handle.read())
    return lines[start:end]


def current_install_tree():
    r"""Returns the tree as it currently appears in the rst, undented."""
    with open(INSTALL_RST) as handle:
        lines, start, end = _install_block(handle.read())
    return [line[3:] if line.startswith('   ') else line
            for line in lines[start:end]]


def write_trees():
    r"""Rewrites both documents from `TREE`.  Returns the paths changed."""
    generated = render_tree()
    changed = []

    with open(README) as handle:
        text = handle.read()
    lines, start, end = _readme_block(text)
    if lines[start:end] != generated:
        lines[start:end] = generated
        with open(README, 'w') as handle:
            handle.write('\n'.join(lines))
        changed.append(README)

    with open(INSTALL_RST) as handle:
        text = handle.read()
    lines, start, end = _install_block(text)
    indented = ['   ' + line if line else line for line in generated]
    if lines[start:end] != indented:
        lines[start:end] = indented
        with open(INSTALL_RST, 'w') as handle:
            handle.write('\n'.join(lines))
        changed.append(INSTALL_RST)

    return changed


def test_tree_matches_git():
    r"""The file entries in `TREE` are exactly the tracked files."""
    tracked = set(subprocess.check_output(
        ['git', 'ls-files'], cwd=ROOT).decode().split())
    listed = {path for path, _ in TREE if not path.endswith('/')}

    missing = sorted(tracked - listed)
    stale = sorted(listed - tracked)
    assert not missing, (
        'tracked but absent from TREE in tests/test_file_tree.py: %s'
        % ', '.join(missing))
    assert not stale, (
        'listed in TREE but not tracked: %s' % ', '.join(stale))


def test_every_directory_in_tree_has_a_parent_entry():
    r"""Each nested entry's parent is listed, and listed before it."""
    seen = set()
    for path, _ in TREE:
        trimmed = path.rstrip('/')
        if '/' in trimmed:
            parent = trimmed.rsplit('/', 1)[0] + '/'
            assert parent in seen, (
                '%s appears before its parent %s' % (path, parent))
        seen.add(path)


def test_readme_tree_is_generated():
    r"""README.md carries exactly what `render_tree` produces."""
    assert current_readme_tree() == render_tree(), (
        'the file tree in README.md is out of date; regenerate it with '
        '`python tests/test_file_tree.py --write`')


def test_installation_rst_tree_is_generated():
    r"""installation.rst carries the same tree, indented for rst."""
    assert current_install_tree() == render_tree(), (
        'the file tree in docs/source/installation.rst is out of date; '
        'regenerate it with `python tests/test_file_tree.py --write`')


def test_tree_is_not_empty():
    r"""The generator produced a tree, rather than nothing."""
    assert len(render_tree()) > 20


def test_every_tree_entry_exists_on_disk():
    r"""The tree does not list files that are not there."""
    for path, _ in TREE:
        assert os.path.exists(os.path.join(ROOT, path.rstrip('/'))), (
            'listed in TREE but absent from disk: %s' % path)


if __name__ == '__main__':
    if '--write' in sys.argv:
        updated = write_trees()
        if updated:
            for name in updated:
                print('rewrote %s' % os.path.relpath(name, ROOT))
        else:
            print('both documents were already up to date')
    else:
        stale_docs = []
        if current_readme_tree() != render_tree():
            stale_docs.append('README.md')
        if current_install_tree() != render_tree():
            stale_docs.append('docs/source/installation.rst')
        if stale_docs:
            print('out of date: %s' % ', '.join(stale_docs))
            print('run `python tests/test_file_tree.py --write` to update')
            sys.exit(1)
        print('both documents are up to date')
