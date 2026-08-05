"""Exact-bit capture of the Earth and slab paths, for refactoring.

**This is not a test and pytest does not collect it.**  It is the
harness for a change that is supposed to leave every number alone ---
hoisting a loop-invariant, removing a temporary, reordering a memory
access --- where "supposed to" needs checking rather than asserting.

Used in pairs, with `bit_compare.py`::

    python tests/bit_capture.py src before.npz
    # ...make the change...
    python tests/bit_capture.py src after.npz
    python tests/bit_compare.py before.npz after.npz

No golden file is committed alongside these, deliberately.  ``acos``,
``cos`` and ``sin`` come from the platform's libm and are not guaranteed
bit-identical between implementations or versions, so a stored baseline
would fail spuriously somewhere in CI while saying nothing about the
change under test.  Comparing two captures taken on the same machine in
the same session removes that variable entirely.

Captures evolution *operators* as well as probabilities: |U|^2 discards
information, so a change that perturbs a phase but not a modulus would
slip past a probabilities-only comparison.  Measured on a deliberately
mutated kernel, the operator arrays showed 17000 ulps where the
probabilities showed 1800.

Both backends, symmetric and asymmetric slab sequences, odd and even
slab counts, and the scalar, array and grid entry points.

Captures evolution *operators* as well as probabilities: |U|^2 discards
information, so a change that perturbs a phase but not a modulus would
slip past a probabilities-only comparison.

Both backends, symmetric and asymmetric slab sequences, odd and even slab
counts, and the scalar, array and grid entry points -- Group A edits
`_build_h_*` (fused path only) and the tmp->acc copy (which appears in
the per-chord, batched, fused and mirrored composers alike), so all of
them have to be pinned.
"""
import sys

import numpy as np

sys.path.insert(0, sys.argv[1])

import earth                                             # noqa: E402
import fastkernels as fk                                 # noqa: E402
import globaldefs as gd                                  # noqa: E402
import hamiltonians2nu as h2                             # noqa: E402
import hamiltonians3nu as h3                             # noqa: E402
import hamiltonians4nu as h4                             # noqa: E402
import slabs                                             # noqa: E402

HV = {
    2: h2.hamiltonian_2nu_vacuum_energy_independent(gd.S12_NO_BF,
                                                    gd.D21_NO_BF),
    3: h3.hamiltonian_3nu_vacuum_energy_independent(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, gd.DCP_NO_BF,
        gd.D21_NO_BF, gd.D31_NO_BF),
    4: h4.hamiltonian_4nu_vacuum_energy_independent(
        gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF, 0.1, 0.1, 0.1,
        gd.D21_NO_BF, gd.D31_NO_BF, 1.0, 0.0, 0.0, 0.0),
}
EARTH = {2: earth.probabilities_2nu_earth, 3: earth.probabilities_3nu_earth,
         4: earth.probabilities_4nu_earth}
SLAB_U = {2: slabs.evolution_operator_2nu_slabs,
          3: slabs.evolution_operator_3nu_slabs,
          4: slabs.evolution_operator_4nu_slabs}
SLAB_P = {2: slabs.probabilities_2nu_slabs,
          3: slabs.probabilities_3nu_slabs,
          4: slabs.probabilities_4nu_slabs}

# Angles chosen to give both odd and even slab counts, short and long
# chords, and a diametric crossing.
ANGLES = (-0.05, -0.2, -0.55, -0.8, -0.95, -1.0)
# 1 and 3 give odd slab counts per segment; 2, 8, 21 exercise the rest.
SUBDIVISIONS = (1, 2, 3, 8, 21)
ENERGIES = np.array([1.0e9, 3.7e9, 1.0e10, 4.0e10, 2.5e11])

captured = {}


def keep(name, value):
    r"""Stores an array under a name, refusing a silent overwrite."""
    assert name not in captured, name
    captured[name] = np.ascontiguousarray(value)


for backend in ('numba', 'numpy'):
    fk.USE_NUMBA = (backend == 'numba')

    for n_flavors in (2, 3, 4):
        hv = HV[n_flavors]
        fn = EARTH[n_flavors]

        for costhz in ANGLES:
            for n_per in SUBDIVISIONS:
                tag = '%s_%dnu_cz%s_n%d' % (backend, n_flavors, costhz, n_per)

                # Scalar entry point, one row per energy
                keep('earth_scalar_' + tag,
                     np.array([np.asarray(fn(hv, e, costhz, n_per),
                                          dtype=float)
                               for e in ENERGIES]))

                # Array entry point: the fused/batched path
                keep('earth_array_' + tag,
                     np.asarray(fn(hv, ENERGIES, costhz, n_per), dtype=float))

                # The operators themselves, which probabilities discard
                h, widths = earth._earth_hamiltonians(
                    hv, ENERGIES[2], costhz, n_per,
                    gd.ELECTRON_FRACTION_EARTH_CRUST, n_flavors)
                keep('slab_u_' + tag,
                     np.asarray(SLAB_U[n_flavors](h, widths), dtype=complex))
                keep('slab_p_' + tag,
                     np.asarray(SLAB_P[n_flavors](h, widths), dtype=float))

        # A grid, which groups angles and reuses geometry
        keep('earth_grid_%s_%dnu' % (backend, n_flavors),
             np.asarray(fn(hv, ENERGIES[None, :],
                           np.array(ANGLES)[:, None]), dtype=float))

        # Asymmetric sequences, so the non-palindromic composer is pinned
        # too, at several lengths including odd ones.
        for n_slabs in (1, 2, 5, 16, 17, 40):
            rng = np.random.default_rng(4200 + n_slabs + 17*n_flavors)
            a = (rng.normal(size=(n_slabs, n_flavors, n_flavors))
                 + 1.0j*rng.normal(size=(n_slabs, n_flavors, n_flavors)))
            h = (a + np.conj(np.swapaxes(a, -1, -2)))/2.0
            w = rng.uniform(0.1, 2.0, size=n_slabs)
            tag = '%s_%dnu_rand%d' % (backend, n_flavors, n_slabs)
            keep('rand_u_' + tag,
                 np.asarray(SLAB_U[n_flavors](h, w), dtype=complex))
            keep('rand_p_' + tag,
                 np.asarray(SLAB_P[n_flavors](h, w), dtype=float))

            # ...and its palindrome, for the mirrored composer
            hs = np.concatenate([h, h[::-1]])
            ws = np.concatenate([w, w[::-1]])
            keep('pal_u_' + tag,
                 np.asarray(SLAB_U[n_flavors](hs, ws), dtype=complex))
            keep('pal_p_' + tag,
                 np.asarray(SLAB_P[n_flavors](hs, ws), dtype=float))

fk.USE_NUMBA = True

np.savez(sys.argv[2], **captured)
total = sum(v.size for v in captured.values())
print('captured %d arrays, %d values, from %s'
      % (len(captured), total, sys.argv[1]))
