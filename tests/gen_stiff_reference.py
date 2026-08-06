# -*- coding: utf-8 -*-
r"""Generates the four-flavor accuracy reference, at fifty decimal digits.

Run from the repository root, after installing ``mpmath``::

    python tests/gen_stiff_reference.py

It writes ``tests/stiff_reference.json``, which `test_oscprob4nu.py` reads.
``mpmath`` is needed only here, so it is not a test dependency.

Why this file exists
--------------------

The four-flavor latent roots used to be checked against
:func:`numpy.linalg.eigvalsh`.  That was an independent oracle for as long as
the code solved the characteristic quartic in closed form --- and stopped
being one the moment the eigensolver became an implementation option, at
which point those tests compared the implementation against itself.  One of
them duly began asserting ``1.06e-15 <= 0.0``, the zero being eigvalsh minus
eigvalsh.

So the oracle is frozen here instead, computed at eighty digits by a route
the library shares no code with: :func:`mpmath.eig` for the roots and
:func:`mpmath.expm` for the probabilities, with `residual` below checking the
roots against :math:`\det(\psi\mathbb{1} - \tilde{H})` taken straight from
the matrix.

Three traps are worth recording, because each cost a day
-------------------------------------------------------

**The quartic built from the invariants is not always the matrix's
characteristic polynomial.**  Forming the traceless part in float64 leaves a
residual trace :math:`\tau` --- small, but not always zero --- while that
quartic pins its cubic coefficient to exactly zero and carries invariants of
a matrix whose trace is not quite zero.  The two formulations then disagree,
by up to 3.8e-07 relative here, and the inconsistent one can even develop a
complex-conjugate pair that a real iteration will orbit forever.  Wherever
:math:`\tau` is exactly zero they agree to 2e-17, which is what makes this so
easy to miss.  An earlier version of this file used the quartic and therefore
shipped a reference that was wrong on three of five stiff cases.

**Two solvers agreeing is not evidence.**  That wrong reference was
corroborated by ``mpmath.eig``, and the agreement was taken as proof --- so
:func:`numpy.linalg.eigvalsh` was blamed for residuals of 2.9e-9 and 1.8e-5
that were the oracle's, not its.  Against the determinant it is accurate to
2--7e-16 throughout.  Only a check that depends on no eigenvalue algorithm at
all settles this, which is why `residual` exists and why it uses ``det``.

**Probability conservation is not an accuracy measure.**  A row of
:math:`|U|^2` sums to one for any unitary :math:`U`, however wrong its
eigenvalues, so a spectral reconstruction conserves probability perfectly
while being badly inaccurate.  The probabilities are therefore stored
outright and compared entry by entry.

Matrices are stored as hexadecimal floats, so a reader of the JSON gets the
exact bits the reference was computed from rather than a decimal rounding of
them, and the cases do not depend on `globaldefs` staying put.
"""

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), 'src'))

import globaldefs as gd            # noqa: E402
import hamiltonians4nu             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'stiff_reference.json')
DIGITS = 50


def traceless(matrix):
    r"""Returns the traceless part, as the expansions use it."""
    matrix = np.asarray(matrix, dtype=complex)
    return matrix - np.trace(matrix)/4.0*np.eye(4)


def cases():
    r"""Returns the named Hamiltonians the reference covers.

    The 3+1 family spans the stiffness that broke the closed form: 0.1 and
    1 eV^2 are the physical range, 100 upwards is where the invariants stop
    holding the cluster apart.  The random and near-degenerate families are
    the ones a change to the root strategy is most likely to regress
    without anybody noticing.
    """
    out = []
    for dm41 in (0.1, 1.0, 10.0, 100.0, 1000.0):
        h = hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
            gd.S12_NO_BF, gd.S23_NO_BF, gd.S13_NO_BF,
            np.sqrt(0.10), np.sqrt(0.10), 0.0,
            gd.DCP_NO_BF, gd.D21_NO_BF, gd.D31_NO_BF, dm41)
        out.append(('3+1, dm41^2 = %g eV^2' % dm41,
                    np.asarray(h, dtype=complex)/1.0e9,
                    1300.0*gd.CONV_KM_TO_INV_EV))

    rng = np.random.default_rng(20260806)
    for k in range(2):
        a = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
        out.append(('random Hermitian %d' % k,
                    ((a + a.conj().T)/2.0)*1.0e-13,
                    1300.0*gd.CONV_KM_TO_INV_EV))

    for separation in (1.0e-8, 1.0e-16):
        w = np.array([1.0, 1.0 + separation, -0.7, -1.3])*1.0e-13
        a = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
        q, _ = np.linalg.qr(a)
        out.append(('pair separated by %.0e' % separation,
                    q @ np.diag(w).astype(complex) @ q.conj().T,
                    1300.0*gd.CONV_KM_TO_INV_EV))
    return out


def reference(matrix, baseline):
    r"""Returns the exact roots and probabilities, at `DIGITS` digits."""
    import mpmath as mp
    mp.mp.dps = DIGITS + 30

    h0 = traceless(matrix)
    m = mp.matrix(4, 4)
    for i in range(4):
        for j in range(4):
            m[i, j] = mp.mpc(h0[i, j].real, h0[i, j].imag)

    # The eigenvalues of the matrix itself.  NOT the roots of the quartic
    # built from I_2, I_3, I_4: float64 traceless-ing leaves a residual
    # trace tau, and that quartic pins its cubic coefficient to exactly zero
    # while carrying invariants of a matrix whose trace is not quite zero.
    # The two disagree by up to 3.8e-7 relative wherever tau != 0 -- which
    # is a defect of the polynomial, not of the matrix, and cost a day of
    # blaming `numpy.linalg.eigvalsh` for an error that was in the oracle.
    # `residual` below checks these against det(psi*I - H~) taken straight
    # from the matrix, which is the one formulation with nothing to be
    # inconsistent with.
    roots = sorted(mp.re(r) for r in mp.eig(m, left=False, right=False))

    # P_ab = |U_ba|^2, the library's convention.  Getting this transposed
    # produces a plausible 3e-2 disagreement that looks like a real defect.
    u = mp.expm(m*mp.mpc(0, -1)*mp.mpf(baseline))
    probabilities = [float(mp.fabs(u[b, a])**2)
                     for a in range(4) for b in range(4)]

    return [float(r) for r in roots], probabilities


def residual(matrix, roots):
    r"""Returns the worst normalised residual of `roots`, as a self-check.

    A candidate root should annihilate the characteristic polynomial.
    Dividing :math:`\chi(\psi_m)` by :math:`\prod_{l \neq m}(\psi_m -
    \psi_l)` --- which is :math:`\chi'(\psi_m)` --- converts the residual
    into the distance to the true root, so the number below is directly
    comparable to a relative root error.

    The determinant is taken from the matrix rather than from the
    invariants, which makes this the only check here that depends on no
    eigenvalue algorithm and on no reformulation of the polynomial.  That is
    what settled which oracle was wrong: the roots of the invariants-quartic
    scored 2.9e-9 at 10 eV\ :sup:`2` and 1.8e-5 at 1000 against it, while
    :func:`numpy.linalg.eigvalsh` scored 2--7e-16 throughout.  Two solvers
    agreeing with each other is not evidence; annihilating the matrix's own
    characteristic polynomial is.
    """
    import mpmath as mp
    mp.mp.dps = DIGITS + 30

    h0 = traceless(matrix)
    m = mp.matrix(4, 4)
    for i in range(4):
        for j in range(4):
            m[i, j] = mp.mpc(h0[i, j].real, h0[i, j].imag)

    def chi(value):
        shifted = mp.matrix(4, 4)
        for i in range(4):
            for j in range(4):
                shifted[i, j] = ((mp.mpf(value) if i == j else mp.mpf(0))
                                 - m[i, j])
        return mp.det(shifted)

    scale = max(abs(r) for r in roots)
    worst = 0.0
    for index, root in enumerate(roots):
        derivative = mp.mpf(1)
        for other in (r for k, r in enumerate(roots) if k != index):
            derivative *= mp.mpf(root) - mp.mpf(other)
        if derivative == 0:
            continue
        worst = max(worst, float(abs(chi(root)/derivative))/scale)
    return worst


def main():
    payload = {
        'digits': DIGITS,
        'note': ('Generated by tests/gen_stiff_reference.py.  Roots from '
                 'mpmath.eig on the traceless matrix, verified against '
                 'det(psi*I - H~); probabilities from mpmath.expm.  Neither '
                 'shares code with the library, and neither goes through '
                 'the invariants.'),
        'cases': [],
    }
    for label, matrix, baseline in cases():
        roots, probabilities = reference(matrix, baseline)
        payload['cases'].append({
            'label': label,
            'baseline': baseline.hex(),
            'matrix_real': [[matrix[i, j].real.hex() for j in range(4)]
                            for i in range(4)],
            'matrix_imag': [[matrix[i, j].imag.hex() for j in range(4)]
                            for i in range(4)],
            'roots': [r.hex() for r in roots],
            'probabilities': [p.hex() for p in probabilities],
        })
        check = residual(matrix, roots)
        assert check < 1.0e-15, (
            '%s: the stored roots do not annihilate the characteristic '
            'polynomial, residual %.2e; the reference is wrong, not the '
            'library' % (label, check))
        print('  %-28s roots %.3e .. %.3e   residual %.1e'
              % (label, min(abs(r) for r in roots),
                 max(abs(r) for r in roots), check))

    with open(OUT, 'w') as handle:
        json.dump(payload, handle, indent=1, sort_keys=True)
        handle.write('\n')
    print('wrote %s (%d cases)' % (os.path.relpath(OUT), len(payload['cases'])))


if __name__ == '__main__':
    main()
