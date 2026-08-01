Methodology
===========

This page sketches how the closed-form probabilities are obtained, and
records the conventions the code uses.  For the derivation, see
:cite:`Bustamante:2019ggq` and the original treatment in
:cite:`Ohlsson:1999xb`.

The idea
--------

For a time-independent Hamiltonian, the evolution operator is the matrix
exponential

.. math:: U(L) = e^{-i H L} .

Rather than evaluating that exponential numerically, expand :math:`H` in the
basis of generators of SU(2) or SU(3) and exponentiate *analytically*, using
the closure of the algebra.  The result is a closed form for :math:`U(L)` in
the same basis, whose coefficients are elementary functions of two
invariants of the Hamiltonian.

Two flavors
-----------

Expand the Hamiltonian in the Pauli matrices,

.. math:: H = h_0 \mathbb{1} + h_k \sigma^k ,

with :math:`h_k` real for Hermitian :math:`H`.  The term :math:`h_0`
multiplies the identity, so it contributes only an overall phase and is
dropped; the code works with the traceless part throughout, which leaves the
probabilities unchanged.

Exponentiating gives

.. math::
   U_2(L) = u_0 \mathbb{1} + i u_k \sigma^k , \qquad
   u_0 = \cos(|h| L) , \qquad
   u_k = -\frac{h_k}{|h|}\sin(|h| L) ,

with :math:`|h|^2 = h_1^2 + h_2^2 + h_3^2`.  The transition probability
follows directly:

.. math::
   P_{e\mu} = |[U_2]_{\mu e}|^2 = u_1^2 + u_2^2
            = \frac{h_1^2 + h_2^2}{|h|^2}\sin^2(|h| L) .

Both :math:`h_1` and :math:`h_2` appear.  :math:`h_2 = -\mathrm{Im}(H_{12})`
vanishes only when the off-diagonal entry is real, which is why the
:math:`h_2` term is easy to lose and hard to notice: it makes no difference
in vacuum or in matter of constant density, and all the difference for a
CP-violating Hamiltonian.

When :math:`|h| = 0` the Hamiltonian is proportional to the identity, there
is no flavor evolution, and the code takes the limit
:math:`\sin(|h|L)/|h| \to L` explicitly.

Three flavors
-------------

Expand in the Gell-Mann matrices,

.. math:: H = h_0 \mathbb{1} + h_k \lambda^k , \qquad k = 1, \ldots, 8 ,

again dropping :math:`h_0`.  SU(3) does not close as simply as SU(2): the
anticommutator of two generators brings in the totally symmetric tensor

.. math::
   d_{ijk} = \tfrac{1}{4}\,\mathrm{Tr}
             \left(\{\lambda_i, \lambda_j\}\lambda_k\right) ,

whose conventions follow :cite:`MacFarlane:1968`.  Two invariants control the
result: the norm :math:`|h|^2 = h_i h_i` and the cubic invariant

.. math:: \langle h \rangle = d_{ijk}\, h_i h_j h_k = h_i (h \star h)_i ,

where :math:`(h \star h)_i = d_{ijk} h_j h_k` is the *star product*.  These
equal :math:`\mathrm{Tr}(H_0^2)/2` and :math:`\mathrm{Tr}(H_0^3)/2`
respectively, with :math:`H_0` the traceless part --- a relation the test
suite checks.

The characteristic equation of :math:`-H_0`,

.. math:: \psi^3 - |h|^2 \psi - \tfrac{2}{3}\langle h \rangle = 0 ,

has three real roots :math:`\psi_m`, given in closed trigonometric form.
The coefficients of

.. math:: U_3(L) = u_0 \mathbb{1} + i u_k \lambda^k

then follow by Lagrange interpolation over those roots:

.. math::
   u_0 = \frac{1}{3}\sum_m e^{i \psi_m L} , \qquad
   u_k = i \sum_m e^{i \psi_m L}\,
         \frac{\psi_m h_k - (h \star h)_k}{3\psi_m^2 - |h|^2} .

Degenerate spectra
^^^^^^^^^^^^^^^^^^

The denominator :math:`3\psi_m^2 - |h|^2` is the derivative of the
characteristic polynomial, so it vanishes exactly at a repeated root.  Two
degenerate cases are handled separately, and exactly:

* :math:`|h|^2 = 0`, where the Hamiltonian is proportional to the identity
  and :math:`U_3 = \mathbb{1}`;
* a doubly degenerate root :math:`\psi_a = \psi_b \neq \psi_c`, where the
  spectral decomposition collapses onto a single projector,

  .. math::
     U_3 = e^{i \psi_a L}\mathbb{1}
         + \left(e^{i \psi_c L} - e^{i \psi_a L}\right) P_c , \qquad
     P_c = \frac{h_k \lambda^k + \psi_a \mathbb{1}}{\psi_a - \psi_c} ,

  which is linear in :math:`h_k` and so needs no matrix algebra.

The argument of the arc cosine that produces the roots lies in
:math:`[-1, 1]` for any Hermitian Hamiltonian, but only up to round-off; it
is clipped, so that a marginally out-of-range value cannot yield complex
roots and a non-unitary evolution operator.

.. _sign-convention:

Sign conventions
----------------

The vacuum Hamiltonians are built so that **a positive matter potential
added to the** :math:`ee` **entry describes neutrinos**:

.. math::
   H^{2\nu}_{\rm vac} = \frac{\Delta m^2}{4E}
   \begin{pmatrix} -\cos 2\theta & \sin 2\theta \\
                    \sin 2\theta & \cos 2\theta \end{pmatrix} ,
   \qquad
   H^{3\nu}_{\rm vac} = \frac{1}{2E}\, U M^2 U^\dagger ,

with :math:`M^2 = \mathrm{diag}(0, \Delta m^2_{21}, \Delta m^2_{31})` and
:math:`U` the PMNS matrix in the standard PDG parametrization.

This is worth stating explicitly because an overall sign flip of the vacuum
term alone is **invisible in vacuum** --- for a real Hamiltonian the
probabilities are invariant under :math:`H \to -H` --- and yet reverses the
sign of the matter potential *relative* to it, turning neutrinos into
antineutrinos and moving the Mikheyev-Smirnov-Wolfenstein resonance to the
other side.  For antineutrinos, pass a negative ``VCC`` (and conjugate the
CP phase).

Ordering of the probabilities
-----------------------------

:func:`oscprob3nu.probabilities_3nu` returns the nine probabilities with the
*initial* flavor varying slowest, so the returned tuple reads

.. math::
   (P_{ee},\, P_{e\mu},\, P_{e\tau},\,
    P_{\mu e},\, P_{\mu\mu},\, P_{\mu\tau},\,
    P_{\tau e},\, P_{\tau\mu},\, P_{\tau\tau})

with :math:`P_{\alpha\beta} \equiv P(\nu_\alpha \to \nu_\beta)
= |[U_3]_{\beta\alpha}|^2`.  Note the index order: the evolution operator is
indexed (final, initial), the probabilities (initial, final).

Cost
----

The :math:`d` tensor is constant and is tabulated once at import time as a
dense :math:`8\times8\times8` array; the star product and the two invariants
are contractions against that table.

A single three-flavor probability evaluation takes about eight
microseconds, and a two-flavor one about one.  For scans, pass arrays rather
than looping: the routines
accept a stack of Hamiltonians, a stack of baselines, or both, and evaluate
the stack in one pass.  Measured against the equivalent Python loop, on 2000
points:

.. list-table::
   :header-rows: 1
   :widths: 46 27 27

   * - Scan
     - Speedup
     - Also comparable to
   * - Versus baseline (one :math:`H`, many :math:`L`)
     - ~30x
     - one ``eigh`` plus phases
   * - Versus energy (many :math:`H`, one :math:`L`)
     - ~25x
     - batched ``numpy.linalg.eigh``
   * - Oscillogram, 100 x 100
     - ~40x
     -
   * - Two flavors, versus baseline
     - ~70x
     -

These ratios have *narrowed* across successive releases even as both sides got
quicker: the scalar path has itself sped up several-fold, so the loop being
compared against is no longer as slow as it was.  In absolute terms the
vectorised scan is faster than it has ever been.

The two scans differ because the latent roots depend on the Hamiltonian
alone.  Scanning one Hamiltonian over many baselines solves the
characteristic equation once and then only evaluates
:math:`e^{i\psi_m L}`, so almost all the work is amortised; scanning over
energy changes the Hamiltonian at every point and must solve it each time.

The vectorised path is as fast as diagonalising with LAPACK, which is the
honest comparison to draw: the SU(3) route's advantage is that it is a
closed form, not that it outruns ``eigh``.  What the vectorisation removes
is the *disadvantage* it used to carry.

Short stacks
^^^^^^^^^^^^

A batched call carries a couple of hundred microseconds of fixed cost
whatever its length --- allocating and reducing a dozen small arrays --- so
for a handful of points it spends more on the machinery than the scalar path
spends on the whole job.  Stacks below
:data:`oscprob3nu.SMALL_BATCH` are therefore evaluated one element at a time.

The thresholds are measured, not guessed, and differ between the two
expansions because the two-flavor one does much less work per element and so
amortises its overhead sooner: eleven elements for three flavors, seven for
two.  Nothing about this is visible from the outside; the answers are the
same either way.

The optional compiled backend
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With :mod:`fastkernels` --- that is, with Numba installed --- the batched
paths are compiled instead.  The NumPy path evaluates the expansion as a
chain of whole-array operations, so a stack of :math:`N` Hamiltonians makes
roughly fifteen passes over :math:`N`-element arrays, each writing a
temporary that the next reads back.  The compiled kernel does the same
arithmetic one element at a time, keeping every intermediate in registers,
and spreads the elements over the available cores.  Against the NumPy path:

.. list-table::
   :header-rows: 1
   :widths: 55 45

   * - Stack
     - Speedup
   * - 200 000 energies, three flavors
     - ~15x
   * - 20 000 energies, three flavors
     - ~13x
   * - 2 000 energies, three flavors
     - ~3x
   * - 200 000 baselines, two flavors
     - ~1.4x
   * - 2 000 baselines, two flavors
     - NumPy is quicker

The two-flavor rows are the honest caveat: that path reduces to a square root
and a sine per element, which NumPy already does about as well as compiled
code can, and the kernel additionally has to materialise the Hamiltonian
stack --- for a scan over baselines, the same matrix repeated, which costs
2.5 ms to copy at two hundred thousand points.

So the backend is not used unconditionally.  :func:`fastkernels.worthwhile`
declines it below the per-flavor thresholds in
:data:`fastkernels.MIN_BATCH`, which were found by alternating the two paths
and taking the best of nine rounds each: for three flavors the kernel wins at
every size, for two it wins from about fifty thousand elements.  A backend
that is sometimes slower than the path it replaces is worse than no backend,
and a test asserts the thresholds are honoured.

The scalar path is deliberately left uncompiled.  One probability takes
about eight microseconds; compiling it would save most of that, at the cost of
a multi-second pause on a user's first call.

Degenerate spectra on the batched path
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The degenerate branch cannot be taken elementwise inside a vectorised
expression.  The general formula is therefore evaluated everywhere, with
vanishing denominators replaced by one, and the affected elements are then
recomputed individually.  Degeneracy is measure-zero among floating-point
Hamiltonians, so that fallback loop is empty in essentially every real use.
