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

A single three-flavor probability evaluation takes about 300 microseconds.
The :math:`d` tensor is constant and is tabulated once at import time as a
dense :math:`8\times8\times8` array; the star product is a baseline-
independent contraction, computed once per call.

For scans over energy or baseline, note that the *Hamiltonian* is the only
thing that changes with energy, and the latent roots depend on it alone ---
so the energy-independent vacuum term should be built once, outside the loop,
as the bundled examples do.  Vectorising the whole expansion over a baseline
axis is a further factor of a few hundred, and would require an interface
that accepts arrays; it is not implemented.
