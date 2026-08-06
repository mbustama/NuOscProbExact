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

Four flavors
------------

Everything above generalizes with :math:`3 \to 4`.  Expand in the fifteen
generalized Gell-Mann matrices,

.. math:: H = h_0 \mathbb{1} + h_a \lambda^a , \qquad a = 1, \ldots, 15 ,

again dropping :math:`h_0`.  Three things are new, and each of them is a
consequence of SU(4) having rank three where SU(3) has rank two.

**A third invariant.**  The traceless part carries

.. math::
   I_2 = \tfrac12 \mathrm{Tr}\,\tilde{H}^2 , \qquad
   I_3 = \tfrac12 \mathrm{Tr}\,\tilde{H}^3 , \qquad
   I_4 = \tfrac12 \left(\mathrm{Tr}\,\tilde{H}^4 - I_2^2\right) ,

the first two being the :math:`|h|^2` and :math:`\langle h \rangle` of the
three-flavor case.  Taking them from traces means the SU(4) :math:`d` tensor
--- a :math:`15\times15\times15` table --- is never built.

**A quartic, which still solves.**  The characteristic equation becomes

.. math::
   \psi^4 - I_2 \psi^2 - \tfrac23 I_3 \psi
   + \tfrac14\left(I_2^2 - 2 I_4\right) = 0 ,

and Euler's reduction turns it into the *resolvent cubic*

.. math:: z^3 - 2 I_2 z^2 + 2 I_4 z - \tfrac49 I_3^2 = 0 ,

whose roots are :math:`z_i = (\psi_i + \psi_j)^2` --- real and non-negative
precisely because :math:`\tilde{H}` is Hermitian.  So the same trigonometric
formula used at three flavors solves it, and then

.. math::
   \psi_m = \tfrac12\left(s_1\sqrt{z_1} + s_2\sqrt{z_2} + s_3\sqrt{z_3}\right),
   \qquad s_1 s_2 s_3 \sqrt{z_1 z_2 z_3} = \tfrac23 I_3 .

The SU(3) machinery is literally nested inside the SU(4) solution.

**A longer star-product tower.**  The three-flavor identity
:math:`(h \star h) \star h = \tfrac13 |h|^2 h` is a Cayley-Hamilton accident
of :math:`n = 3` and is *false* at :math:`n = 4` --- over two thousand random
Hamiltonians the two sides differ by a median of 56%, with the central 90% of
draws between 34% and 116% --- so the third rung enters as independent data.
The median is the stable statistic here: the extremes of a finite draw are
not, and grow with the sample, which is why a range of largest-and-smallest
is not quoted.

.. math::
   u_0 = \frac14 \sum_m e^{-i\psi_m L} , \qquad
   i u_a = \sum_m e^{-i\psi_m L}\,
   \frac{\left(\psi_m^2 - \tfrac12 I_2\right) h_a
         + \psi_m (h \star h)_a + ((h \star h) \star h)_a}{\chi'(\psi_m)} .

That also exposes the general-:math:`n` pattern: a numerator of degree
:math:`n-2` in :math:`\psi_m`, a star tower cut off at length :math:`n-1`,
and always :math:`\chi'(\psi_m)` underneath.

.. _why-four-is-the-end:

Why the method stops at four
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Not for want of effort, and not because :math:`n = 5` is uninteresting.  The
whole construction rests on one thing: that the eigenvalues of the traceless
Hamiltonian --- the roots of its characteristic polynomial --- can be written
down *in radicals*, as an explicit formula in the invariants.

That polynomial has degree :math:`n`.  Quadratics, cubics and quartics are
solvable in radicals; the Abel-Ruffini theorem says the general quintic is
not, and Galois theory says why: the symmetric group :math:`S_5` is not
soluble, while :math:`S_2`, :math:`S_3` and :math:`S_4` are.  At
:math:`n = 5` there is no formula to write, and the shortfall is a theorem
rather than a gap in anyone's algebra.

So the closed-form road ends at four, and it ends for a reason external to
neutrino physics entirely.

What does *not* end there is the philosophy.  Nothing above the eigenvalues
needs radicals: the interpolation over the roots, the fact that no
eigenvectors are ever required, and the whole probability construction go
through for any :math:`n`.  Feed numerically computed eigenvalues into the
same Sylvester sum and the method degrades gracefully rather than breaking
--- which is what a general-SU(:math:`n`) treatment would do, and what codes
carrying SU(:math:`N`) expansions to :math:`n = 6` in fact do.  It would no
longer be a *closed form*, which is this library's reason to exist, so it is
out of scope here.

Four flavors is therefore both the natural stopping point and a useful one:
it is exactly what 3+1 sterile scenarios need.

.. _stiff-spectra:

Stiff spectra, and what they cost
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, the perspective, because the numbers below are small enough to be
misread as a problem.

**None of this is near any measurable effect.**  Oscillation probabilities
are confronted with data at the per-cent level at best, and the systematic
uncertainties of a real experiment dominate long before the fourth decimal
place.  Both figures in this paragraph are errors on a **probability**.  The
table further down measures something else --- the relative error of the
latent roots themselves --- so its numbers are smaller and are not
comparable to these: an error in a root reaches the probability as a phase
error :math:`\delta\psi\,L`, amplified by the baseline.  Even the worst
probability here --- the unrefined four-flavor result at
:math:`5\times10^{-7}` --- sits four or five orders of magnitude below
anything an experiment can resolve, and the refined one at :math:`10^{-9}`
is far beyond any physics requirement.

So why care?  Three reasons, none of them about a single probability:

* **The claim.**  This library says it computes probabilities exactly, with
  no approximation beyond round-off.  A probability wrong by
  :math:`5\times10^{-7}` is still round-off-limited in a sense, but it is not
  the same claim, and the difference should be stated rather than glossed.
* **Composition.**  :mod:`slabs` and :mod:`earth` multiply evolution
  operators across many layers, so a per-layer error accumulates.  What is
  invisible in one probability need not stay invisible across a hundred.
* **Regression testing.**  A suite that pins agreement at :math:`10^{-9}`
  catches a real mistake; one that pins it at :math:`10^{-6}` has room for a
  bug to hide in.

Now the mechanism.  A 3+1 Hamiltonian with :math:`\Delta m^2_{41} \sim 1`
eV\ :sup:`2` has a *stiff* spectrum: the eigenvalues span four orders of
magnitude, with three of them clustered.  The invariants
:math:`I_2, I_3, I_4` are sums over that spectrum, so forming them in double
precision compresses a :math:`4\times4` matrix into three numbers and loses
what separates the cluster.  Perturbing the three invariants at the
:math:`10^{-16}` level --- their own rounding --- moves the roots by
:math:`6\times10^{-11}` relative at :math:`\Delta m^2_{41} \sim 1`
eV\ :sup:`2`, and by more as the splitting grows: the amplification from
coefficients to roots measures :math:`2.3\times10^9`, which on the stiffest
Hamiltonian in the reference puts the closed form at
:math:`2.2\times10^{-7}`.  That is a property of the problem, not of the
solver: it is the classic ill-conditioning of polynomial roots with respect
to their coefficients, and it means no better root-finder helps.  Deflating
the quartic to a cubic first was tried, and does not.

There are two ways out of a bottleneck like that, and the library has
shipped both.

**Stop asking the invariants.**  One Newton step on

.. math:: \chi(\psi) = \det\left(\psi \mathbb{1} - \tilde{H}\right)

refines the roots using the Hamiltonian *entries* at full precision, which
never pass through the three-number bottleneck.  A second step changes
nothing --- Newton doubles the correct digits, and one already reaches that
route's floor --- so exactly one is taken.  This was the whole answer in
1.12.0, and :data:`oscprob4nu.POLISH_ROOTS` still switches it.

**Or widen the bottleneck.**  The compression is only fatal because
:math:`10^{-16}` is what a ``float64`` coefficient carries.  Form
:math:`I_2, I_3, I_4` in *double-double* arithmetic --- each number a pair
of ``float64`` limbs, together holding some 32 digits --- and the same
:math:`2.3\times10^9` amplification acts on :math:`10^{-32}` instead,
landing at :math:`10^{-23}`.  The roots are then limited by the final
rounding of the answer to ``float64`` rather than by the algebra at all.
One Aberth sweep, also in double-double, carries the quartic's roots to that
limit.  This is the default since 1.13.0, and
:data:`oscprob4nu.ROOT_STRATEGY` chooses between the two.

What each gains, and what the alternatives gain
"""""""""""""""""""""""""""""""""""""""""""""""

Measured against ground truth from ``mpmath`` at fifty decimal digits, worst
of the nine Hamiltonians in ``tests/stiff_reference.json``, which run from
the physical 3+1 range out to :math:`\Delta m^2_{41} = 1000` eV\ :sup:`2`
and down to a pair separated by :math:`10^{-16}`.  Cost is a full
:func:`~oscprob4nu.probabilities_4nu` over a 100 000-point stiff stack
through the compiled kernel, relative to what 1.12.0 shipped:

.. list-table::
   :header-rows: 1
   :widths: 40 20 15 25

   * - Strategy for the latent roots
     - Relative error
     - Cost
     - Keeps the closed form?
   * - Closed form alone
     - 2.2e-07
     - ---
     - yes
   * - Closed form + one Newton step
     - 3.9e-16
     - 1.00x
     - yes
   * - ``numpy.linalg.eigvalsh`` alone
     - 6.9e-16
     - 0.82x
     - no
   * - **Invariants in double-double**
     - **3.6e-17**
     - 1.25x
     - yes

The first row has no cost against it because it is no longer a route a
caller can select: :func:`~oscprob4nu.psi_roots_4nu` still solves the
quartic in closed form, and its contract is unchanged, but nothing in the
probability path is built on it any more.  Its error is quoted on the same
nine Hamiltonians as the rest, which is what it is there to show.

Four things in that table are worth reading twice.

The Newton step is **more accurate than LAPACK**, by about a factor of two.
That is not a fluke: ``eigvalsh`` reduces the matrix by Householder and QR
similarity transforms, each carrying a backward error of order
:math:`\epsilon \|H\|`, while the Newton step converges onto the root of
:math:`\det(\psi\mathbb{1} - \tilde{H})` for the matrix it was handed.

Double-double is **the only row that is not fighting the conditioning**.
Every other entry accepts a :math:`10^{-16}` coefficient and then tries to
recover from it, whether by refining against the matrix or by avoiding the
polynomial altogether; this one denies the amplification anything to work
on.  It is also the only row whose accuracy is set by the *output* format
rather than by the method, which is why more Aberth sweeps do not improve
it: one, two and three give bitwise identical answers on all nine
Hamiltonians.  That one suffices is a property of the *start*, not of the
iteration --- the eigensolver's quartet with the residual-trace shift
removed in double-double rather than in ``float64``.  Removing it in
``float64`` throws away the low limb, puts the start an ulp out, and costs
a second sweep to recover; the count was two for exactly as long as that
was how it was done.

Extended precision is still a **poor trade**, and instructively so.  Solving
the same quartic in ``numpy.longdouble`` was measured at 4.5e-11 --- under
one extra digit, rather than the three its wider mantissa suggests, because
the cluster amplifies whatever coefficient error survives --- while being
slower for not being hardware-vectorised, and silently ``float64`` on Apple
Silicon and on Windows, where the "fix" would quietly do nothing at all.
Double-double gets four orders more than that, portably, out of pairs of the
``float64`` every machine has.

The double-double route does call ``eigvalsh``, once per element, and only
for a **starting quartet**.  That is a smaller concession than it looks:
what a start has to supply is not accuracy but four distinct basins, and its
own error is discarded by the first sweep.  Euler's closed form is twice as
fast and exact on every stiff case in the reference, and still cannot be
used, because on a cluster Aberth converges *linearly* at ratio one half ---
from :math:`10^{-7}` five sweeps reach :math:`3.8\times10^{-9}`, and it
would need some thirty.  A backward-stable Hermitian eigensolver separates a
cluster as well as ``float64`` allows, in one call.  Durand-Kerner was
measured in place of Aberth as well, at one division per root against five,
and rejected for being non-monotone in the sweep count: 3.9e-16, then
9.7e-17, then 1.9e-16.

The honest summary on cost is unchanged in shape: four flavors costs more
per point than three, and the accurate routes cost more than the fast ones.

Why the refinement is not applied selectively
"""""""""""""""""""""""""""""""""""""""""""""

The obvious saving is to refine only the elements that need it, the way
:data:`oscprob3nu.SMALL_BATCH` and :data:`fastkernels.MIN_BATCH` dispatch on
a measured threshold.  It was measured, and it does not work.  The result is
recorded here so that it is not rediscovered.

Two criteria were tried, on 6300 Hamiltonians spanning clustered, doubly
paired and generic spectra, against ``eigvalsh``.  A criterion is *safe* at a
given cut only if every element below the cut is more accurate than the
target; the question is how many elements a safe cut can skip.

* **The gap-based amplification** :math:`\max_m
  |\psi|^3_{\max}/|\chi'(\psi_m)|`, which is what perturbation theory
  suggests, since the root sensitivity goes as :math:`1/\chi'`.  It predicts
  the error well for a single cluster and badly for two degenerate *pairs*,
  a family it does not model: those reach :math:`1.7\times10^{-10}` at an
  amplification of ten, where the criterion expects round-off.  The largest
  safe cut is about 2.3, which is close to the smallest value the indicator
  ever takes --- so it skips nothing.
* **A matrix residual**, comparing :math:`\prod_m \psi_m` against
  :math:`\det \tilde{H}`, which is :math:`\chi` evaluated at zero and costs
  one determinant instead of four.  Being built from the matrix it cannot
  lose information the way a gap heuristic does, but it is one scalar
  constraint on four roots, and errors cancel in the product: there are
  samples with a residual of :math:`10^{-17}` and a root error of
  :math:`5\times10^{-5}`.  It also skips nothing safely.

The second failure points at the general reason.  A criterion complete enough
to certify all four roots has to evaluate :math:`\chi` at all four roots ---
and that *is* the refinement.  The check and the fix are the same
computation, so there is nothing to save by doing the check first.

That conclusion held up when it was retested for the double-double route,
where it cuts the other way and hands over a free certificate.  An Aberth
sweep evaluates :math:`\chi` at four roots by construction, so the size of
the step it takes says whether the answer has converged --- around
:math:`10^{-32}` relative once it has, against :math:`10^{-9}` while still
crawling up a cluster, a gap wide enough that any threshold between them
works.  That made a genuinely cheap arrangement possible: start from Euler's
closed form, which costs nothing, and fall back to the eigensolver only for
the elements the certificate rejects.  It was built and measured, and
dropped --- the fallback fires on 80% of stiff spectra and 99% of those with
a tight pair, which leaves it slower than simply starting from the
eigensolver every time, for twice the code.  The criterion was sound; the
saving it was guarding was not there.

Note also who would benefit.  The four-flavor module exists mainly for 3+1,
and a 3+1 scan is stiff at every point, so even a working criterion would
skip nothing on the workload that motivates the module, while adding its own
cost.  Unconditional refinement is therefore not a compromise: it is what the
measurement supports.

Finally, this is specific to four flavors rather than a general caveat.  The
same measurement on :mod:`oscprob3nu` gives :math:`10^{-14}`, because there
:math:`\Delta m^2_{31}/\Delta m^2_{21}` is 34 rather than 13500.

That figure has since been reached independently.  Compared against
`nuSQuIDS <https://github.com/arguelles/nuSQuIDS>`_, which integrates the
density matrix numerically and shares nothing with this expansion, the
four-flavor probabilities agree to :math:`4\times10^{-16}` where the spectrum
is benign and :math:`3\times10^{-10}` where it is stiffest --- the latter
being this limit, arrived at from outside rather than asserted from within.
`Notebook 17
<https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/17_cross_checks.ipynb>`_
attributes that residual: against ``scipy.linalg.expm`` on the same
Hamiltonian our error is the same size, so it is ours and not a
disagreement.

Degenerate spectra at four flavors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Reconstructing :math:`U_4` from its roots divides by their differences, so a
repeated root needs care.  Rather than branch on a tolerance, the exponential
is interpolated over the roots in **Newton form**, with divided differences:
a repeated node is then a derivative, and for the exponential that derivative
is known exactly, :math:`f^{(k)}(\psi)/k! = (-iL)^k e^{-i\psi L}/k!`.

The alternative --- solving the Vandermonde system for the Cayley-Hamilton
coefficients --- is singular the moment two roots coincide, which includes a
Hamiltonian proportional to the identity and any triply degenerate spectrum.
The Newton form handles every degenerate case with no special branch at all.

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
     - ~21x
     - one ``eigh`` plus phases
   * - Versus energy (many :math:`H`, one :math:`L`)
     - ~23x
     - batched ``numpy.linalg.eigh``
   * - Oscillogram, 100 x 100
     - ~37x
     -
   * - Two flavors, versus baseline
     - ~99x
     -

These are the ratios of the same four measurements tabulated on the
:doc:`landing page <index>` and in ``README.md``, which state them as
absolute timings; the numbers here are those timings divided.  They used to
be quoted independently and had drifted apart --- 30x against 21x, and 70x
against 99x --- which is why they are now derived from one set of figures
and guarded together by ``tests/test_documented_figures.py``.

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

The thresholds are measured, not guessed: thirteen elements for three
flavors, twelve for two.  They sit close together even though the two-flavor
expansion does much less work per element, because what has to be amortised
is the array machinery's fixed cost rather than the arithmetic.  Nothing
about this is visible from the outside; the answers are the same either way.

The compiled backend
^^^^^^^^^^^^^^^^^^^^

With :mod:`fastkernels` --- that is, with Numba, which is installed with the
package as of 1.13.0 --- the batched paths are compiled instead.  The NumPy path evaluates the expansion as a
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
   * - 200 000 energies, four flavors
     - ~19x
   * - 20 000 energies, four flavors
     - ~18x
   * - 200 000 energies, three flavors
     - ~15x
   * - 20 000 energies, three flavors
     - ~9x
   * - 100 x 100 oscillogram
     - ~3.5x
   * - 200 000 baselines, two flavors
     - ~1.5x

Four flavors gains the most, and the reason is worth stating, because it is
not that the kernel is cleverer there.  That expansion needs a quartic, a
Newton refinement of its four roots against the matrix, and a Newton-form
reconstruction; as whole-array operations that is some forty passes rather
than fifteen, one of them a batched :math:`4\times4` determinant.  Done one
element at a time none of it leaves the registers, so the path carrying the
most fixed cost is the one with the most to shed.

The two-flavor row is the honest caveat: that path reduces to a square root
and a sine per element, which NumPy already does about as well as compiled
code can, and the kernel additionally has to materialise the Hamiltonian
stack --- for a scan over baselines, the same matrix repeated, which costs
2.5 ms to copy at two hundred thousand points.

So the backend is not used unconditionally.  :func:`fastkernels.worthwhile`
declines it below the per-flavor thresholds in
:data:`fastkernels.MIN_BATCH`, which were found by alternating the two paths
and taking the best of nine rounds each: for three and four flavors the
kernel wins at every size, for two it wins from about fifty thousand
elements.  A backend
that is sometimes slower than the path it replaces is worse than no backend,
and a test asserts the thresholds are honoured.

The scalar path is deliberately left uncompiled.  One probability takes
about eight microseconds; compiling it would save most of that, at the cost of
a multi-second pause on a user's first call.

Checking the input, and what it costs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The expansion assumes a Hermitian Hamiltonian, and one that is not Hermitian
does not fail loudly: the probabilities it returns still sum to one, so the
check a caller would actually apply cannot tell that the answer is
meaningless.  Every entry point therefore verifies it, and the tolerance is
relative to the largest entry, so a matrix assembled in floating point passes
--- everything the sample-Hamiltonian modules build is Hermitian to about
:math:`2 \times 10^{-17}` relative, against a tolerance of :math:`10^{-12}`.

The cost is stated rather than buried, because it is larger than it looks.
Validating a stack is a pass over it, which is the same order of work as
evaluating it, and the compiled kernel has made evaluating it fast; on a
large stack the check therefore dominates.  Interleaved, best of fifteen
rounds each:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Stack
     - Two flavors
     - Four flavors
   * - 2 000 points
     - 1.5x
     - 1.3x
   * - 200 000 points
     - 5.7x
     - 3.2x

Three flavors sits between them, at 1.8x and 3.9x.  Two ways of making it
cheaper were measured.  Comparing real and imaginary parts separately, rather
than forming :math:`H - H^\dagger`, avoids a temporary the size of the stack
and a square root per element, and is what the code does.  Replacing
``np.abs(...).max()`` with reductions that allocate nothing was the obvious
next step and came out 1.4x *slower*, because those views are strided over
the complex array while ``np.abs`` reads it contiguously.

It defaults to on regardless, because a library that silently returns
meaningless numbers costs its user more than the check does.
:data:`oscprob3nu.CHECK_HERMITICITY` turns it off, per module.

Degenerate spectra on the batched path
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The degenerate branch cannot be taken elementwise inside a vectorised
expression.  The general formula is therefore evaluated everywhere, with
vanishing denominators replaced by one, and the affected elements are then
recomputed individually.  Degeneracy is measure-zero among floating-point
Hamiltonians, so that fallback loop is empty in essentially every real use.
