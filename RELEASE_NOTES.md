The first release since **v1.0.0** (April 2019) — the version described in
[arXiv:1904.12391](https://arxiv.org/abs/1904.12391). Ten releases of work are folded in, and
this is the first to reach PyPI. The full history is in
[CHANGELOG.md](https://github.com/mbustama/NuOscProbExact/blob/main/CHANGELOG.md).

## Upgrading from 1.0.0: two fixes change numbers

Both are two-flavor. If you have results from 1.0.0, these are the ones to check.

- **`probabilities_2nu` dropped the h₂ contribution.** It computed
  P<sub>eμ</sub> = |h₁|²/|h|² sin²(|h|L), where the transition probability is
  |U<sub>μe</sub>|² = u₁² + u₂², so the h₂ term is needed as well. Since h₂ = −Im(H₁₂), it
  vanishes whenever the off-diagonal entry is real — so vacuum and constant-density matter
  were unaffected, and a CP-violating two-flavor Hamiltonian was not.
- **The two-flavor vacuum Hamiltonian used the opposite sign convention.** Built from
  M² = diag(Δm², −Δm²), it returned the negative of the textbook Hamiltonian. Invisible in
  vacuum; it matters as soon as a matter potential is added.

Three-flavor results from 1.0.0 are unaffected by both.

## What is new since 1.0.0

**Four flavors** *(1.9.0)* — `oscprob4nu` and `hamiltonians4nu` carry the closed form to
SU(4), bringing 3+1 sterile scenarios in as *closed* four-state systems rather than a leak
out of the three-flavor block. Four is where the method ends: at five flavors the eigenvalues
stop being expressible in radicals, which is a theorem rather than a missing feature.

**The Earth, and layered matter** *(1.8.0)* — `slabs` propagates across a sequence of
adjacent slabs of arbitrary width and density, solving each exactly and multiplying the
operators. `earth` builds those slabs from the Preliminary Reference Earth Model, along a
zenith angle or between two of fifteen named sites.

**Whole scans in one call** *(1.2.0)* — every core routine accepts a stack of Hamiltonians,
an array of baselines, or both broadcast against each other. Roughly 20× to 90× faster than
the equivalent Python loop, with no extra dependency.

**An optional compiled backend** *(1.6.0)* — with `numba` installed, the batched paths run as
compiled kernels. The answers are identical to round-off, and the backend is used only where
it has been measured to win.

**Input validation** *(1.11.0)* — a non-Hermitian Hamiltonian now raises rather than
returning probabilities that still sum to one, which is what made the mistake invisible.

Also since 1.0.0: a [documentation site](https://mbustama.github.io/NuOscProbExact/),
eighteen worked notebooks, and a regression suite of 596 tests at 100% coverage, cross-checked
against [nuSQuIDS](https://github.com/arguelles/nuSQuIDS) and against the Zaglauer–Schwarzer
closed form for the matter spectrum.

## Install

```shell
pip install nuoscprobexact
```

Python 3.9 or newer, tested on 3.9 through 3.13. `numpy` is the only requirement;
`pip install "nuoscprobexact[fast]"` adds the optional compiled backend.

## Links

- [Documentation](https://mbustama.github.io/NuOscProbExact/)
- [Quickstart](https://mbustama.github.io/NuOscProbExact/quickstart.html)
- [What it can compute, with code](https://mbustama.github.io/NuOscProbExact/recipes.html)
- [The paper](https://arxiv.org/abs/1904.12391)
