Tutorial notebooks
==================

Twenty worked notebooks live in `notebooks/
<https://github.com/mbustama/NuOscProbExact/tree/main/notebooks>`_,
numbered in reading order.  Each carries its figures inline, so they can be
read on GitHub without being run, and each ends with a footer pointing at
the previous notebook, the next one, and the API reference.

They are the long form of :doc:`recipes`.  A recipe is a few lines and a
figure; a notebook is the same calculation with the reasoning around it ---
why the convention is what it is, what happens at the edges, and what the
numbers were checked against.  The figures on both are produced by the same
code, so there is no third version to drift out of step.

To run them rather than read them::

   pip install "nuoscprobexact[notebooks]"
   jupyter lab notebooks/

.. note::

   The notebooks are not built into this documentation --- executing twenty
   of them on every docs build would be slow, and they are more useful
   where their outputs are already stored.  The links below go to GitHub,
   which renders them with their figures.

Start here
----------

The conventions everything else assumes, and the two plots every treatment
of oscillations opens with.

`01. Basics <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/01_basics.ipynb>`_
   Units, one probability, and why to pass arrays.

`02. Oscillations in vacuum <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/02_vacuum_oscillations.ipynb>`_
   The probabilities against baseline and against energy.

Matter, and new physics
-----------------------

Once a matter potential is added the Hamiltonian stops being diagonal in the
mass basis, and the conventions start to bite.

`03. Matter, NSI, and LIV <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/03_matter_nsi_liv.ipynb>`_
   Constant-density matter and two kinds of new physics.

`04. Oscillograms <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/04_oscillogram.ipynb>`_
   A two-dimensional map in a single call.

`05. Bi-probability plots <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/05_biprobability.ipynb>`_
   CP violation, as an ellipse.

The Earth, and other profiles
-----------------------------

A varying density is handled by slabbing it: each piece is solved exactly
and the operators are composed.

`06. The Earth: PREM, chords, and slabs <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/06_earth_and_prem.ipynb>`_
   How a varying profile becomes a sequence of exact pieces.

`07. Probabilities through the Earth <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/07_earth_probabilities.ipynb>`_
   Zenith-angle scans, an Earth oscillogram, and real baselines.

`08. Unusual density profiles <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/08_unusual_density_profiles.ipynb>`_
   Castle walls, and why the arrangement of matter matters.

Speed, and the published figures
--------------------------------

What the library costs, and the figures of the paper it implements.

`09. Performance <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/09_performance.ipynb>`_
   Looping versus broadcasting, measured on your machine.

`10. The paper's figures <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/10_paper_figures.ipynb>`_
   The two figures of arXiv:1904.12391, reproduced.

Where the approximations fail
-----------------------------

The reason to want an exact result: three places the familiar formulas do
not reach.

`11. Exact versus the approximations <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/11_exact_vs_approximations.ipynb>`_
   Where the familiar formulas break, and by how much.

`12. Mass ordering and the octant <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/12_ordering_and_octant.ipynb>`_
   Two open questions, and how they show up.

`13. Antineutrinos, done properly <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/13_antineutrinos.ipynb>`_
   Conjugate *and* flip -- and two ways to get it half right.

Hard cases
----------

Resonances, degeneracies, a fourth flavor, and corroboration from codes that
share none of this one.

`14. Solar neutrinos and the MSW resonance <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/14_solar_and_adiabatic_msw.ipynb>`_
   The adiabatic resonance, and the limits of slabbing.

`15. Numerical edge cases <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/15_numerical_edge_cases.ipynb>`_
   Degeneracies, and what returns a number instead of NaN.

`16. Four neutrinos, and a sterile state <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/16_four_neutrinos.ipynb>`_
   A 3+1 system through SU(4), and why the method stops there.

`17. Cross-checks with other codes <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/17_cross_checks.ipynb>`_
   Corroboration from nuSQuIDS and Zaglauer-Schwarzer.

Underneath, and beyond
----------------------

The evolution operator itself, the library in motion, and a Hamiltonian of
your own.

`18. The evolution operator and the SU(n) coefficients <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/18_evolution_operator.ipynb>`_
   The machinery underneath, for composing and extending.

`19. Animated scenes <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/19_animations.ipynb>`_
   Four sweeps, as stills, and the reel they came from.

`20. An arbitrary Hamiltonian, through three profiles <https://github.com/mbustama/NuOscProbExact/blob/main/notebooks/20_arbitrary_hamiltonian.ipynb>`_
   A long-range force carried through an invented body, an exponential one,
   and the Earth.
