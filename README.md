# NuOscProbExact
Code to compute exact two- and three-neutrino oscillation probabilities using SU(2) and SU(3) expansions


## What is NuOscProbExact?

**NuOscProbExact** is a Python implementation of the method to compute exact two-flavor and three-flavor neutrino oscillation probabilities for arbitrary time-independent Hamiltonians presented in 1904.XXXXX.  The method relies on expansions of the Hamiltonian and time-evolution operators in terms of SU(2) and SU(3) matrices in order to obtain concise, analytical, and exact expressions for the probabilities, that are also easy to implement and evaluate.  For details of the method, see the paper above.

NuOscProbExact was developed by Mauricio Bustamante.  If you use NuOscProbExact in your work, please follow the directions on [Citing](#citing).


## Requirements

NuOscProbExact is fully written in Python 3.  It uses standard modules that are available, sometimes by default, as part of most Python installations, either stand-alone or via Anaconda:

* The two core modules (`oscprob2nu.py` and `oscprob3nu.py`) require only `numpy` and `cmath`.  These are the bare minimum requirements.

* The modules containing example Hamiltonians (`hamiltonians2nu.py` and `hamiltonians3nu.py`) require only `numpy`, `cmath`, and `copy`

* The modules containing the test suites (`oscprob2nu_tests.py`, `oscprob3nu_tests.py`, `oscprob3nu_plotpaper.py`, and `oscprob_testsuite.py`) require only `numpy`, `cmath`, `copy`, and `matplotlib`


## Installation



## Quick start


## Bundled test suite


## Documentation


## Citing

