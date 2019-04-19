# NuOscProbExact
Code to compute exact two- and three-neutrino oscillation probabilities using SU(2) and SU(3) expansions


## What is NuOscProbExact?

**NuOscProbExact** is a Python implementation of the method to compute exact two-flavor and three-flavor neutrino oscillation probabilities for arbitrary time-independent Hamiltonians presented in the paper [1904.XXXXX](https://arxiv.org/abs/1904.XXXXX).  The method relies on expansions of the Hamiltonian and time-evolution operators in terms of SU(2) and SU(3) matrices in order to obtain concise, analytical, and exact expressions for the probabilities, that are also easy to implement and evaluate.  For details of the method, see the paper above.

**NuOscProbExact** was developed by Mauricio Bustamante, who also authored the paper [1904.XXXXX](https://arxiv.org/abs/1904.XXXXX).  If you use NuOscProbExact in your work, please follow the directions on [Citing](#citing).


## Requirements

**NuOscProbExact** is fully written in Python 3.  It uses standard modules that are available, sometimes by default, as part of most Python installations, either stand-alone or via Anaconda:

* The two core modules (`oscprob2nu.py` and `oscprob3nu.py`) require only `numpy` and `cmath`.  These are the bare minimum requirements.

* The modules containing example Hamiltonians (`hamiltonians2nu.py` and `hamiltonians3nu.py`) require only `numpy`, `cmath`, and `copy`

* The modules containing the test suites (`oscprob2nu_tests.py`, `oscprob3nu_tests.py`, `oscprob3nu_plotpaper.py`, and `oscprob_testsuite.py`) require only `numpy`, `cmath`, `copy`, and `matplotlib`


## Installation




## Quick start


## Bundled test suite


## Documentation and help

All of the modules provided in **NuOscProbExact** have been documented using Python docstrings.

To view the documentation of a module from within an interactive Python session, run, *e.g.*,
```python
import oscprob3nu

print(oscprob3nu.__doc__)
```
This will print to screen a description of what the module does (in this example, `oscprob3nu`) and a list of the functions that it contains, including a description of each.

To view the documentation of a particular function from within an interactive Python session, run, *e.g.*,
```python
import oscprob3nu

help(oscprob3nu.hamiltonian_3nu_coefficients)
```
This will print to screen a description of what the function does (in this example, `oscprob3nu.hamiltonian_3nu_coefficients`), a list and description of its input parameters, and a description of the values that it returns.


## Citing

If you use **NuOscProbExact** in your work, we ask you that you please cite the following paper: Mauricio Bustamante, *Exact neutrino oscillation probabilities with arbitrary time-independent Hamiltonians*, arXiv:1904.XXXXX.

Please consider using the LaTeX or BibTeX entries provided by INSPIRE.





