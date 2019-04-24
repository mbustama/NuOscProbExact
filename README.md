| arXiv |
|:-----:|
|[![arXiv](https://img.shields.io/badge/arXiv-1904.XXXXX-orange.svg)](https://arXiv.org/abs/1904.XXXXX)


# NuOscProbExact
Code to compute exact two- and three-neutrino oscillation probabilities using SU(2) and SU(3) expansions


## What is NuOscProbExact?

**NuOscProbExact** is a Python implementation of the method to compute exact two-flavor and three-flavor neutrino oscillation probabilities for arbitrary time-independent Hamiltonians presented in the paper [arXiv:1904.XXXXX](https://arxiv.org/abs/1904.XXXXX).  The method relies on expansions of the Hamiltonian and time-evolution operators in terms of SU(2) and SU(3) matrices in order to obtain concise, analytical, and exact expressions for the probabilities, that are also easy to implement and evaluate.  For details of the method, see the paper above.

**NuOscProbExact** was developed by Mauricio Bustamante, who also authored the paper [arXiv:1904.XXXXX](https://arxiv.org/abs/1904.XXXXX).  If you use it in your work, please follow the directions on [Citing](#citing).


## Requirements

**NuOscProbExact** is fully written in Python 3.  It uses standard modules that are available, sometimes by default, as part of most Python installations, either stand-alone or via Anaconda:

* **Bare minimum requirements:** The two core modules (`oscprob2nu.py` and `oscprob3nu.py`) require only `numpy` and `cmath`.  These are the bare minimum requirements.

* **To use the bundled sample Hamiltonians:** The modules containing example Hamiltonians (`hamiltonians2nu.py` and `hamiltonians3nu.py`) require only `numpy`, `cmath`, and `copy`

* **To run the test suite:** The modules containing the test suites (`oscprob2nu_plot.py`, `oscprob3nu_plot.py`, `oscprob3nu_plotpaper.py`, and `run_testsuite.py`) require only `numpy`, `cmath`, `copy`, and `matplotlib`


## Installation

Because **NuOscProbExact** is written fully in Python, no compilation or linking is necessary.  The installation is simple and consists only in fetching the files from GitHub.

> **Python 2 compatibility:** The code was written and tested using Python 3.  Yet, because the core modules `oscprob2nu` and `oscprob3nu` use only native Python functions and popular modules, they might also run in Python 2.  However, this is currently untested.

Instructions:

1. In the file system where you would like to install **NuOscProbExact**, go to the directory where you would like the code to be downloaded, *e.g.*,
   ```shell
   cd /home/MyProjects
   ```

2. From there, fetch the code from the GitHub repository with
   ```shell
   git clone https://github.com/mbustama/NuOscProbExact.git
   ```
   (Alternatively, you can download the zip file from GitHub and uncompress it.)

   Doing this will create the directory `/home/MyProjects/NuOscProbExact`, with the following file structure:
   ```
   /NuOscProbExact/README.md           The file that you are reading
   /NuOscProbExact/run_testsuite.py    Run this to create test plots of the probabilities
   /NuOscProbExact/fig                 Contains plots generate by the test suite (initially empyty)
   /NuOscProbExact/src                 Contains the main source files
       ../hamiltonians2nu.py           Routines to compute example two-flavor Hamiltonians
       ../hamiltonians3nu.py           Routines to compute example three-flavor Hamiltonians
       ../globaldefs.py                Physical constants and unit-conversion constants
       ../oscprob2nu.py                Routines to compute the two-flavor probabilities
       ../oscprob3nu.py                Routines to compute the three-flavor probabilities
   /NuOscProbExact/test                Contains the source files to run the test suite
       ../matplotlibrc                 Customized matplotlib style file
       ../oscprob2nu_plot.py           Routines to generate two-flavor probability test plots
       ../oscprob3nu_plot.py           Routines to generate three-flavor probability test plots
       ../oscprob3nu_plotpaper.py      Routine to generate the plot shown in the paper
   ```
   Now you are ready to start using **NuOscProbExact**.

3. (Optional) Run the test suite
   ```shell
   cd /home/MyProjects/NuOscProbExact
   python run_testsuite.py
   ```
   Doing this will generate plots of the two-neutrino and three-neutrino probabilities *vs.* distance and *vs.* energy, for different oscillation scenarios.  It will also generate the plot of probabilities *vs.* energy that is included in the paper.  The plots are stored in the `NuOscProbExact/fig` directory.  The code `run_testsuite.py` calls routines defined in `oscprob2nu_plot.py`, `oscprob3nu_plot.py`, and `oscprob3nu_plotpaper.py`, located in the `NuOscProbExact/test/` directory.  Inspecting these files may help you in coding your own project.


## Usage and examples

There are only two core modules: `oscprobn2nu.py` and `oscprob3nu.py`.  Each one is stand-alone (except for the dependencies described [above](#requirements)).  To use either module in your code, copy it to your project's working directory, or add their location to the paths where your environment looks for modules, *e.g.*,
```python
import sys

sys.path.append('../src')
```

In the examples below, we focus mostly on `oscprob3nu`, but what we show applies to `oscprob2nu` as well.


### Basics

Most of the time, you will be only interested in computing oscillation probabilities, not in the intermediate steps of the method.  The functions to compute and return the probabilities are `probabilities_3nu`, for the three-neutrino case, and `probabilities_2nu`, for the two-neutrino case.  Below, we show how to use them.

#### Three-neutrino oscillations

The function to compute three-neutrino probabilities is `probabilities_3nu` in the module `oscprob3nu`.  It takes as input parameters the `hamiltonian`, in the form of a 3x3 Hermitian matrix, and the baseline `L`.

This function returns the list of probabilities `Pee` (nu_e --> nu_e), `Pem` (nu_e --> nu_mu), `Pet` (nu_e --> nu_tau), `Pme` (nu_mu --> nu_e), `Pmm` (nu_mu --> nu_mu), `Pmt` (nu_mu --> nu_tau), `Pte` (nu_tau --> nu_e), `Ptm` (nu_tau --> nu_mu), and `Ptt` (nu_tau --> nu_tau).

To use it, call
```python
import oscprob3nu

hamiltonian = [[H11, H12, H13], [H21, H22, H23], [H31, H32, H33]]
Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = oscprob3nu.probabilities_3nu(hamiltonian, L)
```

#### Two-neutrino oscillations

The function to compute two-neutrino probabilities is `probabilities_2nu` in the module `oscprob2nu`.  It takes as input parameters the `hamiltonian`, in the form of a 2x2 Hermitian matrix, and the baseline `L`.

This function returns the list of probabilities `Pee` (nu_e --> nu_e), `Pem` (nu_e --> nu_mu), `Pme` (nu_mu --> nu_e), and `Pmm` (nu_mu --> nu_mu).  (These probabilities could also be `Pmm`, `Pmt`, `Pmt`, and `Ptt` instead, depending on what Hamiltonian you pass to `probabilities_2nu`.)

To use it, call
```python
import oscprob2nu

hamiltonian = [[H11, H12], [H21, H22]]
Pee, Pem, Pme, Pmm = oscprob2nu.probabilities_2nu(hamiltonian, L)
```

> **Important:** If you feed the code a non-Hermitian matrix, it will output nonsensical results

> **About the units:** The code in the modules `oscprob3nu` and `osprob2nu` does not assume units for any of the model parameters, so you need to make sure that you pass values with the correct units.  The module `globaldefs` contains conversion factors which might come in handy for this.


### Trivial example

#### Three-neutrino oscillations

As a first, trivial example, we pass an arbitrary Hamiltonian and baseline to `probabilities_3nu`:
```python
import oscprob3nu

hamiltonian = [
                [1.0+0.0j, 0.0+2.0j, 0.0-1.0j],
                [0.0-2.0j, 3.0+0.0j, 3.0+0.0j],
                [0.0+1.0j, 3.0-0.0j, 5.0+0.0j]
]

L = 1.0

Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = oscprob3nu.probabilities_3nu(hamiltonian, L)

print("Pee = %6.5f, Pem = %6.5f, Pet = %6.5f" % (Pee, Pem, Pet))
print("Pme = %6.5f, Pmm = %6.5f, Pmt = %6.5f" % (Pme, Pmm, Pmt))
print("Pte = %6.5f, Ptm = %6.5f, Ptt = %6.5f" % (Pte, Ptm, Ptt))
```
This returns
```shell
Pee = 0.34273, Pem = 0.41369, Pet = 0.24358
Pme = 0.41369, Pmm = 0.00485, Pmt = 0.58146
Pte = 0.24358, Ptm = 0.58146, Ptt = 0.17497
```

As expected, `Pem == Pme`, `Pet == Pte`, `Pmt == Ptm`, `Pee + Pem + Pet = 1`, `Pme + Pmm + Pmt = 1`, and `Pte + Ptm + Ptt = 1`.

#### Two-neutrino oscillations

In this case, we use `probabilities_2nu`:
```python
import oscprob2nu

hamiltonian = [
                [1.0+0.0j, 0.0+2.0j],
                [0.0-2.0j, 3.0+0.0j]
]

L = 1.0

Pee, Pem, Pme, Pmm = oscprob2nu.probabilities_2nu(hamiltonian, L)

print("Pee = %6.5f, Pem = %6.5f" % (Pee, Pem))
print("Pme = %6.5f, Pmm = %6.5f" % (Pme, Pmm))
```
This returns
```shell
```

As expected, `Pem == Pme` and `Pee + Pem = 1`.


### Oscillations in vacuum: fixed energy and baseline

#### Three-neutrino oscillations

Now we compute the three-neutrino oscillation probabilities in vacuum.  To do this, we can use the routine
```python
hamiltonian_3nu_vacuum_energy_independent(s12, s23, s13, dCP, D21, D31)
```
that is provided in the `hamiltonians3nu` module.  It returns the 3x3 Hamiltonian for oscillations in vacuum.  The input parameters `s12`, `s23`, `s13`, `dCP`, `D21`, and `D31` are, respectively, sin(theta_12), sin(theta_23), sin(theta_13), delta_CP, Delta m_21^2, and Delta m_31^2.  For this example, we set them to their current best-fit values, which we pull from `globaldefs` (inspect that file for more information about these values).

> **Important:** The function `hamiltonian_3nu_vacuum_energy_independent` returns the Hamiltonian in vacuum **without** multiplying it by the *1/E* prefactor, where *E* is the neutrino energy.  It was done in this way so that, if we wish to compute the probabilities at different energies, we need to compute `hamiltonian_3nu_vacuum_energy_independent` only once, and then multiply it by a varying *1/E* prefactor.

```python
import numpy as np

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(  S12_BF, S23_BF,
                                                                                    S13_BF, DCP_BF,
                                                                                    D21_BF, D31_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)

Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = oscprob3nu.probabilities_3nu( h_vacuum,
                                                                            baseline*CONV_KM_TO_INV_EV)

print("Pee = %6.5f, Pem = %6.5f, Pet = %6.5f" % (Pee, Pem, Pet))
print("Pme = %6.5f, Pmm = %6.5f, Pmt = %6.5f" % (Pme, Pmm, Pmt))
print("Pte = %6.5f, Ptm = %6.5f, Ptt = %6.5f" % (Pte, Ptm, Ptt))
````
This returns
```shell
Pee = 0.96711, Pem = 0.01593, Pet = 0.01695
Pme = 0.01823, Pmm = 0.64417, Pmt = 0.33761
Pte = 0.01466, Ptm = 0.33990, Ptt = 0.64544
```

Sometimes, you might be interested also in returning the coefficients `h1`, ..., `h8` of the expansion of the Hamiltonian in terms of Gell-Mann matrices (Table II in the paper), the coefficients `u0`, ..., `u8` of the SU(3) expansion of the associated time-evolution operator (Eqs. (13) and (14) in the paper), or the time-evolution operator `evol_operator` itself, as a 3x3 matrix (Eq. (15) in the paper).  See the paper [arXiv:1904.XXXXX](https://arxiv.org/abs/1904.XXXXX) for details on these quantities.  The module `oscprob3nu` has functions to do this:
```python
import numpy as np

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(  S12_BF, S23_BF,
                                                                                    S13_BF, DCP_BF,
                                                                                    D21_BF, D31_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)

h1, h2, h3, h4, h5, h6, h7, h8 = oscprob3nu.hamiltonian_3nu_coefficients(h_vacuum)
u0, u1, u2, u3, u4, u5, u6, u7, u8 = oscprob3nu.evolution_operator_3nu_u_coefficients(  \
                                                                            h_vacuum,
                                                                            baseline*CONV_KM_TO_INV_EV)
evol_operator = oscprob3nu.evolution_operator_3nu(h_vacuum, baseline*CONV_KM_TO_INV_EV)
```

#### Two-neutrino oscillations

To compute the two-neutrino oscillation probabilities in vacuum, we can use the routine
```python
hamiltonian_2nu_vacuum_energy_independent(sth, Dm2)
```
that is provided in the `hamiltonians2nu` module.  The input parameters `sth`, and `Dm2` are, respectively, sin(theta), and Delta m^2.  For this example, we set them to current best-fit values for atmospheric neutrinos.

```python
import numpy as np

import oscprob2nu
import hamiltonians2nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(S23_BF, D31_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)

Pee, Pem, Pme, Pmm = oscprob3nu.probabilities_2nu(h_vacuum, baseline*CONV_KM_TO_INV_EV)

print("Pee = %6.5f, Pem = %6.5f" % (Pee, Pem))
print("Pme = %6.5f, Pmm = %6.5f" % (Pme, Pmm))
````
This returns
```shell
```

Like in the three-neutrino case, we can also return the coefficients `h1`, `h2`, `h3` of the expansion of the Hamiltonian in terms of Pauli matrices (Table I in the paper), or the time-evolution operator `evol_operator` itself, as a 2x2 matrix (Eq. (5) in the paper).
```python
import numpy as np
from pylab import *
from matplotlib import *
import matplotlib as mpl

import oscprob2nu
import hamiltonians2nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians2nu.hamiltonian_2nu_vacuum_energy_independent(S23_BF, D31_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)

h1, h2, h3 = oscprob2nu.hamiltonian_2nu_coefficients(h_vacuum)
evol_operator = oscprob2nu.evolution_operator_2nu(h_vacuum, baseline*CONV_KM_TO_INV_EV)
```


### Three-neutrino oscillations in vacuum: fixed energy, varying baseline

Now we fix the energy at, say, 10 MeV, and vary the baseline between 1 and 500 km.  We use a fine grid in `L` so that the oscillations are clearly rendered.

```python
import numpy as np

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e7     # Neutrino energy [eV]

# Baselines, L
log10_l_min = 0.0  # log10 [km]
log10_l_max = 3.0  # log10 [km]
log10_l_npts = 1000
log10_l_val = np.linspace(log10_l_min, log10_l_max, log10_l_npts)  # [km]
l_val = [CONV_KM_TO_INV_EV*10.**x for x in log10_l_val]


h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(  S12_BF, S23_BF,
                                                                                    S13_BF, DCP_BF,
                                                                                    D21_BF, D31_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)

# Each element of prob: [Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt]
prob = [oscprob3nu.probabilities_3nu(h_vacuum, l) for l in l_val]
prob_ee = [x[0] for x in prob]  # Pee
prob_em = [x[1] for x in prob]  # Pem
prob_et = [x[2] for x in prob]  # Pet
```

To plot the data:
```python
from pylab import *
from matplotlib import *
import matplotlib as mpl

fig = plt.figure(figsize=[9,9])
ax = fig.add_subplot(1,1,1)
ax.set_xlabel(r'Baseline $L$ [km]', fontsize=25)
ax.set_ylabel(r'Three-neutrino probability', fontsize=25)
ax.legend(loc='center left', frameon=False)

ax.plot(l_val, prob_ee, label=r'$P_{\nu_e \to \nu_e}$', color='C0', zorder=1)
ax.plot(l_val, prob_em, label=r'$P_{\nu_e \to \nu_\mu}$', color='C1', zorder=1)
ax.plot(l_val, prob_et, label=r'$P_{\nu_e \to \nu_\tau}$', color='C2', zorder=1)
````

Alternatively, you can automatically produce plots of probability *vs.* baseline using the following function from the `oscprob3nu_tests` module:
```python
import oscprob3nu_tests

case = 'vacuum'
oscprob3nu_tests.plot_probability_3nu_vs_baseline(
                case, energy=1.e-2,
                log10_l_min=0.0, log10_l_max=3.0, log10_l_npts=1000,
                plot_prob_ee=True, plot_prob_em=True, plot_prob_et=True,
                plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=False,
                plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False,
                output_filename='prob_3nu_vacuum_vs_baseline_ee_em_et', output_format='png',
                legend_loc='center left', legend_ncol=1, path_save='../fig/')
```
The function `plot_probability_3nu_vs_baseline` assumes that `energy` is in GeV and the (log10) of the baselines `log10_l_min` and `log_l_max` are in km.  The function call above produces the following plot:
![prob_3nu_vacuum_vs_baseline_ee_em_et.png](https://github.com/mbustama/NuOscProbExact/blob/master/img/prob_3nu_vacuum_vs_baseline_ee_em_et.png)

The parameter `case` can take any of the following values:
* `vacuum`: for oscillations in vacuum, assuming the default values of mixing parameters from the `globaldefs` module
* `matter`: for oscillations in constant matter, assuming the density of the Earth's crust as set in `globaldefs`
* `nsi`: for oscillations in matter with non-standard interactions, with the NSI strengh parameters fixed to the default values in `globaldefs`
* `liv`: for oscillations in a CPT-odd Lorentz invariance-violating (LIV) background, with the LIV parameters fixed to the default values in `globaldefs`

For more information about these cases, see the paper [arXiv:1904.XXXXX](https://arxiv.org/abs/1904.XXXXX).  For more information about how to use `plot_probability_3nu_vs_baseline`, see the documentation of the function in the `oscprob3nu_tests` module.


### Three-neutrino oscillations in vacuum: fixed baseline, varying energy

Now we fix the baseline at, say, 1300 km, and vary the energy between 100 MeV and 10 GeV.

```python
import numpy as np

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

baseline = 1.3e3                       # Baseline [km]
baseline = baseline*CONV_KM_TO_INV_EV  # [eV^{-1}]

# Neutrino energies
log10_energy_min = -1.0 # [GeV]
log10_energy_max = 1.0  # [GeV]
log10_energy_npts = 200
log10_energy = np.linspace( log10_energy_min,
                            log10_energy_max,
                            log10_energy_npts)
energy = [10.**x for x in log10_energy] # [GeV]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_vacuum_energy_independent(  S12_BF, S23_BF,
                                                                                S13_BF, DCP_BF,
                                                                                D21_BF, D31_BF)

# Each element of prob: [Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt]
prob = [oscprob3nu.probabilities_3nu(np.multiply(1./x/1.e9, h_vacuum_energy_indep), baseline) \
        for x in energy]
prob_ee = [x[0] for x in prob]  # Pee
prob_em = [x[1] for x in prob]  # Pem
prob_et = [x[2] for x in prob]  # Pet
```

To visualize the data:


Alternatively, you can automatically produce plots of probability using the following function from the `oscprob3nu_tests` module:
```python
import oscprob3nu_tests

case = 'vacuum'
oscprob3nu_tests.plot_probability_3nu_vs_energy(
                                    case, baseline=1.3e3,
                                    output_filename='prob_3nu_vacuum_vs_energy', output_format='pdf',
                                    log10_energy_min=-1.0, log10_energy_max=1.0, log10_l_npts=200,
                                    plot_prob_ee=True, plot_prob_em=True, plot_prob_et=True,
                                    plot_prob_me=False, plot_prob_mm=False, plot_prob_mt=False,
                                    plot_prob_te=False, plot_prob_tm=False, plot_prob_tt=False)
```
The parameter `case` can take any of the same values as listed [above](#oscillations-in-vacuum-fixed-energy-varying-baseline).


### Three-neutrino oscillations in matter

For oscillation in matter, we proceed in an analogous way as for oscillations in vacuum.  To compute the Hamiltonian in matter, we can use the routine `hamiltonian_matter` in the module `hamiltonians3nu`.  First, we need to compute `h_vacuum`, and then pass it to `hamiltonian_matter`, together with the neutrino-electron charged-current potential `VCC`, with V_CC = sqrt(2.0) * G_F * n_e.  This routine is called as:
```python
hamiltonian_matter(h_vacuum, VCC)
```

In the example below, we set the matter potential to `VCC_EARTH_CRUST`, which is computed using the electron density of the crust of the Earth, and is read from `globaldefs`.
```python
import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(  S12_BF, S23_BF,
                                                                                    S13_BF, DCP_BF,
                                                                                    D21_BF, D31_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)
h_matter = hamiltonians3nu.hamiltonian_3nu_matter(h_vacuum, VCC_EARTH_CRUST)

Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = oscprob3nu.probabilities_3nu( h_matter,
                                                                            baseline*CONV_KM_TO_INV_EV)

print("Pee = %6.5f, Pem = %6.5f, Pet = %6.5f" % (Pee, Pem, Pet))
print("Pme = %6.5f, Pmm = %6.5f, Pmt = %6.5f" % (Pme, Pmm, Pmt))
print("Pte = %6.5f, Ptm = %6.5f, Ptt = %6.5f" % (Pte, Ptm, Ptt))
````
This returns
```shell
Pee = 0.99981, Pem = 0.00010, Pet = 0.00009
Pme = 0.00011, Pmm = 0.65964, Pmt = 0.34025
Pte = 0.00008, Ptm = 0.34026, Ptt = 0.65966
```

### Three-neutrino oscillations in matter with non-standard interactions (NSI)

For oscillation in matter with NSI, we can use the routine `hamiltonian_nsi` in the module `hamiltonians3nu`.  First, we need to compute `h_vacuum`, and then pass it to `hamiltonian_nsi`, together with `VCC`, and with a vector `eps` containing the NSI strength parameters,
```python
eps = [eps_ee, eps_em, eps_et, eps_mm, eps_mt, eps_tt]
```
This routine is called as
```python
hamiltonian_nsi(h_vacuum, VCC, eps)
```

In the example below, we set `eps` to its default value pulled from `globaldefs`:
```python
import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(  S12_BF, S23_BF,
                                                                                    S13_BF, DCP_BF,
                                                                                    D21_BF, D31_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)
h_nsi = hamiltonians3nu.hamiltonian_3nu_nsi(h_vacuum, VCC_EARTH_CRUST, EPS_TEST)

Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = oscprob3nu.probabilities_3nu( h_nsi,
                                                                            baseline*CONV_KM_TO_INV_EV)

print("Pee = %6.5f, Pem = %6.5f, Pet = %6.5f" % (Pee, Pem, Pet))
print("Pme = %6.5f, Pmm = %6.5f, Pmt = %6.5f" % (Pme, Pmm, Pmt))
print("Pte = %6.5f, Ptm = %6.5f, Ptt = %6.5f" % (Pte, Ptm, Ptt))
````
This returns
```shell
```


### Three-neutrino oscillations in a Lorentz invariance-violating (LIV) background

For oscillation LIV, we can use the routine `hamiltonian_liv` in the module `hamiltonians3nu`.  As before, first, we need to compute `h_vacuum`, and then pass it to `hamiltonian_liv`, together with the vector `liv_params` containing the LIV parameters:
```python
liv_params = [Lambda, b1, b2, b3, s12_liv, s23_liv, s13_liv, dCP_liv]
```
This routine is called as
```python
hamiltonian_liv(h_vacuum, liv_params)
```

In the example below, we set `liv_params` to its default value pulled from `globaldefs`:
```python
import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(  S12_BF, S23_BF,
                                                                                    S13_BF, DCP_BF,
                                                                                    D21_BF, D31_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)
h_liv = hamiltonians3nu.hamiltonian_3nu_liv(h_vacuum, LIV_PARAMS_TEST)

Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = oscprob3nu.probabilities_3nu( h_liv,
                                                                            baseline*CONV_KM_TO_INV_EV)

print("Pee = %6.5f, Pem = %6.5f, Pet = %6.5f" % (Pee, Pem, Pet))
print("Pme = %6.5f, Pmm = %6.5f, Pmt = %6.5f" % (Pme, Pmm, Pmt))
print("Pte = %6.5f, Ptm = %6.5f, Ptt = %6.5f" % (Pte, Ptm, Ptt))
````
This returns
```shell
```

### Arbitrary Hamiltonians

Of course, you can supply your custom Hamiltonian and compute the associated oscillation probabilities; see [Trivial example](#trivial-example) above.  Usually, you will want to add an extra term from your preferred model to the vacuum Hamiltonian.  To do that, take a cue from the examples above.

In the following example, the function `hamiltonian_mymodel`, supplied by you, should return a 3x3 matrix:
```python
import oscprob3nu
import hamiltonians3nu
from globaldefs import *

energy = 1.e9     # Neutrino energy [eV]
baseline = 1.3e3  # Baseline [km]

h_vacuum_energy_indep = hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(  S12_BF, S23_BF,
                                                                                    S13_BF, DCP_BF,
                                                                                    D21_BF, D31_BF)
h_vacuum = np.multiply(1./energy, h_vacuum_energy_indep)
h_mymodel = h_vacuum + hamiltonian_mymodel(mymodel_parameters)

Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = oscprob3nu.probabilities_3nu( h_mymodel,
                                                                            baseline*CONV_KM_TO_INV_EV)

```
Though we do not show it here, `hamiltonian_mymodel` could also depend on `energy`.  The code for two-neutrino oscillations is analogous, but `hamiltonian_mymodel` should return a 2x2 matrix instead.


## Documentation and help

All of the modules provided in **NuOscProbExact** have been documented using Python docstrings.  These are human-readable by opening the source `.py` files.  Alternatively, they can be printed from within an interactive Python session.

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

If you use **NuOscProbExact** in your work, we ask you that you please cite the following paper: Mauricio Bustamante, *Exact neutrino oscillation probabilities with arbitrary time-independent Hamiltonians*, [arXiv:1904.XXXXX](https://arxiv.org/abs/1904.XXXXX).

If you are citing **NuOscProbExact** in a document that will be uploaded to the arXiv, please consider using the LaTeX or BibTeX entries provided by INSPIRE.





