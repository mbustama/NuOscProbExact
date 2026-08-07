# The compiled drivers behind the Earth speed–accuracy figure

Three of the six codes in `prem_speed_accuracy.pdf` are C or C++ and cannot be
driven from Python. Their drivers live here, together with the raw output each
one produced, so that the figure can be reproduced or disputed without anyone
having to reconstruct what was run. `tests/prem_scan.py` reads the `.txt`
files; it does not need the codes themselves.

The previous round of this comparison built its drivers in a session
scratchpad and lost them. That is why these are in the repository.

## What every driver has to get right

All three are handed the *same* Earth as `earth.earth_slabs`: this library's
PREM, evaluated at Y_e = 0.5, on the layering scheme each code natively uses —
four major layers (1221.5, 3480, 5701, 6371 km), each cut into `n` equal
sub-shells held at their midpoint density. Sweeping `n` is the profile dial,
the counterpart of `n_slabs_per_segment`.

The PREM itself is emitted from `earth._PREM_COEFFS` into `our_prem.h` rather
than transcribed, and the generated header reproduces `earth.density_prem` to
5e-14 over 6372 radii.

Each code then needs two conventions matched, and each was verified twice —
once by derivation from the constant in its source, once by scanning for the
residual minimum:

| code | its constant | ratio to ours | scan minimum | vacuum vs this code |
|---|---|---|---|---|
| NuFast-Earth | `YerhoE2a = 1.526493231029146e-4` | 0.9920938 | 1.000000 | 1.2e-14 |
| GLoBES | `GLB_V_FACTOR = 7.63247e-14` | 0.992093 | 1.000000 | 1.8e-14 |
| Prob3++ | `tworttwoGf = 1.52588e-4` | 0.9924922 | — | 5.7e-4 |

All three are three roundings of the same atomic-mass-unit constant, and none
of them equals another; each needs its own factor.

For the length, all three carry hbar c = 1.97327e-7 eV m, so a kilometre is
5.0677302143e9 eV^-1 against this library's 5.06773e9. The chord is
L = -2 R cos(theta_z), so each is handed a cosine shrunk by 4.23e-8. The
vacuum column above is what that buys: without it, all three sit at ~1e-7.

## Per-code traps, all hit

* **GLoBES applies `GLB_Ne_MANTLE = 0.5` itself.** Handing it rho*Y_e
  double-counts the electron fraction and costs a factor of two in V_CC — an
  error of 0.94 in probability, which is how it was noticed.
* **`glbAllocParams()` before `glbInit()` segfaults.** No diagnostic.
* **Prob3++'s `SetMNS` takes Dm2_32, not Dm2_31.** Confirmed from `setmass`
  in `mosc.c`, which builds `mVac[2] = dms21 + dms23`. Passing Dm2_31 leaves
  0.26 in probability.
* **Prob3++'s Earth path uses the `Y_p` column of the profile file, not
  `SetDensityConversion`,** which only reaches `propagateLinear`. The
  mass-defect factor therefore goes into `Y_p`.
* **Prob3++ does not reach round-off in vacuum.** Its 5.7e-4 is not a
  convention error: the path length matches ours exactly and the mixing and
  splittings check out from source. It is the same floor as its 4.5e-5 at
  1300 km in the constant-density figure, scaled by this chord's 8.8x longer
  phase. It is why Prob3++ flattens at 3e-4 no matter how fine the shells.
* **NuFast-Earth's `Earth_Density` subclasses use `reserve` where they mean
  `resize`** and then index the empty vector. The subclass here uses `resize`.

## Building

```bash
# our_prem.h, from the library's own coefficients
python tests/external_drivers/gen_prem_header.py > our_prem.h

# NuFast-Earth -- https://github.com/PeterDenton/NuFast-Earth
git clone --depth 1 https://github.com/PeterDenton/NuFast-Earth.git
g++ -O3 -std=c++17 -c NuFast-Earth/src/*.cpp -INuFast-Earth/include
g++ -O3 -std=c++17 -INuFast-Earth/include -I. nufast_drv.cpp *.o -o nufast_drv
./nufast_drv > nufast_earth_prem.txt

# Prob3++ -- mosc.c and mosc3.c MUST be compiled as C, not C++
git clone --depth 1 https://github.com/rogerwendell/prob3plusplus.git
cd prob3plusplus && gcc -O2 -c mosc.c mosc3.c \
  && g++ -O2 -c EarthDensity.cc BargerPropagator.cc && cd ..
g++ -O2 -Iprob3plusplus -I. prob3_drv.cpp prob3plusplus/*.o -o prob3_drv
./prob3_drv 2>/dev/null | grep -v Loading > prob3_prem.txt

# GLoBES 3.2.18 -- needs GSL
curl -O https://www.mpi-hd.mpg.de/personalhomes/globes/download/globes-3.2.18.tar.gz
tar xzf globes-3.2.18.tar.gz && cd globes-3.2.18
./configure --prefix=$PREFIX --with-gsl-prefix=/usr && make && make install && cd ..
gcc -O2 -I$PREFIX/include -I. globes_drv.c -L$PREFIX/lib \
    -lglobes -lgsl -lgslcblas -lm -o globes_drv
LD_LIBRARY_PATH=$PREFIX/lib ./globes_drv > globes_prem.txt
```

Each driver also takes `--vacuum` and `--scan lo hi steps`, which are the two
checks tabulated above; `nufast_drv --profile` dumps the PREM it is using.

## The constant-density performance sweep

`globes_perf.c`, `prob3_perf.cpp` and `perf_nusquids.py` are the other
figure: cost against the number of energies in a scan, at L = 1300 km
through 3 g/cm^3, which is `performance.pdf`. They exist because none of
those three codes batches over energies, so their cost per probability is
flat in N, and a figure about batching needs that shown rather than
asserted. Build them the same way as the drivers above; each prints
`N seconds` per line, best of five.

nuCraft is not among them and cannot be: its interface takes a zenith
angle and propagates through its own Earth, with no constant-density mode.

## Output format

One row per shell count: `n  microseconds_per_probability  P(1) ... P(12)`,
for the twelve energies `numpy.logspace(log10(3), log10(40), 12)` GeV at
cos(theta_z) = -0.9, channel P(nu_mu -> nu_mu).
