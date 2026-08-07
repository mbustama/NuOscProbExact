// NuFast-LBL on the constant-density scan of the performance figure, at two
// Newton settings.
//
// The frozen rows for this code were taken at N_Newton = 0, which is its
// default and a truncation rather than an error: it stops at 1.5e-5, where
// two Newton steps reach 8.3e-12.  Timing only the default understates what
// the accurate configuration costs, so both are measured here.
//
// Same problem as the other rows: P(nu_mu -> nu_e) at L = 1300 km through
// constant matter of 3 g/cm^3, three flavors, N energies from 0.1 to
// 31.6 GeV.  rhoYe = 3 * 0.5 * 0.99209 carries the mass defect.  The
// baseline is 1299.999945 km, which is the length in km that reproduces
// this library's L in eV^-1 given NuFast's hbar*c.

#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <chrono>

// NuFast_LBL.cpp ships its own main(); rename it out of the way rather
// than copy the routine out of the file, so that what is timed is exactly
// the released source.
#define main nufast_lbl_example_main
#include "NuFast_LBL.cpp"
#undef main

int main()
{
    const int sizes[] = {1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000};
    const int newtons[] = {0, 2};

    // Timed in alternated pairs, not in blocks: run as two blocks this
    // reported N_Newton = 2 as *faster* than N_Newton = 0, which is not
    // something more Newton steps can be, and is the same ordering artefact
    // that turned up between the two four-flavor root strategies.
    for (int k = 0; k < 10; k++)
    {
        const int n = sizes[k];
        std::vector<double> E(n);
        for (int i = 0; i < n; i++)
            E[i] = pow(10.0, -1.0 + 2.5*(n == 1 ? 0.0 : (double)i/(n-1)));
        std::vector<std::array<std::array<double, 3>, 3>> probs(n);

        double best[2] = {1e30, 1e30};
        for (int rep = 0; rep < 9; rep++)
        {
            for (int j = 0; j < 2; j++)
            {
                const int w = (rep % 2 == 0) ? j : 1 - j;
                auto t0 = std::chrono::high_resolution_clock::now();
                NuFast::Probability_Matter_LBL(
                    0.310, 2.240e-2, 0.582, 217.0*M_PI/180.0,
                    7.39e-5, 2.525e-3, 1299.999945, E,
                    3.0*0.5*0.99209238, newtons[w], probs);
                auto t1 = std::chrono::high_resolution_clock::now();
                double s = std::chrono::duration<double>(t1-t0).count();
                if (s < best[w]) best[w] = s;
            }
        }
        printf("%d %d %.9g\n", newtons[0], n, best[0]);
        printf("%d %d %.9g\n", newtons[1], n, best[1]);
        fflush(stdout);
    }
    return 0;
}
