// Prob3++ on the constant-density scan of the performance figure.
//
// Same problem as the frozen rows: P(nu_mu -> nu_e) at L = 1300 km through
// constant matter of 3 g/cm^3, three flavors, N energies from 0.1 to
// 31.6 GeV.  propagateLinear takes one energy, so the scan is a loop.
//
// SetMNS wants Dm2_32, not Dm2_31.  Density carries the 0.9924922 ratio
// between this library's normalisation constant and Prob3++'s 1.52588e-4;
// on the linear path it is SetDensityConversion that applies the electron
// fraction, so it is set to 0.5 times that ratio.

#include <cstdio>
#include <cmath>
#include <chrono>

#include "BargerPropagator.h"

static const double POT_RATIO = 1.514423e-4/1.52588e-4;

int main()
{
    const int sizes[] = {1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000};
    BargerPropagator bnu;
    bnu.SetWarningSuppression(true);
    bnu.SetDensityConversion(0.5*POT_RATIO);

    for (int k = 0; k < 10; k++)
    {
        const int n = sizes[k];
        double best = 1e30;
        for (int rep = 0; rep < 5; rep++)
        {
            auto t0 = std::chrono::high_resolution_clock::now();
            volatile double sink = 0.0;
            for (int i = 0; i < n; i++)
            {
                double E = pow(10.0, -1.0 + 2.5*(n == 1 ? 0.0
                                                 : (double)i/(n-1)));
                bnu.SetMNS(0.310, 2.240e-2, 0.582, 7.39e-5,
                           2.525e-3 - 7.39e-5, 217.0*M_PI/180.0,
                           E, true, 1);
                bnu.propagateLinear(1, 1300.0, 3.0);
                sink += bnu.GetProb(2, 1);
            }
            auto t1 = std::chrono::high_resolution_clock::now();
            (void)sink;
            double s = std::chrono::duration<double>(t1-t0).count();
            if (s < best) best = s;
        }
        printf("%d %.9g\n", n, best);
        fflush(stdout);
    }
    return 0;
}
