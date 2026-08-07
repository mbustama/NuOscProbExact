/* GLoBES on the constant-density scan of the performance figure.
 *
 * Same problem as the rows already frozen in tests/timing_other_codes.json:
 * P(nu_mu -> nu_e) at L = 1300 km through constant matter of 3 g/cm^3,
 * three flavors, over N energies from 0.1 to 31.6 GeV.
 *
 * GLoBES exposes no batched entry point -- glbConstantDensityProbability
 * takes one energy -- so the scan is a loop, and the interest of the curve
 * is precisely that its cost per probability does not fall with N.
 *
 * Density carries the 0.992093 mass defect; GLoBES applies its own
 * GLB_Ne_MANTLE = 0.5, so no electron fraction is folded in here.
 */

#include <stdio.h>
#include <math.h>
#include <time.h>

#include <globes/globes.h>

#define POT_SCALE 0.992093

int main(int argc, char *argv[])
{
    const int sizes[] = {1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000};
    glb_params p;
    int i, k, rep;

    glbInit(argv[0]);
    p = glbAllocParams();
    glbDefineParams(p, asin(sqrt(0.310)), asin(sqrt(2.240e-2)),
                    asin(sqrt(0.582)), 217.0*M_PI/180.0, 7.39e-5, 2.525e-3);
    glbSetDensityParams(p, 1.0, GLB_ALL);
    glbSetOscillationParameters(p);

    for (k = 0; k < 10; k++)
    {
        const int n = sizes[k];
        double best = 1e30;
        for (rep = 0; rep < 5; rep++)
        {
            struct timespec t0, t1;
            volatile double sink = 0.0;
            clock_gettime(CLOCK_MONOTONIC, &t0);
            for (i = 0; i < n; i++)
            {
                double E = pow(10.0, -1.0 + 2.5*(n == 1 ? 0.0
                                                 : (double)i/(n-1)));
                sink += glbConstantDensityProbability(2, 1, +1, E, 1300.0,
                                                      3.0*POT_SCALE);
            }
            clock_gettime(CLOCK_MONOTONIC, &t1);
            (void)sink;
            double s = (t1.tv_sec-t0.tv_sec)
                       + (t1.tv_nsec-t0.tv_nsec)*1e-9;
            if (s < best) best = s;
        }
        printf("%d %.9g\n", n, best);
        fflush(stdout);
    }
    glbFreeParams(p);
    return 0;
}
