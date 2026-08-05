// NuFast-Earth on the same PREM chord as everything else.
//
// The comparison is made fair in the same three places as for nuSQuIDS and
// nuCraft:
//
//   Profile.  NuFast-Earth's own PREM_NDiscontinuityLayer carries its own
//   PREM and its own Y_e.  The subclass below keeps its *layering scheme*
//   exactly -- four major layers, each cut into n equal sub-shells, each
//   held at its midpoint density -- but evaluates the density with this
//   library's PREM at Y_e = 0.5.  So n means what it means in NuFast-Earth,
//   and the Earth is ours.
//
//   Matter potential.  Its YerhoE2a = 1.526493231029146e-4 is the
//   atomic-mass-unit constant; this library's equivalent is 1.514423e-4.
//   The ratio is the nuclear mass defect, 0.99209.
//
//   Length.  Its km is 1e3/1.97327e-7 = 5.0677302143e9 eV^-1 against this
//   library's 5.06773e9, an excess of 4.23e-8 relative.  The chord is
//   L = -2 R cos(theta_z), so the cosine handed over is shrunk by that.
//
// Modes:  (none) sweep n;  --profile dump rho(r);  --vacuum;  --scan.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <chrono>

#include "NuFastEarth.h"
#include "Matrix.h"
#include "our_prem.h"

using namespace NuFast;

static const double COSZ          = -0.9;
static const double YE            = 0.5;
static const double DENSITY_SCALE = 0.99209238;
static const double COSZ_SCALE    = 5.0677300000e9/(1.0e3/1.97327e-7);

class OurPREM : public Earth_Density
{
    public:
        OurPREM(int n_per_layer, double scale)
        {
            layers[0] = 1221.5; layers[1] = 3480.0;
            layers[2] = 5701.0; layers[3] = 6371.0;
            n_per = n_per_layer;
            n_discontinuities = 4*n_per;
            discontinuities.resize(n_discontinuities);
            rhoYes.resize(n_discontinuities);
            int count = 0;
            double lo = 0.0;
            for (int L = 0; L < 4; L++)
            {
                const double hi = layers[L];
                for (int i = 0; i < n_per; i++)
                {
                    discontinuities[count] = lo + (i+1)*(hi-lo)/n_per;
                    const double r         = lo + (i+0.5)*(hi-lo)/n_per;
                    rhoYes[count]          = our_prem_rho(r)*YE*scale;
                    count++;
                }
                lo = hi;
            }
            constant_shells = true;
        }
        double rhoYe(double r)
        {
            if (r > layers[3]) return 0.0;
            for (int k = 0; k < n_discontinuities; k++)
                if (r <= discontinuities[k]) return rhoYes[k];
            return rhoYes[n_discontinuities-1];
        }
    private:
        int n_per;
        double layers[4];
        std::vector<double> rhoYes;
};

static std::vector<double> energies()
{
    std::vector<double> Es;
    for (int i = 0; i < 12; i++) Es.push_back(3.0*pow(40.0/3.0, i/11.0));
    return Es;
}

static void run(int n, double scale, std::vector<double> &out, double *us)
{
    std::vector<double> Es = energies();
    std::vector<double> coszs = {COSZ*COSZ_SCALE};
    OurPREM prem(n, scale);
    Probability_Engine engine;
    engine.Set_Oscillation_Parameters(0.310, 2.240e-2, 0.582,
                                      217.0*M_PI/180.0, 7.39e-5,
                                      2.525e-3, true);
    engine.Set_Earth(0.0, &prem);
    engine.Set_Production_Height(0.0);
    engine.Set_Eigenvalue_Precision(-1);
    engine.Set_Spectra(Es, coszs);
    std::vector<std::vector<Matrix3r>> p = engine.Get_Probabilities();
    out.clear();
    for (unsigned int i = 0; i < Es.size(); i++)
        out.push_back(p[i][0].arr[1][1]);

    if (!us) return;
    double best = 1e30;
    for (int rep = 0; rep < 5; rep++)
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        OurPREM prem2(n, scale);
        Probability_Engine e2;
        e2.Set_Oscillation_Parameters(0.310, 2.240e-2, 0.582,
                                      217.0*M_PI/180.0, 7.39e-5,
                                      2.525e-3, true);
        e2.Set_Earth(0.0, &prem2);
        e2.Set_Production_Height(0.0);
        e2.Set_Eigenvalue_Precision(-1);
        e2.Set_Spectra(Es, coszs);
        volatile double sink = e2.Get_Probabilities()[0][0].arr[1][1];
        (void)sink;
        auto t1 = std::chrono::high_resolution_clock::now();
        double t = std::chrono::duration<double, std::micro>(t1-t0).count()
                   / Es.size();
        if (t < best) best = t;
    }
    *us = best;
}

int main(int argc, char **argv)
{
    const char *mode = (argc > 1) ? argv[1] : "";

    if (strcmp(mode, "--profile") == 0)
    {
        for (double r = 0.0; r <= 6371.0; r += 1.0)
            printf("%.6f %.15g\n", r, our_prem_rho(r));
        return 0;
    }

    if (strcmp(mode, "--vacuum") == 0)
    {
        std::vector<double> p;
        run(64, 0.0, p, nullptr);            // zero density: vacuum
        std::vector<double> Es = energies();
        for (unsigned int i = 0; i < Es.size(); i++)
            printf("%.15g %.15g\n", Es[i], p[i]);
        return 0;
    }

    if (strcmp(mode, "--scan") == 0)
    {
        // extra multiplicative factor on the density, at fixed n
        double lo = atof(argv[2]), hi = atof(argv[3]);
        int steps = atoi(argv[4]);
        std::vector<double> p;
        for (int k = 0; k <= steps; k++)
        {
            double extra = lo + (hi-lo)*k/steps;
            run(256, DENSITY_SCALE*extra, p, nullptr);
            printf("%.10f", extra);
            for (unsigned int i = 0; i < p.size(); i++)
                printf(" %.15g", p[i]);
            printf("\n");
        }
        return 0;
    }

    const int ns[] = {1, 2, 4, 8, 16, 32, 64, 128, 256};
    for (int idx = 0; idx < 9; idx++)
    {
        std::vector<double> p;
        double us = 0.0;
        run(ns[idx], DENSITY_SCALE, p, &us);
        printf("%d %.6f", ns[idx], us);
        for (unsigned int i = 0; i < p.size(); i++) printf(" %.15g", p[i]);
        printf("\n");
    }
    return 0;
}
