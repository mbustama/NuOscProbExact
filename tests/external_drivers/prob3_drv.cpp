// Prob3++ on the same PREM chord as everything else.
//
//   Profile.  EarthDensity's own table is the SK three-flavor profile, with
//   Y_p = 0.468 in the core and 0.497 in the mantle.  Its documented file
//   route takes "radius rho Yp" per shell, so this library's PREM is written
//   out instead, on the same four-major-layer, n-sub-shell scheme handed to
//   NuFast-Earth, each sub-shell at its midpoint density.  Note that the
//   Earth path in BargerPropagator::propagate uses the Y_p from the file and
//   NOT SetDensityConversion, which only reaches propagateLinear.
//
//   Matter potential.  tworttwoGf = 1.52588e-4 against this library's
//   equivalent 1.514423e-4; the ratio 0.9924922 is folded into Y_p, so that
//   the density column of the profile file stays literally PREM.
//
//   Length.  lovere = 1.26693281 L/E implies km = 1e3/1.97327e-7 eV^-1,
//   4.23e-8 longer than this library's, so the cosine is shrunk by that.
//
//   Dm2.  SetMNS takes Dm2_32, not Dm2_31 -- checked in vacuum below.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <chrono>

#include "BargerPropagator.h"
#include "our_prem.h"

static const double COSZ       = -0.9;
static const double YE         = 0.5;
static const double POT_RATIO  = 1.514423e-4/1.52588e-4;   // refined below
static const double COSZ_SCALE = 5.0677300000e9/(1.0e3/1.97327e-7);

static const double S12SQ = 0.310, S13SQ = 2.240e-2, S23SQ = 0.582;
static const double DCP   = 217.0*M_PI/180.0;
static const double DM21  = 7.39e-5, DM31 = 2.525e-3;

static std::vector<double> energies()
{
    std::vector<double> Es;
    for (int i = 0; i < 12; i++) Es.push_back(3.0*pow(40.0/3.0, i/11.0));
    return Es;
}

// Write this library's PREM as a Prob3++ profile: four major layers, each
// cut into n equal sub-shells held at their midpoint density.
static void write_profile(const char *path, int n, double yp)
{
    const double layers[4] = {1221.5, 3480.0, 5701.0, 6371.0};
    FILE *f = fopen(path, "w");
    double lo = 0.0;
    for (int L = 0; L < 4; L++)
    {
        const double hi = layers[L];
        for (int i = 0; i < n; i++)
        {
            const double outer = lo + (i+1)*(hi-lo)/n;
            const double mid   = lo + (i+0.5)*(hi-lo)/n;
            fprintf(f, "%.10f %.12f %.12f\n", outer, our_prem_rho(mid), yp);
        }
        lo = hi;
    }
    fclose(f);
}

static void run(int n, double yp, double dm_atm,
                std::vector<double> &out, double *us)
{
    char path[256];
    snprintf(path, sizeof(path), "prob3_prem_%d.dat", n);
    write_profile(path, n, yp);

    std::vector<double> Es = energies();
    BargerPropagator bnu(path);
    bnu.SetWarningSuppression(true);
    out.clear();
    for (unsigned int i = 0; i < Es.size(); i++)
    {
        bnu.SetMNS(S12SQ, S13SQ, S23SQ, DM21, dm_atm, DCP, Es[i], true, 1);
        bnu.DefinePath(COSZ*COSZ_SCALE, 0.0);
        bnu.propagate(1);
        out.push_back(bnu.GetProb(2, 2));      // 1 = e, 2 = mu, 3 = tau
    }

    if (!us) return;
    double best = 1e30;
    for (int rep = 0; rep < 5; rep++)
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        BargerPropagator b2(path);
        b2.SetWarningSuppression(true);
        for (unsigned int i = 0; i < Es.size(); i++)
        {
            b2.SetMNS(S12SQ, S13SQ, S23SQ, DM21, dm_atm, DCP, Es[i], true, 1);
            b2.DefinePath(COSZ*COSZ_SCALE, 0.0);
            b2.propagate(1);
            volatile double sink = b2.GetProb(2, 2);
            (void)sink;
        }
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
    std::vector<double> p;

    if (strcmp(mode, "--vacuum") == 0)
    {
        // Y_p = 0 is vacuum.  Print both Dm2 conventions so the right one
        // is chosen by measurement rather than by reading the header.
        std::vector<double> Es = energies();
        std::vector<double> p31, p32;
        run(16, 0.0, DM31, p31, nullptr);
        run(16, 0.0, DM31-DM21, p32, nullptr);
        for (unsigned int i = 0; i < Es.size(); i++)
            printf("%.15g %.15g %.15g\n", Es[i], p31[i], p32[i]);
        return 0;
    }

    if (strcmp(mode, "--geom") == 0)
    {
        for (int n = 1; n <= 256; n *= 4)
        {
            char path[256];
            snprintf(path, sizeof(path), "prob3_prem_%d.dat", n);
            write_profile(path, n, YE);
            BargerPropagator b(path);
            b.SetWarningSuppression(true);
            b.SetMNS(S12SQ, S13SQ, S23SQ, DM21, DM31-DM21, DCP, 10.0, true, 1);
            b.DefinePath(COSZ*COSZ_SCALE, 0.0);
            b.propagate(1);
            printf("n=%3d  GetPathLength=%.9f  GetVacuumProb(2,2)=%.15g\n",
                   n, b.GetPathLength(), b.GetVacuumProb(2, 2, 10.0, b.GetPathLength()));
        }
        return 0;
    }

    if (strcmp(mode, "--dmscan") == 0)
    {
        double lo = atof(argv[2]), hi = atof(argv[3]);
        int steps = atoi(argv[4]);
        for (int k = 0; k <= steps; k++)
        {
            double dm = lo + (hi-lo)*k/steps;
            run(16, 0.0, dm, p, nullptr);
            printf("%.12g", dm);
            for (unsigned int i = 0; i < p.size(); i++)
                printf(" %.15g", p[i]);
            printf("\n");
        }
        return 0;
    }

    if (strcmp(mode, "--scan") == 0)
    {
        double lo = atof(argv[2]), hi = atof(argv[3]);
        int steps = atoi(argv[4]);
        for (int k = 0; k <= steps; k++)
        {
            double extra = lo + (hi-lo)*k/steps;
            run(256, YE*POT_RATIO*extra, DM31-DM21, p, nullptr);
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
        double us = 0.0;
        run(ns[idx], YE*POT_RATIO, DM31-DM21, p, &us);
        printf("%d %.6f", ns[idx], us);
        for (unsigned int i = 0; i < p.size(); i++) printf(" %.15g", p[i]);
        printf("\n");
    }
    return 0;
}
