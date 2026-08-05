/* GLoBES on the same PREM chord as everything else.
 *
 *   Profile.  GLoBES ships no Earth model at this level: glb_probability_matrix
 *   takes an arbitrary layered profile, (length[], density[]) along the path.
 *   It is handed the chord decomposition of the SAME radial shells given to
 *   NuFast-Earth and Prob3++ -- four major layers, each cut into n equal
 *   sub-shells held at their midpoint density -- evaluated with this
 *   library's PREM.  So all three layered codes see the identical Earth.
 *
 *   Matter potential.  GLB_V_FACTOR = 7.63247e-14 is sqrt(2) G_F/u with u the
 *   atomic mass unit, and GLB_Ne_MANTLE = 0.5.  Against this library's
 *   3.7860575e-14 per (g/cm^3) at Y_e = 0.5 the ratio is 0.992093, the
 *   nuclear mass defect, so the density is scaled by it.
 *
 *   Length.  GLB_EV_TO_KM_FACTOR = 1.97327e-10 gives km = 5.0677302143e9
 *   eV^-1 against this library's 5.06773e9; the cosine is shrunk by 4.23e-8.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <globes/globes.h>
#include "our_prem.h"

/* Declared in the private header source/glb_probability.h, exported by the
   library; repeated here so the driver needs only the installed headers. */
int glb_probability_matrix(double P[3][3], int cp_sign, double E,
                           int psteps, const double *length,
                           const double *density,
                           double filter_sigma, void *user_data);

#define COSZ        (-0.9)
/* GLoBES applies GLB_Ne_MANTLE = 0.5 itself, so the density handed over
   carries only the mass-defect factor and no Y_e. */
#define POT_SCALE   0.992093              /* refined by the scan below */
#define COSZ_SCALE  (5.0677300000e9/(1.0/1.97327e-10))
#define NMAX        8192

static double Es[12];

static void set_energies(void)
{
    int i;
    for (i = 0; i < 12; i++) Es[i] = 3.0*pow(40.0/3.0, i/11.0);
}

/* Chord decomposition of the radial shells: four major layers, n sub-shells
   each.  r(s)^2 = R^2 + s^2 - 2 s R |cosz| along the chord. */
static int chord(int n, double scale, double *length, double *density)
{
    const double layers[4] = {1221.5, 3480.0, 5701.0, 6371.0};
    double edge[NMAX], rho[NMAX], cross[2*NMAX+2];
    int nshell = 0, ncross = 0, i, j, k, L;
    double lo = 0.0, R = OUR_EARTH_RADIUS;
    double cz = COSZ*COSZ_SCALE, acz = fabs(cz);
    double rmin = R*sqrt(1.0 - cz*cz), total = 2.0*R*acz;

    for (L = 0; L < 4; L++)
    {
        double hi = layers[L];
        for (i = 0; i < n; i++)
        {
            edge[nshell] = lo + (i+1)*(hi-lo)/n;
            rho[nshell]  = our_prem_rho(lo + (i+0.5)*(hi-lo)/n)*scale;
            nshell++;
        }
        lo = hi;
    }

    cross[ncross++] = 0.0;
    for (j = 0; j < nshell; j++)
    {
        double r = edge[j];
        if (r <= rmin || r >= R) continue;
        double d = sqrt(r*r - rmin*rmin);
        cross[ncross++] = acz*R - d;
        cross[ncross++] = acz*R + d;
    }
    cross[ncross++] = total;

    for (i = 0; i < ncross; i++)          /* insertion sort, ncross is small */
        for (j = i+1; j < ncross; j++)
            if (cross[j] < cross[i])
            { double t = cross[i]; cross[i] = cross[j]; cross[j] = t; }

    k = 0;
    for (i = 0; i < ncross-1; i++)
    {
        double a = cross[i], b = cross[i+1], s, r;
        if (b - a <= 1.0e-12) continue;
        s = 0.5*(a+b);
        r = sqrt(R*R + s*s - 2.0*s*R*acz);
        for (j = 0; j < nshell; j++) if (r <= edge[j]) break;
        if (j >= nshell) j = nshell-1;
        length[k]  = b - a;
        density[k] = rho[j];
        k++;
    }
    return k;
}

static void run(int n, double scale, double *out, double *us)
{
    double length[2*NMAX+2], density[2*NMAX+2], P[3][3];
    int psteps = chord(n, scale, length, density), i;
    struct timespec t0, t1;
    double best = 1e30;
    int rep;

    for (i = 0; i < 12; i++)
    {
        glb_probability_matrix(P, +1, Es[i], psteps, length, density,
                               -1.0, NULL);
        out[i] = P[1][1];
    }
    if (!us) return;
    for (rep = 0; rep < 5; rep++)
    {
        double t;
        clock_gettime(CLOCK_MONOTONIC, &t0);
        psteps = chord(n, scale, length, density);
        for (i = 0; i < 12; i++)
            glb_probability_matrix(P, +1, Es[i], psteps, length, density,
                                   -1.0, NULL);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        t = ((t1.tv_sec-t0.tv_sec)*1e6 + (t1.tv_nsec-t0.tv_nsec)*1e-3)/12.0;
        if (t < best) best = t;
    }
    *us = best;
}

int main(int argc, char *argv[])
{
    glb_params p;
    double out[12];
    int i, idx;
    const int ns[] = {1, 2, 4, 8, 16, 32, 64, 128, 256};
    const char *mode = (argc > 1) ? argv[1] : "";

    glbInit(argv[0]);
    p = glbAllocParams();
    set_energies();
    glbDefineParams(p, asin(sqrt(0.310)), asin(sqrt(2.240e-2)),
                    asin(sqrt(0.582)), 217.0*M_PI/180.0, 7.39e-5, 2.525e-3);
    glbSetDensityParams(p, 1.0, GLB_ALL);
    glbSetOscillationParameters(p);
    /* no experiment is loaded, so glbSetRates() would dereference nothing */

    if (strcmp(mode, "--vacuum") == 0)
    {
        run(16, 0.0, out, NULL);
        for (i = 0; i < 12; i++) printf("%.15g %.15g\n", Es[i], out[i]);
        return 0;
    }

    if (strcmp(mode, "--scan") == 0)
    {
        double lo = atof(argv[2]), hi = atof(argv[3]);
        int steps = atoi(argv[4]), k;
        for (k = 0; k <= steps; k++)
        {
            double extra = lo + (hi-lo)*k/steps;
            run(256, POT_SCALE*extra, out, NULL);
            printf("%.10f", extra);
            for (i = 0; i < 12; i++) printf(" %.15g", out[i]);
            printf("\n");
        }
        return 0;
    }

    for (idx = 0; idx < 9; idx++)
    {
        double us = 0.0;
        run(ns[idx], POT_SCALE, out, &us);
        printf("%d %.6f", ns[idx], us);
        for (i = 0; i < 12; i++) printf(" %.15g", out[i]);
        printf("\n");
    }
    glbFreeParams(p);
    return 0;
}
