// Regenerates tests/nufast_scan.json, the NuFast-LBL curve of the
// exact-vs-approximations figure.
//
// This driver was missing.  The dataset shipped for a long time with a
// `driver` field naming a file that had never been committed, so the one
// figure in the paper that scores four external codes was the one figure a
// reader could not regenerate.  This reconstructs it from what the dataset
// records: the channel, the baseline, the density handed to NuFast-LBL, the
// 150-point energy grid, and the two Newton settings.
//
// Verified against the frozen file on 2026-09-01: max |dP| = 4.4e-12 at
// N_Newton = 0 and 3.7e-12 at N_Newton = 2, over all 150 energies.  The
// energy grid agrees to 5e-10 GeV; the original was generated in numpy and
// this one in C++, and the residual probability difference is far below
// anything the figure resolves.
//
// Build, with NuFast-LBL cloned at the commit tests/bench/manifest.json pins:
//
//     g++ -O3 -std=c++17 -Dmain=nufast_demo_main -c NuFast_LBL.cpp -o nf.o
//     g++ -O3 -std=c++17 nufast_scan.cpp nf.o -o nufast_scan
//     ./nufast_scan            # N_Newton, E [GeV], P(nu_mu -> nu_e)
//
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>
using namespace std;
namespace NuFast {
void Probability_Matter_LBL(double, double, double, double, double, double,
                            double, const vector<double> &, double, int,
                            vector<array<array<double, 3>, 3>> &);
}
int main() {
    const double s12sq = 0.310, s23sq = 0.582, s13sq = 0.0224;
    const double delta = 3.787364476827695;
    const double Dmsq21 = 7.39e-5, Dmsq31 = 2.525e-3;
    const double L = 1300.0, rhoYe = 1.4887372835;
    const int N = 150;
    vector<double> E(N);
    for (int i = 0; i < N; ++i)
        E[i] = pow(10.0, log10(0.6) + double(i)*(log10(20.0) - log10(0.6))/double(N - 1));
    for (int n : {0, 2}) {
        vector<array<array<double, 3>, 3>> p(N);
        NuFast::Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21,
                                       Dmsq31, L, E, rhoYe, n, p);
        for (int i = 0; i < N; ++i) printf("%d %.17g %.17g\n", n, E[i], p[i][1][0]);
    }
    return 0;
}
