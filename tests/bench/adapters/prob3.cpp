// Prob3++, driven the way its interface is shaped: one (energy, zenith) at
// a time.  SetMNS takes the energy as an argument and propagate() takes no
// grid, so evaluate() loops -- in C++, inside this translation unit -- and
// capabilities() says batches_energy = false.  The loop is the code's own
// interface, not a handicap imposed on it.
//
// Everything invariant under the delta_CP scan is hoisted into setup(): the
// profile file (this library's PREM cut into the same four-major-layer,
// n-sub-shell scheme handed to NuFast-Earth), the BargerPropagator built
// from it, and the trajectory (DefinePath only stores the cosine and path
// length; propagate() re-derives the per-trajectory profile itself on every
// call, which is Prob3++'s own design and is timed as such).
//
// Conventions, all from conversions.h and the Problem -- no typed constants:
//   * SetMNS takes Dmsq32, not Dmsq31: OSC_DMSQ32, derived in conversions.py.
//   * On the sphere path the electron fraction comes from the profile file's
//     Y_p column and NOT SetDensityConversion (which only reaches
//     propagateLinear), so the mass defect OUR_PREM_MASS_DEFECT_PROB3 is
//     folded into Y_p and the density column stays literally PREM.
//   * Prob3++'s km is built on hbar c = 1.97327e-7 eV m; the chord is
//     L = -2 R cos th_z, so the cosine carries OUR_COSZ_HBARC_SCALE.
#include "../bench.hpp"
#include "BargerPropagator.h"
#include "our_prem.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

BargerPropagator   *g_prop = nullptr;
std::vector<double> g_e, g_z;
double g_s12sq, g_s13sq, g_s23sq, g_dm21;
double g_L, g_rho, g_dcp;
bool   g_sphere = false;    // chord through the Earth vs constant density

// This library's PREM as a Prob3++ radial profile: the four major layers
// (radii from our_prem.h, not typed), each cut into n equal sub-shells held
// at their midpoint density.  Rows are "outer_radius rho Y_p"; Y_p carries
// the electron fraction times the mass defect.
std::string write_profile(int n, double yp) {
    const double layers[4] = {OUR_PREM_B[0], OUR_PREM_B[1], OUR_PREM_B[2],
                              OUR_EARTH_RADIUS};
    const char *tmp = std::getenv("TMPDIR");
    char path[512];
    std::snprintf(path, sizeof path, "%s/bench_prob3_prem_%d.dat",
                  tmp && tmp[0] ? tmp : "/tmp", n);
    FILE *f = std::fopen(path, "w");
    double lo = 0.0;
    for (int L = 0; L < 4; ++L) {
        const double hi = layers[L];
        for (int i = 0; i < n; ++i) {
            const double outer = lo + (i + 1) * (hi - lo) / n;
            const double mid   = lo + (i + 0.5) * (hi - lo) / n;
            std::fprintf(f, "%.10f %.12f %.12f\n", outer, our_prem_rho(mid),
                         yp);
        }
        lo = hi;
    }
    std::fclose(f);
    return path;
}

}  // namespace

namespace driver {

const char *name() { return "Prob3++"; }

bench::Capabilities capabilities() {
    bench::Capabilities c;
    c.batches_energy = false;   // SetMNS/propagate take one energy; the loop
    c.batches_zenith = false;   // in evaluate() is the interface, in C++,
    c.batch_symbol   = "none";  // not a handicap imposed from outside.
    c.knob_name      = "n_shells";      // sub-shells per major PREM layer,
    c.knob_domain    = {1, 2, 4, 8, 16, 32, 64, 128, 256};  // 256 -> 1024 total
    return c;
}

void setup(const bench::Problem &p) {
    const int n = p.knob >= 1 ? p.knob : 256;    // default: the dense grid

    g_e      = p.energies_gev;
    g_sphere = !p.costhz.empty();
    g_s12sq = p.s12sq; g_s13sq = p.s13sq; g_s23sq = p.s23sq;
    g_dm21  = p.dm21;
    g_dcp   = p.dcp;
    g_L     = p.L_km;
    g_rho   = p.density;

    delete g_prop;
    if (g_sphere) {
        g_z.clear();
        for (double cz : p.costhz) g_z.push_back(cz * OUR_COSZ_HBARC_SCALE);
        const std::string path =
            write_profile(n, p.ye * OUR_PREM_MASS_DEFECT_PROB3);
        // EarthDensity announces every profile load on stdout, which would
        // sit in front of the harness's JSON; silence cout for the one
        // untimed constructor call, then remove the file it has now read.
        std::ostringstream sink;
        std::streambuf *old = std::cout.rdbuf(sink.rdbuf());
        g_prop = new BargerPropagator(path.c_str());
        std::cout.rdbuf(old);
        std::remove(path.c_str());
        // The trajectory is invariant under the scan; propagate() rebuilds
        // the per-trajectory profile from it on every call by itself.
        g_prop->DefinePath(g_z[0], 0.0);
    } else {
        g_prop = new BargerPropagator();
        // On the linear path the electron fraction enters through
        // SetDensityConversion, and the mass defect rides with it.
        g_prop->SetDensityConversion(p.ye * OUR_PREM_MASS_DEFECT_PROB3);
    }
    g_prop->SetWarningSuppression(true);
}

// The one thing a fit moves.  Prob3++ takes delta through SetMNS, which
// also takes the energy, so the actual call sits in evaluate()'s loop.
void configure(double dcp) { g_dcp = dcp; }

double evaluate() {
    double sink = 0.0;
    if (g_sphere) {
        for (std::size_t j = 0; j < g_z.size(); ++j) {
            if (g_z.size() > 1) g_prop->DefinePath(g_z[j], 0.0, false);
            for (double e : g_e) {
                g_prop->SetMNS(g_s12sq, g_s13sq, g_s23sq, g_dm21, OSC_DMSQ32,
                               g_dcp, e, true, 1);
                g_prop->propagate(1);
                sink += g_prop->GetProb(2, 2);      // 1 = e, 2 = mu, 3 = tau
            }
        }
        return sink;
    }
    for (double e : g_e) {
        g_prop->SetMNS(g_s12sq, g_s13sq, g_s23sq, g_dm21, OSC_DMSQ32,
                       g_dcp, e, true, 1);
        g_prop->propagateLinear(1, g_L, g_rho);
        sink += g_prop->GetProb(2, 2);
    }
    return sink;
}

}  // namespace driver
