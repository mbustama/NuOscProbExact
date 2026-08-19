// GLoBES, driven one energy at a time because that is its interface:
// glb_probability_matrix takes an arbitrary layered profile and ONE energy,
// so evaluate() loops -- in C++, inside this translation unit -- and
// capabilities() says batches_energy = false.  The loop is the code's own
// interface, not a handicap imposed on it.
//
// Everything invariant under the delta_CP scan is hoisted into setup():
// glbInit, the parameter block, and the chord decomposition of the profile
// (per zenith), so the timed region is glbSetOscParams + the probability
// calls and nothing of ours.
//
// Profile.  GLoBES ships no Earth model at this level; it is handed the
// chord decomposition of the SAME radial shells given to NuFast-Earth and
// Prob3++: four major layers (radii from our_prem.h, not typed), each cut
// into n equal sub-shells held at their midpoint density, evaluated with
// this library's PREM.
//
// Matter potential.  GLoBES computes density * GLB_V_FACTOR * GLB_Ne_MANTLE
// internally, so the density handed over carries the mass defect
// OUR_PREM_MASS_DEFECT_GLOBES (derived in conversions.py from the pinned
// GLB_V_FACTOR) and the ratio of our Y_e to its GLB_Ne_MANTLE -- both from
// conversions.h, neither typed here.
//
// Length.  GLoBES's km is built on hbar c = 1.97327e-7 eV m; the chord is
// L = -2 R cos th_z, so the cosine carries OUR_COSZ_HBARC_SCALE, exactly as
// tests/external_drivers/globes_drv.c did.  The constant-density baseline is
// handed through unscaled, matching the frozen constant-density drivers.
#include "../bench.hpp"
#include "our_prem.h"

#include <globes/globes.h>

#include <algorithm>
#include <cmath>
#include <vector>

// Declared in the private header source/glb_probability.h and exported by
// the library; repeated here (C linkage) so the adapter needs only the
// installed headers, exactly as tests/external_drivers/globes_drv.c did.
extern "C" int glb_probability_matrix(double P[3][3], int cp_sign, double E,
                                      int psteps, const double *length,
                                      const double *density,
                                      double filter_sigma, void *user_data);

namespace {

struct Profile {
    std::vector<double> length;    // km, per traversed segment
    std::vector<double> density;   // g/cm^3, already carrying the factors
};

glb_params           g_params = nullptr;
std::vector<double>  g_e;
std::vector<Profile> g_profiles;   // one per zenith; one constant slab else

// Chord decomposition of the radial shells for one zenith cosine (already
// carrying OUR_COSZ_HBARC_SCALE): r(s)^2 = R^2 + s^2 - 2 s R |cos|.
Profile chord(double cz, int n, double rho_factor) {
    const double layers[4] = {OUR_PREM_B[0], OUR_PREM_B[1], OUR_PREM_B[2],
                              OUR_EARTH_RADIUS};
    const double R = OUR_EARTH_RADIUS, acz = std::fabs(cz);
    const double rmin = R * std::sqrt(1.0 - cz * cz), total = 2.0 * R * acz;

    std::vector<double> edge, rho;
    double lo = 0.0;
    for (int L = 0; L < 4; ++L) {
        const double hi = layers[L];
        for (int i = 0; i < n; ++i) {
            edge.push_back(lo + (i + 1) * (hi - lo) / n);
            rho.push_back(our_prem_rho(lo + (i + 0.5) * (hi - lo) / n)
                          * rho_factor);
        }
        lo = hi;
    }

    std::vector<double> cross = {0.0, total};
    for (double r : edge) {
        if (r <= rmin || r >= R) continue;
        const double d = std::sqrt(r * r - rmin * rmin);
        cross.push_back(acz * R - d);
        cross.push_back(acz * R + d);
    }
    std::sort(cross.begin(), cross.end());

    Profile prof;
    for (std::size_t i = 0; i + 1 < cross.size(); ++i) {
        const double a = cross[i], b = cross[i + 1];
        if (b - a <= 1.0e-12) continue;              // degenerate crossing
        const double s = 0.5 * (a + b);
        const double r = std::sqrt(R * R + s * s - 2.0 * s * R * acz);
        std::size_t j = 0;
        while (j + 1 < edge.size() && r > edge[j]) ++j;
        prof.length.push_back(b - a);
        prof.density.push_back(rho[j]);
    }
    return prof;
}

}  // namespace

namespace driver {

const char *name() { return "GLoBES"; }

bench::Capabilities capabilities() {
    bench::Capabilities c;
    c.batches_energy = false;   // glb_probability_matrix takes one energy;
    c.batches_zenith = false;   // the loop in evaluate() is the interface,
    c.batch_symbol   = "none";  // in C++, not a handicap imposed from outside.
    c.knob_name      = "n_shells";      // sub-shells per major PREM layer,
    c.knob_domain    = {1, 2, 4, 8, 16, 32, 64, 128, 256};  // 256 -> 1024 total
    return c;
}

void setup(const bench::Problem &p) {
    static char argv0[] = "bench_globes";
    static bool inited = false;
    if (!inited) { glbInit(argv0); inited = true; }

    const int n = p.knob >= 1 ? p.knob : 256;    // default: the dense grid
    // GLoBES multiplies by GLB_Ne_MANTLE itself, so the handed density
    // carries the mass defect and our Y_e over its electron fraction.
    const double rho_factor =
        OUR_PREM_MASS_DEFECT_GLOBES * (p.ye / GLOBES_NE_MANTLE);

    g_e = p.energies_gev;
    g_profiles.clear();
    if (!p.costhz.empty()) {
        for (double cz : p.costhz)
            g_profiles.push_back(chord(cz * OUR_COSZ_HBARC_SCALE, n,
                                       rho_factor));
    } else {
        Profile prof;
        prof.length.push_back(p.L_km);
        prof.density.push_back(p.density * rho_factor);
        g_profiles.push_back(prof);
    }

    if (!g_params) g_params = glbAllocParams();
    glbDefineParams(g_params,
                    std::asin(std::sqrt(p.s12sq)),
                    std::asin(std::sqrt(p.s13sq)),
                    std::asin(std::sqrt(p.s23sq)),
                    p.dcp, p.dm21, p.dm31);
    glbSetDensityParams(g_params, 1.0, GLB_ALL);
    glbSetOscillationParameters(g_params);
}

// The one thing a fit moves; GLoBES recomputes its mixing matrix here,
// which is the cost the AMORTIZED protocol is defined to include.
void configure(double dcp) {
    glbSetOscParams(g_params, dcp, GLB_DELTA_CP);
    glbSetOscillationParameters(g_params);
}

double evaluate() {
    double P[3][3], sink = 0.0;
    for (const Profile &prof : g_profiles) {
        const int psteps = static_cast<int>(prof.length.size());
        for (double e : g_e) {
            glb_probability_matrix(P, +1, e, psteps, prof.length.data(),
                                   prof.density.data(), -1.0, nullptr);
            sink += P[1][1];                        // P(nu_mu -> nu_mu)
        }
    }
    return sink;
}

}  // namespace driver
