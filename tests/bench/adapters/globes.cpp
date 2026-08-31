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
// internally, so the density handed over carries only the ratio of our Y_e
// to its GLB_Ne_MANTLE -- a shared physical input, not a convention.  Its
// GLB_V_FACTOR is absorbed into this code's own reference instead.
//
// Length.  Handed through unscaled.  GLoBES's km is built on its own
// GLB_EV_TO_KM_FACTOR, read by conversions.hbar_c('GLoBES') and absorbed
// into its reference, so the geometry it is given is the honest one.
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

// Chord decomposition of the radial shells for one zenith cosine, unscaled:
// r(s)^2 = R^2 + s^2 - 2 s R |cos|.
// Every PREM boundary, not four of them.
//
// The four-layer cut {1221.5, 3480, 5701, R} leaves SIX PREM
// discontinuities -- 5771, 5971, 6151, 6346.6, 6356, 6368 -- inside
// its outermost layer, so every sub-shell spanning one of them
// averages across a density jump.  How much that costs depends on
// where the shell edges happen to fall relative to the jump, which
// is why the error did not fall monotonically: 9.6e-5 at 128 shells,
// 3.2e-6 at 256, then WORSE at 1.1e-5 at 512.  This library cuts at
// all sixteen crossings the chord makes, so it converged cleanly to
// 5.5e-8 while three external codes sat at 2e-6 agreeing with each
// other to four figures -- which was the tell: they were not meeting
// their own limits, they were both reproducing our profile.
//
// The paper's appendix describes the four-layer cut.  It has to be
// corrected with this.
static const int kNLayer = OUR_PREM_NB + 1;
static double layer_edge(int i) {
    return i < OUR_PREM_NB ? OUR_PREM_B[i] : OUR_EARTH_RADIUS;
}

Profile chord(double cz, int n, double rho_factor) {
    // Subdivided along the PATH, not in radius.
    //
    // GLoBES takes a chord profile -- a list of lengths and densities -- so
    // nothing here requires a radial shell stack.  Building one anyway was
    // this adapter's choice and it cost GLoBES dearly: near the turning
    // point the chord runs tangent to the shells, so a thin radial shell
    // spans a long stretch of path whose density its radial midpoint does
    // not represent.  The resulting error does not converge, it OSCILLATES
    // with where the shell edges fall relative to r_min -- measured at
    // n = 100..180, worst at 1.2e-4 and best at 9.2e-6 with no monotone
    // trend between them, and 512 shells worse than 256.
    //
    // Cutting the chord at every PREM crossing and then subdividing each
    // piece uniformly in path length is what this library does for itself,
    // and it is what a careful user of GLoBES would hand it.
    const double R = OUR_EARTH_RADIUS, acz = std::fabs(cz);
    const double rmin = R * std::sqrt(1.0 - cz * cz), total = 2.0 * R * acz;

    // Every PREM boundary the chord crosses, as positions along the path,
    // plus the two ends and the turning point itself.
    std::vector<double> cut = {0.0, acz * R, total};
    for (int i = 0; i < OUR_PREM_NB; ++i) {
        const double r = OUR_PREM_B[i];
        if (r <= rmin || r >= R) continue;
        const double d = std::sqrt(r * r - rmin * rmin);
        cut.push_back(acz * R - d);
        cut.push_back(acz * R + d);
    }
    std::sort(cut.begin(), cut.end());

    Profile prof;
    for (std::size_t k = 0; k + 1 < cut.size(); ++k) {
        const double a = cut[k], b = cut[k + 1];
        if (b - a <= 1.0e-12) continue;
        for (int i = 0; i < n; ++i) {
            const double s0 = a + i * (b - a) / n;
            const double s1 = a + (i + 1) * (b - a) / n;
            const double s = 0.5 * (s0 + s1);
            const double r = std::sqrt(R * R + s * s - 2.0 * s * R * acz);
            prof.length.push_back(s1 - s0);
            prof.density.push_back(our_prem_rho(r) * rho_factor);
        }
    }
    return prof;
}

// Which parameter the scan turns.  Objection Earth-2 named Dmsq31 as a
// realistic thing for a fit to move, and it is not delta_CP's equal:
// codes that cache do not invalidate the same things for both.
bool g_scan_dmsq31 = false;

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
    g_scan_dmsq31 = (p.scan == "dmsq31");
    static char argv0[] = "bench_globes";
    static bool inited = false;
    if (!inited) { glbInit(argv0); inited = true; }

    const int n = p.knob >= 1 ? p.knob : p.n_layers;  // Earth-1: sweepable
    // GLoBES multiplies by GLB_Ne_MANTLE itself, so the handed density
    // carries the mass defect and our Y_e over its electron fraction.
    // Only the Y_e ratio remains: GLoBES multiplies by its own
    // GLB_Ne_MANTLE, so dividing by it hands GLoBES OUR electron fraction,
    // which is a shared physical input rather than an absorbed convention.
    // The mass defect is gone -- it belongs to this code's reference.
    const double rho_factor = p.ye / GLOBES_NE_MANTLE;

    g_e = p.energies_gev;
    g_profiles.clear();
    if (!p.costhz.empty()) {
        for (double cz : p.costhz)
            g_profiles.push_back(chord(cz, n, rho_factor));
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
void configure(double v) {
    glbSetOscParams(g_params, v, g_scan_dmsq31 ? GLB_DM_31 : GLB_DELTA_CP);
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

// Nothing to reset: glb_probability_matrix rebuilds every S-matrix from the
// profile on each call, so each repetition is already cold.
void reset() {}

// Untimed.  Mirrors evaluate()'s traversal exactly -- one profile per cosz
// outer, energy inner -- so both describe the same walk of the grid.
void probabilities(std::vector<double> &out) {
    double P[3][3];
    out.reserve(out.size() + 3*g_profiles.size()*g_e.size());
    for (const Profile &prof : g_profiles) {
        const int psteps = static_cast<int>(prof.length.size());
        for (double e : g_e) {
            glb_probability_matrix(P, +1, e, psteps, prof.length.data(),
                                   prof.density.data(), -1.0, nullptr);
            for (int b = 0; b < 3; ++b) out.push_back(P[1][b]);  // nu_mu row
        }
    }
}

}  // namespace driver
