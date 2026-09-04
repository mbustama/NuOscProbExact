// NuFast-Earth, driven the way its author drives it, in whichever of its two
// modes the problem calls for: the Earth model for a chord, and
// single_trajectory_mode -- Set_E_Spectra + Set_Trajectory, added in v1.1.0 --
// for constant density.
//
// Everything invariant under the scan is built in setup() and therefore never
// timed: the Earth model, the engine, the profile, the production height, the
// eigenvalue precision, and the energy-and-zenith spectra.  configure() moves
// delta_CP alone.  evaluate() makes ONE call, because Set_Spectra batches over
// energy *and* zenith together and Get_Probabilities returns the whole grid.
//
// The previous driver put the PREM object, the engine and all five setters
// inside its own timed loop, per repetition, over twelve energies -- which is
// what the author's email objected to, and it inflated his point by more than
// the caching it defeated.
//
// The density lookup below is O(1) by arithmetic.  Its predecessor scanned a
// vector of up to 1024 discontinuities on every call, from inside the engine's
// inner loop, so our harness contributed O(N^2) work to his measurement.
#include "../bench.hpp"
#include "Earth.h"
#include "Geometry.h"
#include "NuFastEarth.h"
#include "our_prem.h"
#include "conversions.h"

#include <cassert>
#include <vector>

namespace {

// The four major PREM layers, from our_prem.h rather than typed: prob3.cpp
// and globes.cpp cut the same boundaries, and a second copy is how "the same
// Earth for every code" stops being true.
// Four major layers, as the paper's appendix states; see prob3.cpp for why
// the ten-boundary cut was tried and reverted.
const int kNLayer = 4;
static double kLayer(int i) {
    return i < 3 ? OUR_PREM_B[i] : OUR_EARTH_RADIUS;
}

// This library's PREM at Y_e = 0.5, cut the way the paper's appendix says:
// four major layers, each into n equal sub-shells held at their midpoint
// density.  n is shells PER LAYER, so n = 256 is 1024 shells in total.
class OurPREM : public NuFast::Earth_Density {
  public:
    OurPREM(int n_per_layer, double ye, bool mean)
        : n_(n_per_layer), ye_(ye), mean_(mean) {
        rho_.resize(kNLayer * n_);
        // Earth_Density's contract: the engine reads these three fields
        // directly and validates none of them.  n_discontinuities and
        // constant_shells have no default member initialiser, so leaving
        // them unset is not "the default" -- it is indeterminate, and
        // Geometry.cpp then indexes an empty `discontinuities` under a
        // garbage bound.  Every upstream subclass sets all three.
        // The discontinuity list IS the discretisation: Mean_Density makes
        // ONE slab per interval between consecutive entries, whatever
        // constant_shells says.  Declaring only the nine PREM boundaries
        // therefore hands the engine ten slabs for the whole chord and its
        // error against continuous PREM was 3e-2 -- five orders WORSE than
        // the stepped configuration.  The fine sub-shell cut stays; the flag
        // below only chooses what density each thin slab carries.
        discontinuities.resize(kNLayer * n_);
        n_discontinuities = kNLayer * n_;
        // Each sub-shell is held at one density, which is what puts the
        // engine on its cached path: one eigendecomposition per (energy,
        // shell), reused across every zenith angle.  That reuse is the
        // advantage objection Earth-3 exists to measure.
        // FALSE: the density varies continuously inside each PREM shell and
        // this engine can integrate it -- Mean_Density samples rhoYe(r) along
        // the path when constant_shells is false.  Every code is judged
        // against the continuous PREM, so handing this one a pre-stepped
        // profile would score it against an Earth it never solved.
        //
        // NOT NuFast::PREM_Full, which would have been the obvious route: its
        // polynomial is identical to ours coefficient for coefficient, but its
        // PREM_Full_Ye returns 0.466 below 3480 km and 0.494 above, where
        // every other code in this comparison is given 0.5.  Using it would
        // hand this code a different matter potential from every other.  (It
        // also carries the reserve-then-index defect: discontinuities.reserve
        // followed by discontinuities[i] = , writing past size() == 0.)
        //
        // This costs the engine its eigenvalue caching -- constant_shells
        // false takes the eigens_varying branch, one decomposition per
        // (energy, cosz, layer) rather than per (energy, shell) reused across
        // zenith.  That caching is what objection Earth-3 measures, so the
        // SPEED axis must use the stepped configuration and the two must
        // never share an axes.
        constant_shells = !mean_;
        for (int L = 0; L < kNLayer; ++L) {
            const double lo = L ? kLayer(L - 1) : 0.0, hi = kLayer(L);
            for (int i = 0; i < n_; ++i) {
                const double r = lo + (i + 0.5) * (hi - lo) / n_;
                rho_[L * n_ + i] = our_prem_rho(r) * ye;
                discontinuities[L * n_ + i] = lo + (i + 1) * (hi - lo) / n_;
                // Shell k is the region BELOW discontinuities[k]: the
                // convention Calculate_Eigens (which samples
                // rhoYe(discontinuities[j] - 1e-8)) and
                // Calculate_Internal_Amplitudes (which indexes with
                // i_discontinuity + 1) already agree on.  Ascending, with
                // the surface last, as Mean_Density's downward scan needs.

            }
        }
    }

    // O(1): the shell index is arithmetic, because the sub-shells are equal
    // width within each major layer.  No scan.
    double rhoYe(double r) override {
        // Midpoint: the engine reads the stepped shell value and caches one
        // eigendecomposition per shell, reused across every zenith angle.
        // Mean: it samples this along the path and averages per slab, which
        // defeats that caching.  The same Earth either way.
        if (!mean_) return rhoYe_stepped(r);
        if (r >= OUR_EARTH_RADIUS) return our_prem_rho(OUR_EARTH_RADIUS)*ye_;
        return our_prem_rho(r)*ye_;
    }

    double rhoYe_stepped(double r) {
        if (r >= kLayer(kNLayer - 1)) return rho_.back();
        int L = 0;
        while (L < kNLayer - 1 && r > kLayer(L)) ++L;
        const double lo = L ? kLayer(L - 1) : 0.0;
        int i = static_cast<int>((r - lo) * n_ / (kLayer(L) - lo));
        if (i < 0) i = 0;
        if (i >= n_) i = n_ - 1;
        return rho_[L * n_ + i];
    }

  private:
    int                 n_;
    double              ye_ = 0.5;
    bool                mean_ = false;
    std::vector<double> rho_;
};

OurPREM                     *g_prem   = nullptr;
NuFast::Probability_Engine  *g_engine = nullptr;
std::vector<double>          g_e, g_z;
double                       g_s12sq, g_s13sq, g_s23sq, g_dm21, g_dm31;
double                       g_delta = 0.0;
int                          g_knob  = 0;
bool                         g_constant = false;   // single-trajectory mode

// Which parameter the scan turns.  Objection Earth-2 named Dmsq31 as a
// realistic thing for a fit to move, and it is not delta_CP's equal:
// codes that cache do not invalidate the same things for both.
bool g_scan_dmsq31 = false;

}  // namespace

namespace driver {

const char *name() { return "NuFast-Earth"; }

bench::Capabilities capabilities() {
    bench::Capabilities c;
    c.batches_energy = true;
    c.batches_zenith = true;   // Set_Spectra takes both vectors at once
    c.batch_symbol   = "Probability_Engine::Set_Spectra + Get_Probabilities";
    c.knob_name      = "eigenvalue_precision";
    c.knob_domain    = {-1, 0, 1, 2, 3};   // -1 is the exact-eigenvalue mode
    return c;
}

void setup(const bench::Problem &p) {
    g_scan_dmsq31 = (p.scan == "dmsq31");
    const int n_per_layer = p.n_layers;   // Earth-1: sweepable
    delete g_engine; delete g_prem;
    g_prem   = nullptr;
    g_engine = new NuFast::Probability_Engine();

    g_s12sq = p.s12sq; g_s13sq = p.s13sq; g_s23sq = p.s23sq;
    g_dm21  = p.dm21;  g_dm31  = p.dm31;
    g_delta = p.dcp;
    g_knob  = p.knob;
    g_e     = p.energies_gev;

    g_engine->Set_Oscillation_Parameters(p.s12sq, p.s13sq, p.s23sq, p.dcp,
                                         p.dm21, p.dm31, true);
    g_engine->Set_Eigenvalue_Precision(p.knob);

    // The engine's two modes are mutually exclusive -- Set_Earth and
    // Set_Spectra assert `not single_trajectory_mode`, Set_E_Spectra and
    // Set_Trajectory assert `not earth_mode` -- so the choice is made once,
    // here, and never switched.
    g_constant = p.costhz.empty();
    if (g_constant) {
        // Objection Earth-4.  This code has had a constant-density mode since
        // v1.1.0 and the previous adapter never used it: it substituted a
        // chord at cosz = -0.9 and stamped the artifact CONST, which is a
        // wrong answer that looks right.  One slab, honest rhoYe.
        // Set_Production_Height must NOT be called here -- it asserts
        // `not single_trajectory_mode`, and the constructor already defaults
        // production_height to zero.
        g_z.clear();
        g_engine->Set_E_Spectra(g_e);
        g_engine->Set_Trajectory({{p.L_km, p.density*p.ye}});
    } else {
        // Honest rhoYe: no mass-defect factor.  This code's YerhoE2a is
        // absorbed into its own reference instead of being applied here.
        g_prem = new OurPREM(n_per_layer, p.ye, p.mean_density);
        g_engine->Set_Earth(0.0, g_prem);
        g_engine->Set_Production_Height(0.0);
        g_z = p.costhz;
        g_engine->Set_Spectra(g_e, g_z);  // batched over BOTH axes, once
    }
}

// The one thing a fit moves.  Everything else stays cached.
void configure(double v) {
    // Set_Dmsq31 is the setter the objection itself named.  Going through
    // it rather than rebuilding the engine is the whole point: what is
    // being measured is how much this code has to recompute.
    if (g_scan_dmsq31) { g_engine->Set_Dmsq31(v); return; }
    g_delta = v;
    g_engine->Set_delta(v);
}

// The only adapter with anything to reset.  Set_Eigenvalue_Precision clears
// eigens, internal amplitudes and probabilities together and -- the reason it
// is the right lever rather than a convenient one -- leaves
// trajectories_calculated alone, so the chord geometry and Mean_Densities stay
// hoisted exactly as every other code's profile does.  The engine's own
// invalidation graph draws the setup/request line in the same place bench.hpp
// asks for.  The value handed back is the one already in force, so this
// changes no physics; it only marks the cached stages stale.
void reset() { g_engine->Set_Eigenvalue_Precision(g_knob); }

double evaluate() {
    auto probs = g_engine->Get_Probabilities();   // whole grid, one call
    double sink = 0.0;
    for (const auto &per_e : probs)
        for (const auto &m : per_e) sink += m.arr[1][1];
    return sink;
}

// Untimed.  Get_Probabilities() is indexed [energy][cosz]; grid order here is
// cosz outer, energy inner, which is what the four Python-side adapters and
// the two looping C++ ones already produce.  Transposing here rather than
// re-ordering six adapters keeps one definition of grid order.
void probabilities(std::vector<double> &out) {
    auto probs = g_engine->Get_Probabilities();
    const std::size_t n_e = probs.size();
    const std::size_t n_z = n_e ? probs[0].size() : 0;
    out.reserve(out.size() + 3 * n_e * n_z);
    for (std::size_t j = 0; j < n_z; ++j)
        for (std::size_t i = 0; i < n_e; ++i)
        {                                           // the nu_mu row
            out.push_back(probs[i][j].arr[1][0]);
            out.push_back(probs[i][j].arr[1][1]);
            out.push_back(probs[i][j].arr[1][2]);
        }
}

void teardown() { delete g_engine; delete g_prem; g_engine = nullptr; g_prem = nullptr; }

}  // namespace driver
