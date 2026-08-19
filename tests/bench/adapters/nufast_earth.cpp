// NuFast-Earth, driven the way its author drives it.
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

constexpr double kLayers[4] = {1221.5, 3480.0, 5701.0, 6371.0};

// This library's PREM at Y_e = 0.5, cut the way the paper's appendix says:
// four major layers, each into n equal sub-shells held at their midpoint
// density.  n is shells PER LAYER, so n = 256 is 1024 shells in total.
class OurPREM : public NuFast::Earth_Density {
  public:
    OurPREM(int n_per_layer, double scale) : n_(n_per_layer) {
        rho_.resize(4 * n_);
        for (int L = 0; L < 4; ++L) {
            const double lo = L ? kLayers[L - 1] : 0.0, hi = kLayers[L];
            for (int i = 0; i < n_; ++i) {
                const double r = lo + (i + 0.5) * (hi - lo) / n_;
                rho_[L * n_ + i] = our_prem_rho(r) * 0.5 * scale;
            }
        }
    }

    // O(1): the shell index is arithmetic, because the sub-shells are equal
    // width within each major layer.  No scan.
    double rhoYe(double r) override {
        if (r >= kLayers[3]) return rho_.back();
        int L = 0;
        while (L < 3 && r > kLayers[L]) ++L;
        const double lo = L ? kLayers[L - 1] : 0.0;
        int i = static_cast<int>((r - lo) * n_ / (kLayers[L] - lo));
        if (i < 0) i = 0;
        if (i >= n_) i = n_ - 1;
        return rho_[L * n_ + i];
    }

  private:
    int                 n_;
    std::vector<double> rho_;
};

OurPREM                     *g_prem   = nullptr;
NuFast::Probability_Engine  *g_engine = nullptr;
std::vector<double>          g_e, g_z;
double                       g_s12sq, g_s13sq, g_s23sq, g_dm21, g_dm31;

}  // namespace

namespace driver {

const char *name() { return "NuFast-Earth"; }

bench::Capabilities capabilities() {
    bench::Capabilities c;
    c.batches_energy = true;
    c.batches_zenith = true;   // Set_Spectra takes both vectors at once
    c.batch_symbol   = "Probability_Engine::Set_Spectra + Get_Probabilities";
    c.knob_name      = "eigenvalue_precision";
    c.knob_domain    = {-1, 0, 1, 2};   // -1 is the exact-eigenvalue mode
    return c;
}

void setup(const bench::Problem &p) {
    // n_shells arrives through the knob's sibling; the runner passes it as the
    // shell count when the sweep is over shells rather than precision.
    const int n_per_layer = 256;
    const double scale = OUR_PREM_MASS_DEFECT_NUFAST_EARTH;  // from conversions.h

    delete g_engine; delete g_prem;
    g_prem   = new OurPREM(n_per_layer, scale);
    g_engine = new NuFast::Probability_Engine();

    g_s12sq = p.s12sq; g_s13sq = p.s13sq; g_s23sq = p.s23sq;
    g_dm21  = p.dm21;  g_dm31  = p.dm31;

    g_engine->Set_Oscillation_Parameters(p.s12sq, p.s13sq, p.s23sq, p.dcp,
                                         p.dm21, p.dm31, true);
    g_engine->Set_Earth(0.0, g_prem);
    g_engine->Set_Production_Height(0.0);
    g_engine->Set_Eigenvalue_Precision(p.knob);

    g_e = p.energies_gev;
    g_z = p.costhz.empty() ? std::vector<double>{-0.9} : p.costhz;
    g_engine->Set_Spectra(g_e, g_z);      // batched over BOTH axes, once
}

// The one thing a fit moves.  Everything else stays cached.
void configure(double dcp) {
    g_engine->Set_delta(dcp);
}

double evaluate() {
    auto probs = g_engine->Get_Probabilities();   // whole grid, one call
    double sink = 0.0;
    for (const auto &per_e : probs)
        for (const auto &m : per_e) sink += m.arr[1][1];
    return sink;
}

void teardown() { delete g_engine; delete g_prem; g_engine = nullptr; g_prem = nullptr; }

}  // namespace driver
