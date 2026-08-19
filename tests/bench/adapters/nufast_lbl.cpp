// NuFast-LBL, driven through its one batched entry point.
//
// Probability_Matter_LBL is a free function with no persistent engine, so
// there is nothing to cache across the delta_CP scan: setup() honestly only
// builds the energy grid and the output array, and configure() stores delta.
// The artifact will therefore show its oneshot/amortized ratio near 1, which
// is a property of the code's interface and not of this driver.
//
// evaluate() makes ONE call for the whole energy vector, because that is the
// interface: the function takes const vector<double>& E and fills a vector
// of 3x3 probability matrices.  Constant density only -- there is no zenith.
//
// The released NuFast_LBL.cpp ships no header, so the prototype is declared
// here exactly as defined at NuFast_LBL.cpp:32.  The object this links
// against was compiled with -Dmain=nufast_lbl_demo_main, so its demo main()
// cannot collide with the harness's -- the same mechanism that keeps this
// file clockless.
#include "../bench.hpp"

#include <array>
#include <vector>

namespace NuFast {
void Probability_Matter_LBL(
    double s12sq, double s13sq, double s23sq, double delta, double Dmsq21,
    double Dmsq31, double L, const std::vector<double> &E, double rhoYe,
    int N_Newton,
    std::vector<std::array<std::array<double, 3>, 3>> &probs_returned);
}  // namespace NuFast

namespace {

std::vector<double> g_e;
std::vector<std::array<std::array<double, 3>, 3>> g_probs;
double g_s12sq, g_s13sq, g_s23sq, g_dm21, g_dm31;
double g_L, g_rhoYe, g_delta;
int    g_newton;

}  // namespace

namespace driver {

const char *name() { return "NuFast-LBL"; }

bench::Capabilities capabilities() {
    bench::Capabilities c;
    c.batches_energy = true;
    c.batches_zenith = false;   // constant density; there is no zenith axis
    c.batch_symbol   = "NuFast::Probability_Matter_LBL(double,double,double,"
                       "double,double,double,double,const std::vector<double>&,"
                       "double,int,std::vector<std::array<std::array<double,3>,3>>&)";
    c.knob_name      = "N_Newton";
    c.knob_domain    = {-1, 0, 1, 2, 3};   // -1 is the exact-eigenvalue mode
    return c;
}

void setup(const bench::Problem &p) {
    g_e = p.energies_gev;
    g_probs.assign(g_e.size(), {});     // "the vector should be allocated first"

    g_s12sq = p.s12sq; g_s13sq = p.s13sq; g_s23sq = p.s23sq;
    g_dm21  = p.dm21;  g_dm31  = p.dm31;
    g_delta = p.dcp;
    g_L     = p.L_km;
    // rhoYe carries the mass-defect factor between this library's matter
    // constant and the YerhoE2a this code hard-codes; from conversions.h.
    g_rhoYe  = p.density * p.ye * OUR_PREM_MASS_DEFECT_NUFAST_LBL;
    g_newton = p.knob;
}

// No engine survives between calls, so this honestly stores delta and
// nothing else; the per-call setup cost lives inside evaluate(), where the
// code itself pays it.
void configure(double dcp) { g_delta = dcp; }

double evaluate() {
    NuFast::Probability_Matter_LBL(g_s12sq, g_s13sq, g_s23sq, g_delta,
                                   g_dm21, g_dm31, g_L, g_e, g_rhoYe,
                                   g_newton, g_probs);
    double sink = 0.0;
    for (const auto &m : g_probs) sink += m[1][1];   // P(nu_mu -> nu_mu)
    return sink;
}

}  // namespace driver
