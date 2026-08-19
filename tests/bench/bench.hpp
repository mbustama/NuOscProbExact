// The benchmark harness for the compiled codes.  It owns main(), every clock,
// and the statistics; a driver supplies physics and nothing else.
//
// That division is the point rather than a convenience.  The previous
// generation of drivers each wrote its own timing loop, and one of them ended
// up with the engine construction, the Earth object and five setters inside
// the timed region, amortized over twelve energies -- which measured our
// harness rather than the code.  Here a driver cannot time anything: it has no
// clock, and if it defines main() it will not link.
//
// A driver defines, in namespace `driver`:
//
//   const char*   name();
//   Capabilities  capabilities();
//   void          setup(const Problem&);   // hoisted; never timed
//   void          configure(double dcp);   // timed, once per scan step
//   double        evaluate();              // timed; returns a checksum
//   void          probabilities(std::vector<double>&);  // untimed; the answer
//
// `probabilities` exists because a harness that can only return a checksum can
// only measure speed, and the accuracy axis is where most of the disputed
// claims live.  It is called outside every clock, once, and must append the
// scored channel for each grid point in grid order.  A driver that does not
// implement it will not link -- which is deliberate: an adapter that cannot be
// checked for accuracy should not be usable for speed either.
//
// `setup` receives the whole problem and must do everything that does not
// change as the scan parameter moves: build the body or Earth model, construct
// the engine, install the profile, allocate the grid, and call every setter
// except the one `configure` moves.  `evaluate` must consume the WHOLE grid in
// one call wherever the code offers a batched entry point -- what
// capabilities() promises is checked against tests/bench/manifest.json.
//
// `evaluate` returns a checksum so that no optimizer can delete the work.
#pragma once

#include "conversions.h"   // generated; carries the shared parameter set

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace bench {

struct Capabilities {
    bool        batches_energy   = false;
    bool        batches_zenith   = false;
    const char *batch_symbol     = "";     // the exact entry point, for the registry test
    const char *knob_name        = "";     // the precision dial this driver exposes
    std::vector<int> knob_domain = {};     // every value, endpoints included
};

struct Problem {
    std::vector<double> energies_gev;
    std::vector<double> costhz;            // empty for constant-density problems
    double L_km    = 1300.0;
    double density = 3.0;                  // g/cm^3
    double ye      = 0.5;
    int    knob    = 0;                    // the precision setting for this run
    // The shared parameter set, from conversions.h, which is generated from
    // tests/bench/manifest.json.  Not typed here: every code must be handed
    // the same numbers, and a second copy is how that stops being true.
    double s12sq = OSC_S12SQ, s13sq = OSC_S13SQ, s23sq = OSC_S23SQ;
    double dcp   = OSC_DCP_RAD;
    double dm21  = OSC_DMSQ21, dm31 = OSC_DMSQ31;

    std::size_t points() const {
        return energies_gev.size() * (costhz.empty() ? 1 : costhz.size());
    }
};

struct Stats {
    double mean = 0.0, sd = 0.0, min = 0.0;
    int    n    = 0;
};

}  // namespace bench

namespace driver {
const char          *name();
bench::Capabilities  capabilities();
void                 setup(const bench::Problem &);
void                 configure(double dcp);
double               evaluate();
void                 probabilities(std::vector<double> &out);
}  // namespace driver

namespace bench {
namespace detail {

inline Stats reduce(const std::vector<double> &v) {
    Stats s;
    s.n = static_cast<int>(v.size());
    if (s.n == 0) return s;
    double sum = 0.0, sumsq = 0.0;
    s.min = v[0];
    for (double x : v) { sum += x; sumsq += x * x; if (x < s.min) s.min = x; }
    s.mean = sum / s.n;
    // Sample standard deviation, guarded: the variance can come out very
    // slightly negative from cancellation when every sample is equal.
    double var = (sumsq / s.n) - s.mean * s.mean;
    s.sd = var > 0.0 ? std::sqrt(var * s.n / (s.n - 1 > 0 ? s.n - 1 : 1)) : 0.0;
    return s;
}

using clk = std::chrono::high_resolution_clock;
inline double seconds(clk::time_point a, clk::time_point b) {
    return std::chrono::duration<double>(b - a).count();
}

// AMORTIZED: the scan is the timed region.  Everything invariant under the
// scan was hoisted into setup(), which ran before the clock started.  Adopted
// from the NuFast-Earth author's own Atmospheric_Speed, so the definition of
// "fair" here is his rather than ours.
inline Stats amortized(const Problem &p, int samples, int steps, double *sink) {
    const double d0 = p.dcp, dd = 0.2 / steps;
    for (int k = 0; k < steps; ++k) {           // untimed warm-up pass
        driver::configure(d0 + k * dd);
        *sink += driver::evaluate();
    }
    std::vector<double> per_point;
    per_point.reserve(samples);
    for (int s = 0; s < samples; ++s) {
        auto t0 = clk::now();
        for (int k = 0; k < steps; ++k) {
            driver::configure(d0 + k * dd);
            *sink += driver::evaluate();
        }
        auto t1 = clk::now();
        per_point.push_back(seconds(t0, t1) / (steps * p.points()) * 1e6);
    }
    return reduce(per_point);
}

// THROUGHPUT: one request for the whole grid, repeated.  Batched codes make
// one call; a code without a batched entry point loops inside evaluate(), in
// its own language, and says so through capabilities().
inline Stats throughput(const Problem &p, int samples, double min_block,
                        double *sink) {
    driver::configure(p.dcp);
    *sink += driver::evaluate();                // untimed warm-up
    int reps = 1;                               // autorange to a stable block
    for (;;) {
        auto t0 = clk::now();
        for (int r = 0; r < reps; ++r) *sink += driver::evaluate();
        double dt = seconds(t0, clk::now());
        if (dt >= min_block || reps > (1 << 24)) break;
        reps = dt > 0.0 ? static_cast<int>(reps * (min_block / dt) * 1.25) + 1
                        : reps * 8;
    }
    std::vector<double> per_point;
    per_point.reserve(samples);
    for (int s = 0; s < samples; ++s) {
        auto t0 = clk::now();
        for (int r = 0; r < reps; ++r) *sink += driver::evaluate();
        auto t1 = clk::now();
        per_point.push_back(seconds(t0, t1) / (double(reps) * p.points()) * 1e6);
    }
    Stats st = reduce(per_point);
    return st;
}

inline std::vector<double> logspace(double lo, double hi, int n) {
    std::vector<double> v;
    v.reserve(n);
    for (int i = 0; i < n; ++i)
        v.push_back(n == 1 ? lo : lo * std::pow(hi / lo, double(i) / (n - 1)));
    return v;
}

inline std::vector<double> linspace(double lo, double hi, int n) {
    std::vector<double> v;
    v.reserve(n);
    for (int i = 0; i < n; ++i)
        v.push_back(n == 1 ? lo : lo + (hi - lo) * double(i) / (n - 1));
    return v;
}

}  // namespace detail
}  // namespace bench

// The harness owns main().  A driver that defines its own will fail to link,
// which is the mechanism that keeps every clock in here.
int main(int argc, char **argv) {
    using namespace bench;

    std::string protocol = "amortized", grid = "CHORD/12x1", out;
    int samples = 30, steps = 25, n_e = 0, n_z = 0, knob = 0;
    double min_block = 0.05;

    for (int i = 1; i < argc; ++i) {
        auto eq = [&](const char *f) { return std::strcmp(argv[i], f) == 0; };
        if      (eq("--protocol")  && i + 1 < argc) protocol  = argv[++i];
        else if (eq("--grid")      && i + 1 < argc) grid      = argv[++i];
        else if (eq("--knob")      && i + 1 < argc) knob      = std::atoi(argv[++i]);
        else if (eq("--samples")   && i + 1 < argc) samples   = std::atoi(argv[++i]);
        else if (eq("--steps")     && i + 1 < argc) steps     = std::atoi(argv[++i]);
        else if (eq("--n-energies")&& i + 1 < argc) n_e       = std::atoi(argv[++i]);
        else if (eq("--n-zenith")  && i + 1 < argc) n_z       = std::atoi(argv[++i]);
        else if (eq("--json")      && i + 1 < argc) out       = argv[++i];
    }

    Problem p;
    p.knob = knob;
    if (grid == "CHORD/12x1") {
        p.energies_gev = detail::logspace(3.0, 40.0, n_e ? n_e : 12);
        p.costhz       = {-0.9};
    } else if (grid == "OSC/100x100") {
        p.energies_gev = detail::logspace(1.0, 100.0, n_e ? n_e : 100);
        p.costhz       = detail::linspace(-1.0, -0.1, n_z ? n_z : 100);
    } else if (grid == "CONST/60E") {
        p.energies_gev = detail::logspace(0.6, 20.0, n_e ? n_e : 60);
    } else {  // CONST/N-sweep
        p.energies_gev = detail::logspace(0.1, 31.6, n_e ? n_e : 1000);
    }

    driver::setup(p);
    Capabilities cap = driver::capabilities();

    // ACCURACY is untimed by construction: it configures once, asks for the
    // probabilities, and prints them.  Nothing here reads a clock, so a
    // number from this protocol can never be mistaken for a speed.
    if (protocol == "accuracy") {
        driver::configure(p.dcp);
        std::vector<double> probs;
        driver::probabilities(probs);
        std::printf("{\n  \"code\": \"%s\",\n"
                    "  \"protocol\": {\"name\": \"accuracy\", \"grid\": \"%s\"},\n"
                    "  \"knob\": {\"%s\": %d},\n  \"n_points\": %zu,\n"
                    "  \"probabilities\": [",
                    driver::name(), grid.c_str(),
                    cap.knob_name[0] ? cap.knob_name : "none", knob,
                    probs.size());
        for (std::size_t i = 0; i < probs.size(); ++i)
            std::printf("%s%.17g", i ? ", " : "", probs[i]);
        std::printf("]\n}\n");
        return 0;
    }

    double sink = 0.0;
    Stats st = (protocol == "throughput")
                   ? detail::throughput(p, samples, min_block, &sink)
                   : detail::amortized(p, samples, steps, &sink);

    char buf[2048];
    std::snprintf(buf, sizeof buf,
        "{\n  \"code\": \"%s\",\n  \"protocol\": {\"name\": \"%s\", \"grid\": \"%s\"},\n"
        "  \"knob\": {\"%s\": %d},\n  \"n_points\": %zu,\n"
        "  \"us_per_point\": {\"mean\": %.6g, \"sd\": %.6g, \"min\": %.6g, \"n\": %d},\n"
        "  \"batched\": {\"energy\": %s, \"zenith\": %s, \"symbol\": \"%s\"},\n"
        "  \"checksum\": %.17g\n}\n",
        driver::name(), protocol.c_str(), grid.c_str(),
        cap.knob_name[0] ? cap.knob_name : "none", knob, p.points(),
        st.mean, st.sd, st.min, st.n,
        cap.batches_energy ? "true" : "false",
        cap.batches_zenith ? "true" : "false", cap.batch_symbol, sink);

    std::fputs(buf, stdout);
    if (!out.empty()) { if (FILE *f = std::fopen(out.c_str(), "w")) { std::fputs(buf, f); std::fclose(f); } }
    return 0;
}
