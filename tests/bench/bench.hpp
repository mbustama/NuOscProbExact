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
//   void          reset();                 // timed; makes the next evaluate cold
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
// `reset` exists because THROUGHPUT repeats `evaluate` only to accumulate a
// block long enough to time: the repetition is measurement scaffold, not a
// workload, since no user asks the same question five thousand times against
// unchanged state.  For that mean to be the cost of a request, every
// repetition must cost what the first one cost.  Six of the seven codes are
// stateless enough that this is automatic and their `reset` is empty; the
// seventh, NuFast-Earth, caches its eigenvalues and internal amplitudes
// across calls, so without this its second and later repetitions timed a
// rotation of amplitudes an earlier untimed call had left lying around.
//
// `reset` must return the driver to the state a fresh request would find, and
// NOT undo what `setup` legitimately hoisted.  The line is the same for every
// code: geometry and profile installation belong to `setup`, and anything
// depending on the oscillation parameters is per request.  It is called only
// under THROUGHPUT.  AMORTIZED never calls it, because there the caching is
// exactly what is being measured -- that protocol is the author of
// NuFast-Earth's own, and defeating his caching inside it would measure
// nothing he would recognise.
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
    // Objection Earth-1 asks for the shell count to be swept and recorded, and
    // a single `knob` could not carry it: NuFast-Earth spends that one on its
    // eigenvalue precision, so the manifest promised a sweep the harness could
    // not run.  `n_layers` is the second dial -- sub-shells per major PREM
    // layer -- and `shells_total` is what the artifact must state, because
    // "256 shells" and "1024 shells" are the same configuration described two
    // ways and the objection was precisely about which.
    int    n_layers = 256;
    // WHICH parameter the amortized scan turns.  Objection Earth-2 asked for
    // a realistic loop and named two candidates, and they are not
    // interchangeable: NuFast-Earth caches by parameter, so delta_CP -- the
    // one this harness scanned -- invalidates the least of it and shows that
    // code at its best.  Scanning Dmsq31 invalidates more.  Measuring only
    // the favourable one and reporting it as "the" amortized cost would be
    // the same error as timing a batched code in a loop.
    std::string scan = "dcp";
    double scan_base()  const { return scan == "dmsq31" ? dm31 : dcp; }
    // Two per cent of Dmsq31 is far inside its measured uncertainty and far
    // outside anything a cache could treat as unchanged; 0.2 rad of delta is
    // what this scan has always used.
    double scan_width() const { return scan == "dmsq31" ? 0.02*dm31 : 0.2; }
    // Objection LBL-1 asks for a looped control beside the batched series, to
    // show what batching buys for a code that offers both.  An adapter that
    // batches honours this by calling its entry point once per point instead.
    bool   force_loop = false;
    // NuFast-Earth alone has two legitimate configurations, and which one is
    // in force changes what is measured.  Midpoint densities let its engine
    // cache eigenvalues across zenith angles -- worth roughly eightfold per
    // probability on an oscillogram, and the whole subject of one objection;
    // mean densities defeat that caching but represent a varying profile
    // differently.  This was a compile-time constant, which meant a speed run
    // could silently measure the code with its main optimisation switched
    // off.  It is a run-time choice now and every artifact records it, so the
    // two can never be mixed on one axes unnoticed.
    bool   mean_density = false;
    // The shared parameter set, from conversions.h, which is generated from
    // tests/bench/manifest.json.  Not typed here: every code must be handed
    // the same numbers, and a second copy is how that stops being true.
    double s12sq = OSC_S12SQ, s13sq = OSC_S13SQ, s23sq = OSC_S23SQ;
    double dcp   = OSC_DCP_RAD;
    double dm21  = OSC_DMSQ21, dm31 = OSC_DMSQ31;

    std::size_t points() const {
        return energies_gev.size() * (costhz.empty() ? 1 : costhz.size());
    }

    // Four major PREM layers, each cut into `n_layers` sub-shells.  Stated
    // rather than left for a reader to multiply.
    int shells_total() const { return 4*n_layers; }
};

struct Stats {
    double mean = 0.0, sd = 0.0, min = 0.0;
    int    n    = 0;
    // What the autorange settled on.  Reported because a timing whose block
    // was too short to measure is not a slow code, and the only way a reader
    // can tell those apart is to be told how long the timed region was.
    int    reps        = 0;      // scan steps (amortized) or grid calls (throughput)
    double block_s     = 0.0;    // seconds in one timed block
};

}  // namespace bench

namespace driver {
const char          *name();
bench::Capabilities  capabilities();
void                 setup(const bench::Problem &);
void                 configure(double dcp);
double               evaluate();
void                 reset();
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
inline int block_samples(double block, double min_block, int lo, int hi) {
    // Same wall-clock budget for every cell: hi blocks at the floor length.
    // A cell cheap enough to be autoranged gets hi; one whose single pass
    // already outlasts the budget gets lo.  The rule reads the clock and
    // nothing else -- no code is named in it, so none can be favoured by it.
    if (block <= 0.0) return hi;
    int n = static_cast<int>(double(hi) * min_block / block);
    return n < lo ? lo : (n > hi ? hi : n);
}

inline Stats amortized(const Problem &p, int samples_lo, int samples_hi,
                       int steps, double min_block, double *sink) {
    const double d0 = p.scan_base(), span = p.scan_width();

    // Autorange the scan until the timed region is long enough to measure.
    // Fixed steps handed the cheapest codes a block of 90 microseconds --- a
    // duration in which the clock's own resolution and one scheduler tick are
    // not small corrections --- and those cells came back with a 35% spread
    // that was the harness, not the code.
    //
    // The scan spans the same delta_CP interval however many steps it takes,
    // so a longer block is a FINER scan and never a repeated one.  That is the
    // point: no code is ever handed the same delta twice, so none can be made
    // to look fast by caching a value it has already seen.  It also means the
    // measured quantity is unchanged, since microseconds-per-point divides by
    // the steps actually run.
    int total = steps > 0 ? steps : 1;
    double dt = 0.0;
    // A warm-up pass of its own, thrown away before anything is measured.
    //
    // This used to be the autorange's first iteration, which quietly broke
    // any code whose FIRST call is dear: NuFast-Earth builds its cache
    // then, so that pass ran past min_block, the loop concluded one step
    // was enough and stopped, and the thirty blocks that followed were a
    // warm cache timed over half a millisecond -- the very defect this
    // autorange exists to prevent, reintroduced by measuring the pass that
    // was supposed to be discarded.
    {
        const double dd0 = 0.2 / total;
        for (int k = 0; k < total; ++k) {
            driver::configure(d0 + k * dd0);
            *sink += driver::evaluate();
        }
    }
    for (;;) {
        const double dd = span / total;
        auto t0 = clk::now();
        for (int k = 0; k < total; ++k) {
            driver::configure(d0 + k * dd);
            *sink += driver::evaluate();
        }
        dt = seconds(t0, clk::now());
        if (dt >= min_block || total > (1 << 24)) break;
        total = dt > 0.0 ? static_cast<int>(total * (min_block / dt) * 1.25) + 1
                         : total * 8;
    }

    const int samples = block_samples(dt, min_block, samples_lo, samples_hi);

    // The scan CONTINUES across blocks rather than restarting in each one.
    // That matters where one pass over the grid already outlasts the block,
    // so a block holds a single step: restarting would hand the code the
    // same delta_CP in every block, and a code that cached on it would be
    // timed doing nothing.  Advancing instead means the whole cell sweeps
    // the interval exactly once and no delta is ever evaluated twice.
    //
    // It is also what lets the step floor be one.  A floor of twenty-five
    // forced a 25-pass block on codes whose single pass takes a minute --
    // nuCraft on the oscillogram would have spent twelve hours in this
    // function -- and bought nothing, since a block that long is already
    // far past the point where the clock is the limitation.
    const double dd = span / (double(total) * samples);
    long long step = 0;
    std::vector<double> per_point;
    per_point.reserve(samples);
    for (int s = 0; s < samples; ++s) {
        auto t0 = clk::now();
        for (int k = 0; k < total; ++k, ++step) {
            driver::configure(d0 + step * dd);
            *sink += driver::evaluate();
        }
        auto t1 = clk::now();
        per_point.push_back(seconds(t0, t1) / (double(total) * p.points()) * 1e6);
    }
    Stats st = reduce(per_point);
    st.reps = total;
    st.block_s = st.mean * 1e-6 * double(total) * p.points();
    return st;
}

// THROUGHPUT: one request for the whole grid, repeated.  Batched codes make
// one call; a code without a batched entry point loops inside evaluate(), in
// its own language, and says so through capabilities().
inline Stats throughput(const Problem &p, int samples_lo, int samples_hi,
                        double min_block, double *sink) {
    driver::configure(p.scan_base());
    driver::reset();
    *sink += driver::evaluate();                // untimed warm-up
    int reps = 1;                               // autorange to a stable block
    double block = 0.0;
    for (;;) {
        auto t0 = clk::now();
        for (int r = 0; r < reps; ++r) { driver::reset(); *sink += driver::evaluate(); }
        double dt = seconds(t0, clk::now());
        block = dt;
        if (dt >= min_block || reps > (1 << 24)) break;
        reps = dt > 0.0 ? static_cast<int>(reps * (min_block / dt) * 1.25) + 1
                        : reps * 8;
    }
    const int samples = block_samples(block, min_block, samples_lo, samples_hi);
    std::vector<double> per_point;
    per_point.reserve(samples);
    for (int s = 0; s < samples; ++s) {
        auto t0 = clk::now();
        for (int r = 0; r < reps; ++r) { driver::reset(); *sink += driver::evaluate(); }
        auto t1 = clk::now();
        per_point.push_back(seconds(t0, t1) / (double(reps) * p.points()) * 1e6);
    }
    Stats st = reduce(per_point);
    st.reps = reps;
    st.block_s = st.mean * 1e-6 * double(reps) * p.points();
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

    std::string protocol = "amortized", grid = "CHORD/12x1", out, scan = "dcp";
    int samples = 30, max_samples = 100, steps = 1;
    int n_e = 0, n_z = 0, knob = 0, n_layers = 256;
    bool force_loop = false, mean_density = false;
    double min_block = 0.25;

    for (int i = 1; i < argc; ++i) {
        auto eq = [&](const char *f) { return std::strcmp(argv[i], f) == 0; };
        if      (eq("--protocol")  && i + 1 < argc) protocol  = argv[++i];
        else if (eq("--grid")      && i + 1 < argc) grid      = argv[++i];
        else if (eq("--knob")      && i + 1 < argc) knob      = std::atoi(argv[++i]);
        else if (eq("--n-layers")  && i + 1 < argc) n_layers  = std::atoi(argv[++i]);
        else if (eq("--loop"))                      force_loop = true;
        else if (eq("--mean-density"))              mean_density = true;
        else if (eq("--samples")   && i + 1 < argc) samples   = std::atoi(argv[++i]);
        else if (eq("--steps")     && i + 1 < argc) steps     = std::atoi(argv[++i]);
        else if (eq("--min-block") && i + 1 < argc) min_block = std::atof(argv[++i]);
        else if (eq("--max-samples") && i + 1 < argc) max_samples = std::atoi(argv[++i]);
        else if (eq("--scan")      && i + 1 < argc) scan      = argv[++i];
        else if (eq("--n-energies")&& i + 1 < argc) n_e       = std::atoi(argv[++i]);
        else if (eq("--n-zenith")  && i + 1 < argc) n_z       = std::atoi(argv[++i]);
        else if (eq("--json")      && i + 1 < argc) out       = argv[++i];
    }

    Problem p;
    p.scan = scan;
    p.knob = knob;
    p.n_layers = n_layers;
    p.force_loop = force_loop;
    p.mean_density = mean_density;
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
        driver::configure(p.scan_base());
        std::vector<double> probs;
        driver::probabilities(probs);
        std::printf("{\n  \"code\": \"%s\",\n"
                    "  \"protocol\": {\"name\": \"accuracy\", \"grid\": \"%s\"},\n"
                    "  \"knob\": {\"%s\": %d},\n"
                    "  \"n_layers\": %d,\n  \"shells_total\": %d,\n"
                    "  \"conventions\": \"own-reference\",\n"
                    "  \"profile_basis\": \"continuous\",\n"
                    "  \"looped\": %s,\n"
                    "  \"shell_density\": \"%s\",\n"
                    "  \"manifest_sha256\": \"%s\",\n"
                    "  \"channels\": [\"numu->nue\", \"numu->numu\", \"numu->nutau\"],\n"
                    "  \"n_points\": %zu,\n"
                    "  \"probabilities\": [",
                    driver::name(), grid.c_str(),
                    cap.knob_name[0] ? cap.knob_name : "none", knob,
                    p.n_layers, p.shells_total(),
                    p.force_loop ? "true" : "false",
                    p.mean_density ? "mean" : "midpoint", MANIFEST_SHA256,
                    p.points());   // grid points, NOT probs.size(): there are
                                   // three channels per point now, and this
                                   // field is what downstream divides by.
        for (std::size_t i = 0; i < probs.size(); ++i)
            std::printf("%s%.17g", i ? ", " : "", probs[i]);
        std::printf("]\n}\n");

        // --json writes the artifact.  The timed branch has always honoured
        // it and this one never did: it printed to stdout and returned, so an
        // orchestrator that asked for a file got a clean exit status and no
        // file.  Found by running the orchestrator on three cells before
        // running it on two hundred.
        if (!out.empty()) {
            if (FILE *f = std::fopen(out.c_str(), "w")) {
                std::fprintf(f, "{\n  \"code\": \"%s\",\n"
                    "  \"protocol\": {\"name\": \"accuracy\", \"grid\": \"%s\"},\n"
                    "  \"knob\": {\"%s\": %d},\n"
                    "  \"n_layers\": %d,\n  \"shells_total\": %d,\n"
                    "  \"conventions\": \"own-reference\",\n"
                    "  \"profile_basis\": \"continuous\",\n"
                    "  \"looped\": %s,\n"
                    "  \"shell_density\": \"%s\",\n"
                    "  \"manifest_sha256\": \"%s\",\n"
                    "  \"channels\": [\"numu->nue\", \"numu->numu\", "
                    "\"numu->nutau\"],\n"
                    "  \"n_points\": %zu,\n  \"probabilities\": [",
                    driver::name(), grid.c_str(),
                    cap.knob_name[0] ? cap.knob_name : "none", knob,
                    p.n_layers, p.shells_total(),
                    p.force_loop ? "true" : "false",
                    p.mean_density ? "mean" : "midpoint", MANIFEST_SHA256,
                    p.points());
                for (std::size_t i = 0; i < probs.size(); ++i)
                    std::fprintf(f, "%s%.17g", i ? ", " : "", probs[i]);
                std::fprintf(f, "]\n}\n");
                std::fclose(f);
            }
        }
        return 0;
    }

    double sink = 0.0;
    Stats st = (protocol == "throughput")
                   ? detail::throughput(p, samples, max_samples, min_block, &sink)
                   : detail::amortized(p, samples, max_samples, steps,
                                       min_block, &sink);

    char buf[4096];
    std::snprintf(buf, sizeof buf,
        "{\n  \"code\": \"%s\",\n  \"protocol\": {\"name\": \"%s\", \"grid\": \"%s\"},\n"
        "  \"knob\": {\"%s\": %d},\n  \"n_points\": %zu,\n"
        "  \"us_per_point\": {\"mean\": %.6g, \"sd\": %.6g, \"min\": %.6g, \"n\": %d},\n"
        "  \"n_layers\": %d,\n  \"shells_total\": %d,\n"
        "  \"conventions\": \"own-reference\",\n"
        "  \"profile_basis\": \"continuous\",\n"
        "  \"looped\": %s,\n"
        "  \"shell_density\": \"%s\",\n"
        "  \"manifest_sha256\": \"%s\",\n"
        "  \"scan\": \"%s\",\n"
        "  \"timing\": {\"min_block_s\": %.6g, \"block_reps\": %d, "
        "\"block_seconds\": %.6g, \"samples\": %d},\n"
        "  \"batched\": {\"energy\": %s, \"zenith\": %s, \"symbol\": \"%s\"},\n"
        "  \"checksum\": %.17g\n}\n",
        driver::name(), protocol.c_str(), grid.c_str(),
        cap.knob_name[0] ? cap.knob_name : "none", knob, p.points(),
        st.mean, st.sd, st.min, st.n,
        p.n_layers, p.shells_total(), p.force_loop ? "true" : "false",
        p.mean_density ? "mean" : "midpoint", MANIFEST_SHA256,
        p.scan.c_str(), min_block, st.reps, st.block_s, st.n,
        cap.batches_energy ? "true" : "false",
        cap.batches_zenith ? "true" : "false", cap.batch_symbol, sink);

    std::fputs(buf, stdout);
    if (!out.empty()) { if (FILE *f = std::fopen(out.c_str(), "w")) { std::fputs(buf, f); std::fclose(f); } }
    return 0;
}
