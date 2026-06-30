#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <limits>
#include <cstdint>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "params_cli.hpp"
#include "runRecorder.hpp"
#include "colony.hpp"  // for RunResult
#include "solution.hpp"  // for Solution

namespace aco 
{

// elapsedMs = vrijeme od početka TOG runa do kraja iteracije it
void RunRecorder :: recordIteration(int runIdx, int it, 
                                   double iterBestCost, 
                                   double globalBestCost, 
                                   double elapsedMs) 
{
    const size_t idx = index(runIdx, it);
    iterBest_[idx]   = iterBestCost;
    globalBest_[idx] = globalBestCost;
    elapsedMs_[idx]  = elapsedMs;
}
void RunRecorder::recordIteration2(int runIdx, int it,
                                   double iterBestPreCost,
                                   double iterBestPostCost,
                                   double globalBestCost,
                                   double elapsedMs,
                                   double weibullK)
{
    const size_t idx = index(runIdx, it);
    iterBest_[idx]     = iterBestPreCost;   // postojeći niz, tretiramo kao PRE
    iterBestPost_[idx] = iterBestPostCost;  // novi niz, POST
    
    // globalBest_ treba biti kumulativni minimum (najbolji cost do te iteracije)
    // globalBestCost koji dolazi iz colony.cpp već JEST kumulativni minimum (rr.best.cost se ažurira samo kada se nađe bolje rješenje)
    // Ali za sigurnost, provjerimo je li globalBestCost bolji od prethodnog globalBest_
    if (it == 0) {
        globalBest_[idx] = globalBestCost;
    } else {
        const size_t prevIdx = index(runIdx, it - 1);
        const double prevGlobalBest = globalBest_[prevIdx];
        // Kumulativni minimum: uzmi najbolji od prethodnog i trenutnog
        // globalBestCost bi trebao biti >= prevGlobalBest (jer je kumulativni minimum),
        // ali za sigurnost uzimamo minimum
        globalBest_[idx] = std::min(prevGlobalBest, globalBestCost);
    }
    
    elapsedMs_[idx]    = elapsedMs;
    weibullK_[idx]     = weibullK;
}

void RunRecorder::recordIterationParams(int runIdx, int it, double q0, int LSW, bool smoothingActive)
{
    const size_t idx = index(runIdx, it);
    q0_[idx] = q0;
    LSW_[idx] = LSW;
    smoothing_[idx] = smoothingActive;
}

void RunRecorder::setCacheStats(int runIdx, double evalHitRate, double reasHitRate)
{
    cacheEvalHitRate_.at(runIdx) = evalHitRate;
    cacheReasHitRate_.at(runIdx) = reasHitRate;
}

void RunRecorder :: writeIterationBestCsv(const std::string& filename, char sep) const
{
    writeCsv(filename, sep, /*useGlobal=*/false);
}

void RunRecorder :: writeGlobalBestCsv(const std::string& filename, char sep) const
{
    writeCsv(filename, sep, /*useGlobal=*/true);
}

size_t RunRecorder :: index(int runIdx, int it) const 
{
    return static_cast<size_t>(runIdx) * static_cast<size_t>(iters_) + static_cast<size_t>(it);
}

StatsRow RunRecorder :: calcStats(std::vector<double> v) 
{
    StatsRow out;
    v.erase(std::remove_if(v.begin(), v.end(), [](double x){ return std::isnan(x); }), v.end());
    if (v.empty()) return out;

    std::sort(v.begin(), v.end());
    out.min = v.front();
    out.max = v.back();
    out.mean = std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size());

    auto q = [&](double p)->double 
    {
        if (v.size() == 1) return v[0];
        double pos = (static_cast<double>(v.size()) - 1.0) * p;
        auto i = static_cast<size_t>(pos);
        double frac = pos - static_cast<double>(i);
        if (i + 1 >= v.size()) return v.back();
        return (1.0 - frac) * v[i] + frac * v[i + 1];
    };

    out.p10 = q(0.10);
    out.q1  = q(0.25);
    out.med = q(0.50);
    out.q3  = q(0.75);
    out.p90 = q(0.90);
    return out;
}

void RunRecorder::writeParamsHeader(std::ofstream& f, char sep) const
{
    f << "param" << sep << "value\n";

    // ---- Basics / run setup
    f << "filename"    << sep << params_.filename    << "\n";
    f << "base_seed"   << sep << params_.seed        << "\n";
    f << "nRuns"       << sep << params_.nRuns       << "\n";
    f << "iterations"  << sep << params_.iterations  << "\n";
    f << "antsN"       << sep << params_.antsN       << "\n";

    // ---- General controls
    f << "stagnation"  << sep << params_.stagnation  << "\n";
    f << "favorites"   << sep << params_.favorites   << "\n";
    f << "gbPeriod"    << sep << params_.gbPeriod    << "\n";

    // ---- ACO core (evaporation)
    f << "rhoNodes"    << sep << params_.rhoNodes    << "\n";
    f << "rhoCars"     << sep << params_.rhoCars     << "\n";
    f << "rhoPass"     << sep << params_.rhoPass     << "\n";

    // ---- ACO core (pheromone importance)
    f << "alphaNodes"  << sep << params_.alphaNodes  << "\n";
    f << "alphaCars"   << sep << params_.alphaCars   << "\n";
    f << "alphaPass"   << sep << params_.alphaPass   << "\n";

    // ---- ACO core (heuristic importance)
    f << "betaNodes"   << sep << params_.betaNodes   << "\n";
    f << "betaCars"    << sep << params_.betaCars    << "\n";
    f << "betaPass"    << sep << params_.betaPass    << "\n";

    // ---- Limits / MMAS-ish controls
    f << "maxMinCars"  << sep << params_.maxMinCars  << "\n";
    f << "maxMinNodes" << sep << params_.maxMinNodes << "\n";
    f << "pBestCars"   << sep << params_.pBestCars   << "\n";

    // ---- Pheromone model selection
    f << "useTBAS"     << sep << (params_.useTBAS ? 1 : 0) << "\n";

    // ---- Construction / policy toggles (NEW A/B/C)
    f << "sampleCars"       << sep << (params_.sampleCars ? 1 : 0) << "\n";
    f << "enforceNoRerent"  << sep << (params_.enforceNoRerent ? 1 : 0) << "\n";
    f << "returnFactorCars" << sep << params_.returnFactorCars << "\n";

    // ---- Cars: heuristic / eta
    f << "epsilonEtaCars" << sep << params_.epsilonEtaCars << "\n";

    // ---- Cars: stagnation / restart / explore (MMAS extras)
    f << "smGammaCars"        << sep << params_.smoothingGammaCars << "\n";
    f << "restartTargetCars"  << sep << aco_cli::toString(params_.restartTargetCars) << "\n";

    f << "exploreCars"        << sep << (params_.exploreEnabledCars ? 1 : 0) << "\n";
    f << "exploreItersCars"   << sep << params_.exploreItersCars << "\n";
    f << "exploreRhoMultCars" << sep << params_.exploreRhoMultiplierCars << "\n";

    // ---- Construction / policy (q0)
    f << "policyQ0"      << sep << params_.policyQ0 << "\n";
    f << "adaptiveQ0"    << sep << (params_.adaptiveQ0 ? 1 : 0) << "\n";
    f << "q0Start"       << sep << params_.q0Start << "\n";
    f << "q0End"          << sep << params_.q0End << "\n";

    // ---- Elite and local search
    f << "eliteKAnts"     << sep << params_.eliteKAnts << "\n";
    f << "lsTopW"         << sep << params_.lsTopW << "\n";

    // ---- LS Quality parameters
    {
        // Opća kvaliteta (1=Low, 2=Medium, 3=High, 4=Extreme)
        const char* qualityNames[] = {"", "Low", "Medium", "High", "Extreme"};
        const int mode = std::clamp(params_.lsQualityMode, 1, 4);
        f << "lsQuality" << sep << qualityNames[mode];
        
        // Override parametri ako su postavljeni (nije -1)
        bool hasOverrides = false;
        if (params_.lsTwoOpt1Passes > 0) hasOverrides = true;
        if (params_.lsRelocPasses > 0) hasOverrides = true;
        if (params_.lsTwoOpt2Passes > 0) hasOverrides = true;
        if (params_.lsDpRetries > 0) hasOverrides = true;
        
        if (hasOverrides) {
            f << " (";
            bool first = true;
            if (params_.lsTwoOpt1Passes > 0) {
                if (!first) f << ",";
                f << "twoOpt1=" << params_.lsTwoOpt1Passes;
                first = false;
            }
            if (params_.lsRelocPasses > 0) {
                if (!first) f << ",";
                f << "reloc=" << params_.lsRelocPasses;
                first = false;
            }
            if (params_.lsTwoOpt2Passes > 0) {
                if (!first) f << ",";
                f << "twoOpt2=" << params_.lsTwoOpt2Passes;
                first = false;
            }
            if (params_.lsDpRetries > 0) {
                if (!first) f << ",";
                f << "dpRetries=" << params_.lsDpRetries;
            }
            f << ")";
        }
        f << "\n";
    }

    f << "useReinsertion"     << sep << (params_.useReinsertion ? 1 : 0) << "\n";
    f << "reinsertionMaxPasses" << sep << params_.reinsertionMaxPasses << "\n";

    // ---- DP cache and tuning
    f << "dpCacheCapacity"        << sep << params_.dpCacheCapacity << "\n";
    f << "dpMaxSegmentLenOffset"  << sep << params_.dpMaxSegmentLenOffset << "\n";

    // ---- Intensifier parameters
    f << "intensifierEnabled" << sep << (params_.intensifierEnabled ? 1 : 0) << "\n";
    f << "intensifierMaxStagnationIterations" << sep << params_.intensifierMaxStagnationIterations << "\n";
    f << "intensifierResetPheromonesCount"    << sep << params_.intensifierResetPheromonesCount << "\n";
    f << "intensifierPheromoneReductionGamma" << sep << params_.intensifierPheromoneReductionGamma << "\n";
    f << "intensifierArchiveAndResetBest"     << sep << (params_.intensifierArchiveAndResetBest ? 1 : 0) << "\n";
    f << "ablationDisableAntDepositsDuringIntensification" << sep
      << (params_.ablationDisableAntDepositsDuringIntensification ? 1 : 0) << "\n";
    f << "ablationResetPheromonesOnIntensifierActivation" << sep
      << (params_.ablationResetPheromonesOnIntensifierActivation ? 1 : 0) << "\n";
    f << "ablationIntensifierOnlyNoAnts" << sep << (params_.ablationIntensifierOnlyNoAnts ? 1 : 0) << "\n";
    f << "ablationIntensifierOnlyRandomInitialTour" << sep
      << (params_.ablationIntensifierOnlyRandomInitialTour ? 1 : 0) << "\n";

    // ---- Run seeds
    f << "runSeeds" << sep;
    for (int r = 0; r < runs_; ++r) {
        if (r) f << ",";
        f << "run" << r << "=" << seeds_[r];
    }
    f << "\n";

    // ---- Run counters
    f << "stagnationEvents" << sep;
    for (int r = 0; r < runs_; ++r) {
        if (r) f << ",";
        f << "run" << r << "=" << stagnationEvents_[r];
    }
    f << "\n";

    f << "lkCalls" << sep;
    for (int r = 0; r < runs_; ++r) {
        if (r) f << ",";
        f << "run" << r << "=" << lkCalls_[r];
    }
    f << "\n";

    f << "lkSuccess" << sep;
    for (int r = 0; r < runs_; ++r) {
        if (r) f << ",";
        f << "run" << r << "=" << lkSuccess_[r];
    }
    f << "\n";

    f << "tOptCalls" << sep;
    for (int r = 0; r < runs_; ++r) {
        if (r) f << ",";
        f << "run" << r << "=" << toptCalls_[r];
    }
    f << "\n";

    f << "tOptSuccess" << sep;
    for (int r = 0; r < runs_; ++r) {
        if (r) f << ",";
        f << "run" << r << "=" << toptSuccess_[r];
    }
    f << "\n";

    f << "kickCalls" << sep;
    for (int r = 0; r < runs_; ++r) {
        if (r) f << ",";
        f << "run" << r << "=" << kickCalls_[r];
    }
    f << "\n";

    f << "kickSuccess" << sep;
    for (int r = 0; r < runs_; ++r) {
        if (r) f << ",";
        f << "run" << r << "=" << kickSuccess_[r];
    }
    f << "\n";

    // ---- Cache statistics (final hit rates per run)
    f << "cacheEvalHitRate" << sep;
    for (int r = 0; r < runs_; ++r) {
        if (r) f << ",";
        if (std::isfinite(cacheEvalHitRate_[r])) {
            f << "run" << r << "=" << std::fixed << std::setprecision(2) << cacheEvalHitRate_[r] << "%";
        } else {
            f << "run" << r << "=N/A";
        }
    }
    f << "\n";

    f << "cacheReasHitRate" << sep;
    for (int r = 0; r < runs_; ++r) {
        if (r) f << ",";
        if (std::isfinite(cacheReasHitRate_[r])) {
            f << "run" << r << "=" << std::fixed << std::setprecision(2) << cacheReasHitRate_[r] << "%";
        } else {
            f << "run" << r << "=N/A";
        }
    }
    f << "\n";

    // ---- LK and 3-opt success rates (with cache eval hit rate if available)
    f << "lkSuccessRate" << sep;
    for (int r = 0; r < runs_; ++r) {
        if (r) f << ",";
        if (lkCalls_[r] > 0) {
            double rate = 100.0 * (double)lkSuccess_[r] / (double)lkCalls_[r];
            f << "run" << r << "=" << std::fixed << std::setprecision(2) << rate << "%";
            if (std::isfinite(cacheEvalHitRate_[r])) {
                f << " (cache:" << std::fixed << std::setprecision(1) << cacheEvalHitRate_[r] << "%)";
            }
        } else {
            f << "run" << r << "=N/A";
        }
    }
    f << "\n";

    f << "tOptSuccessRate" << sep;
    for (int r = 0; r < runs_; ++r) {
        if (r) f << ",";
        if (toptCalls_[r] > 0) {
            double rate = 100.0 * (double)toptSuccess_[r] / (double)toptCalls_[r];
            f << "run" << r << "=" << std::fixed << std::setprecision(2) << rate << "%";
            if (std::isfinite(cacheEvalHitRate_[r])) {
                f << " (cache:" << std::fixed << std::setprecision(1) << cacheEvalHitRate_[r] << "%)";
            }
        } else {
            f << "run" << r << "=N/A";
        }
    }
    f << "\n";

    // separator empty line
    f << "\n";
}


void RunRecorder :: writeCsv(const std::string& filename, char sep, bool useGlobal) const 
{
    std::ofstream f(filename);
    f.imbue(std::locale::classic()); // stabilno; ako hoćeš decimal comma, javi pa dodamo custom facet

    writeParamsHeader(f, sep);

    // Stupci:
    // iter; time_ms_mean; time_ms_min; time_ms_max; run0;...; stats...
    f << "iter"
        << sep << "time_ms_mean"
        << sep << "time_ms_min"
        << sep << "time_ms_max";

    for (int r = 0; r < runs_; ++r) f << sep << "run" << r;

    f << sep << "min"
        << sep << "p10"
        << sep << "q1"
        << sep << "median"
        << sep << "q3"
        << sep << "p90"
        << sep << "max"
        << sep << "mean"
        << "\n";

    std::vector<double> costBuf; costBuf.reserve(static_cast<size_t>(runs_));
    std::vector<double> timeBuf; timeBuf.reserve(static_cast<size_t>(runs_));

    for (int it = 0; it < iters_; ++it) 
    {
        f << it;

        costBuf.clear();
        timeBuf.clear();

        for (int r = 0; r < runs_; ++r) 
        {
            const size_t idx = index(r, it);
            // Za iteration_best koristi iterBestPost_ (nakon svih akcija u iteraciji)
            // Za global_best koristi globalBest_
            const double c = useGlobal ? globalBest_[idx] : iterBestPost_[idx];
            const double t = elapsedMs_[idx];
            costBuf.push_back(c);
            timeBuf.push_back(t);
        }

        const StatsRow tSt = calcStats(timeBuf);
        f << sep << tSt.mean << sep << tSt.min << sep << tSt.max;

        for (int r = 0; r < runs_; ++r) 
        {
            const size_t idx = index(r, it);
            // Za iteration_best koristi iterBestPost_ (nakon svih akcija u iteraciji)
            // Za global_best koristi globalBest_
            const double c = useGlobal ? globalBest_[idx] : iterBestPost_[idx];
            f << sep << c;
        }

        const StatsRow cSt = calcStats(costBuf);
        f << sep << cSt.min
            << sep << cSt.p10
            << sep << cSt.q1
            << sep << cSt.med
            << sep << cSt.q3
            << sep << cSt.p90
            << sep << cSt.max
            << sep << cSt.mean
            << "\n";
    }
};

void RunRecorder::writeWeibullCsv(const std::string& filename, char sep) const
{
    std::ofstream f(filename);
    f.imbue(std::locale::classic());

    writeParamsHeader(f, sep);

    // Stupci: iter; run0;...; stats za k
    f << "iter";
    for (int r = 0; r < runs_; ++r) f << sep << "weibull_k_run" << r;

    f << sep << "min"
      << sep << "p10"
      << sep << "q1"
      << sep << "median"
      << sep << "q3"
      << sep << "p90"
      << sep << "max"
      << sep << "mean"
      << "\n";

    std::vector<double> kBuf; 
    kBuf.reserve(static_cast<size_t>(runs_));

    for (int it = 0; it < iters_; ++it)
    {
        f << it;

        kBuf.clear();
        for (int r = 0; r < runs_; ++r)
        {
            const size_t idx = index(r, it);
            const double k = weibullK_[idx];
            f << sep << k;
            kBuf.push_back(k);
        }

        const StatsRow kSt = calcStats(kBuf);
        f << sep << kSt.min
          << sep << kSt.p10
          << sep << kSt.q1
          << sep << kSt.med
          << sep << kSt.q3
          << sep << kSt.p90
          << sep << kSt.max
          << sep << kSt.mean
          << "\n";
    }
}

void RunRecorder::recordGlobalBestImprovement(int runIdx, int iteration)
{
    improvementIters_.at(runIdx).push_back(iteration);
}

// MLE for Weibull(k, lambda) from positive waiting times.
// Profile likelihood: h(k) = [sum(Y^k ln Y)] / [sum(Y^k)] - 1/k - (1/M) sum(ln Y) = 0
// h(k) is monotonically increasing ⇒ bisection is safe.
static std::pair<double, double> weibullMLE(const std::vector<double>& Y)
{
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    const size_t M = Y.size();
    if (M < 3) return {NaN, NaN};

    double sumLogY = 0.0;
    for (double y : Y) {
        if (y <= 0.0) return {NaN, NaN};
        sumLogY += std::log(y);
    }
    const double C = sumLogY / static_cast<double>(M);

    auto hFunc = [&](double k) -> double {
        double A = 0.0, B = 0.0;
        for (double y : Y) {
            const double yk = std::pow(y, k);
            A += yk * std::log(y);
            B += yk;
        }
        if (B <= 0.0) return NaN;
        return A / B - 1.0 / k - C;
    };

    double lo = 0.01, hi = 200.0;
    double hLo = hFunc(lo), hHi = hFunc(hi);
    if (!std::isfinite(hLo) || !std::isfinite(hHi)) return {NaN, NaN};
    if (hLo * hHi > 0.0) return {NaN, NaN};

    for (int step = 0; step < 100; ++step) {
        const double mid = 0.5 * (lo + hi);
        const double hMid = hFunc(mid);
        if (!std::isfinite(hMid)) break;
        if (std::abs(hMid) < 1e-12) { lo = hi = mid; break; }
        if (hMid * hLo < 0.0) { hi = mid; hHi = hMid; }
        else                   { lo = mid; hLo = hMid; }
    }
    const double k = 0.5 * (lo + hi);

    double sumYk = 0.0;
    for (double y : Y) sumYk += std::pow(y, k);
    const double lambda = std::pow(sumYk / static_cast<double>(M), 1.0 / k);

    if (!std::isfinite(k) || !std::isfinite(lambda)) return {NaN, NaN};
    return {k, lambda};
}

void RunRecorder::writeWeibullWaitingTimeCsv(const std::string& filename, char sep) const
{
    std::ofstream f(filename);
    f.imbue(std::locale::classic());

    writeParamsHeader(f, sep);

    f << "run" << sep << "k_mle" << sep << "lambda_mle" << sep << "M"
      << sep << "num_improvements" << sep << "mean_gap" << sep << "median_gap"
      << sep << "improvements" << sep << "waiting_times" << "\n";

    std::vector<double> allWT;
    std::vector<double> kPerRun;
    kPerRun.reserve(static_cast<size_t>(runs_));

    for (int r = 0; r < runs_; ++r)
    {
        const auto& imps = improvementIters_[r];

        std::vector<double> wt;
        if (!imps.empty()) {
            if (imps[0] > 0) wt.push_back(static_cast<double>(imps[0]));
            for (size_t i = 1; i < imps.size(); ++i) {
                const double gap = static_cast<double>(imps[i] - imps[i - 1]);
                if (gap > 0.0) wt.push_back(gap);
            }
        }

        const auto [k, lam] = weibullMLE(wt);
        kPerRun.push_back(k);

        double meanGap = StatsRow::nan(), medGap = StatsRow::nan();
        if (!wt.empty()) {
            meanGap = std::accumulate(wt.begin(), wt.end(), 0.0) / static_cast<double>(wt.size());
            auto wtSorted = wt;
            std::sort(wtSorted.begin(), wtSorted.end());
            const size_t n = wtSorted.size();
            medGap = (n % 2 == 1) ? wtSorted[n / 2]
                                   : 0.5 * (wtSorted[n / 2 - 1] + wtSorted[n / 2]);
        }

        f << r << sep << k << sep << lam << sep << wt.size()
          << sep << imps.size() << sep << meanGap << sep << medGap;

        f << sep;
        for (size_t i = 0; i < imps.size(); ++i) {
            if (i) f << ",";
            f << imps[i];
        }

        f << sep;
        for (size_t i = 0; i < wt.size(); ++i) {
            if (i) f << ",";
            f << wt[i];
        }
        f << "\n";

        allWT.insert(allWT.end(), wt.begin(), wt.end());
    }

    // Aggregated Weibull fit across all runs
    const auto [kAll, lamAll] = weibullMLE(allWT);
    double meanAll = StatsRow::nan(), medAll = StatsRow::nan();
    if (!allWT.empty()) {
        meanAll = std::accumulate(allWT.begin(), allWT.end(), 0.0)
                / static_cast<double>(allWT.size());
        auto sorted = allWT;
        std::sort(sorted.begin(), sorted.end());
        const size_t n = sorted.size();
        medAll = (n % 2 == 1) ? sorted[n / 2]
                               : 0.5 * (sorted[n / 2 - 1] + sorted[n / 2]);
    }
    f << "all" << sep << kAll << sep << lamAll << sep << allWT.size()
      << sep << "" << sep << meanAll << sep << medAll
      << sep << "" << sep << "" << "\n";

    // Per-run k summary statistics
    const StatsRow kSt = calcStats(kPerRun);
    f << "\n"
      << "k_stats" << sep << "min" << sep << "p10" << sep << "q1"
      << sep << "median" << sep << "q3" << sep << "p90" << sep << "max"
      << sep << "mean" << "\n"
      << "" << sep << kSt.min << sep << kSt.p10 << sep << kSt.q1
      << sep << kSt.med << sep << kSt.q3 << sep << kSt.p90 << sep << kSt.max
      << sep << kSt.mean << "\n";
}

void RunRecorder::recordDiversity(int runIdx, int it, double costStd,
                                   double costRange, double costMean,
                                   double costCV, bool improved)
{
    const size_t idx = index(runIdx, it);
    costStd_[idx]         = costStd;
    costRange_[idx]       = costRange;
    costMean_[idx]        = costMean;
    costCV_[idx]          = costCV;
    improvedThisIter_[idx] = improved;
}

void RunRecorder::recordAntQuantiles(int runIdx, int it, double minCost, double q1Cost, double medianCost,
                                     double q3Cost, double maxCost, double meanCost, double stdCost, int sampleN)
{
    const size_t idx = index(runIdx, it);
    antCostMin_[idx] = minCost;
    antCostQ1_[idx] = q1Cost;
    antCostMedian_[idx] = medianCost;
    antCostQ3_[idx] = q3Cost;
    antCostMax_[idx] = maxCost;
    antCostMean_[idx] = meanCost;
    antCostStd_[idx] = stdCost;
    antCostN_[idx] = sampleN;
}

void RunRecorder::recordIntensifierDeactivationEvent(int runIdx, int it, int activationCount, double globalBest)
{
    intensifierEvents_.push_back(IntensifierEvent{
        runIdx, it, "intensifier_deactivated", activationCount, -1, globalBest
    });
}

void RunRecorder::recordStagnationResetGlobalImprovementEvent(int runIdx, int it, int noImproveBefore, double globalBest)
{
    intensifierEvents_.push_back(IntensifierEvent{
        runIdx, it, "stagnation_reset_global_improvement", -1, noImproveBefore, globalBest
    });
}

void RunRecorder::writeWeibullHypothesesCsv(const std::string& filename, char sep) const
{
    std::ofstream f(filename);
    f.imbue(std::locale::classic());

    writeParamsHeader(f, sep);

    // Per-iteration, per-run "long format": one row = one (run, iteration) pair
    f << "run" << sep << "iter"
      << sep << "weibull_k"
      << sep << "cost_std"
      << sep << "cost_range"
      << sep << "cost_mean"
      << sep << "cost_cv"
      << sep << "global_best"
      << sep << "improved"
      << sep << "cumul_improvements"
      << "\n";

    for (int r = 0; r < runs_; ++r)
    {
        int cumulImp = 0;
        for (int it = 0; it < iters_; ++it)
        {
            const size_t idx = index(r, it);
            const bool imp = improvedThisIter_[idx];
            if (imp) ++cumulImp;

            f << r << sep << it
              << sep << weibullK_[idx]
              << sep << costStd_[idx]
              << sep << costRange_[idx]
              << sep << costMean_[idx]
              << sep << costCV_[idx]
              << sep << globalBest_[idx]
              << sep << (imp ? 1 : 0)
              << sep << cumulImp
              << "\n";
        }
    }
}

void RunRecorder::writeIntensifierEventsCsv(const std::string& filename, char sep) const
{
    std::ofstream f(filename);
    f.imbue(std::locale::classic());

    writeParamsHeader(f, sep);
    f << "run" << sep << "iter" << sep << "event_type" << sep
      << "activation_count" << sep << "no_improve_before" << sep
      << "global_best" << "\n";

    for (const auto& ev : intensifierEvents_)
    {
        f << ev.run << sep << ev.iter << sep << ev.eventType << sep
          << ev.activationCount << sep << ev.noImproveBefore << sep
          << ev.globalBest << "\n";
    }
}

void RunRecorder::writeAntsQuantilesCsv(const std::string& filename, char sep) const
{
    std::ofstream f(filename);
    f.imbue(std::locale::classic());

    writeParamsHeader(f, sep);
    f << "run" << sep << "iter"
      << sep << "n_ants"
      << sep << "ant_cost_min"
      << sep << "ant_cost_q1"
      << sep << "ant_cost_median"
      << sep << "ant_cost_q3"
      << sep << "ant_cost_max"
      << sep << "ant_cost_mean"
      << sep << "ant_cost_std"
      << "\n";

    for (int r = 0; r < runs_; ++r)
    {
        for (int it = 0; it < iters_; ++it)
        {
            const size_t idx = index(r, it);
            f << r << sep << it
              << sep << antCostN_[idx]
              << sep << antCostMin_[idx]
              << sep << antCostQ1_[idx]
              << sep << antCostMedian_[idx]
              << sep << antCostQ3_[idx]
              << sep << antCostMax_[idx]
              << sep << antCostMean_[idx]
              << sep << antCostStd_[idx]
              << "\n";
        }
    }
}

double RunRecorder::lastFiniteFromBack(const std::vector<double>& arr, int runIdx) const
{
    for (int it = iters_ - 1; it >= 0; --it) 
    {
        const double x = arr[index(runIdx, it)];
        if (std::isfinite(x)) return x;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

double RunRecorder::minFiniteAll(const std::vector<double>& arr, int runIdx) const
{
    double best = std::numeric_limits<double>::infinity();
    bool ok = false;

    for (int it = 0; it < iters_; ++it) 
    {
        const double x = arr[index(runIdx, it)];
        if (std::isfinite(x)) {
            best = std::min(best, x);
            ok = true;
        }
    }
    return ok ? best : std::numeric_limits<double>::quiet_NaN();
}


void RunRecorder::writeRunSummaryCsv(const std::string& filename, char sep) const
{
    std::ofstream f(filename);
    f.imbue(std::locale::classic());

    writeParamsHeader(f, sep);

    f << "run"
      << sep << "seed"
      << sep << "final_global_best"
      << sep << "min_global_best"
      << sep << "final_iter_best_pre"
      << sep << "min_iter_best_pre"
      << sep << "final_iter_best_post"
      << sep << "min_iter_best_post"
      << sep << "total_elapsed_ms"
      << sep << "stagnation_events"
      << sep << "lk_calls"
      << sep << "lk_success"
      << sep << "tOpt_calls"
      << sep << "tOpt_success"
      << sep << "kick_calls"
      << sep << "kick_success"
      << "\n";

    for (int r = 0; r < runs_; ++r)
    {
        const double finalGlobal = lastFiniteFromBack(globalBest_, r);
        const double minGlobal   = minFiniteAll(globalBest_, r);

        const double finalPre    = lastFiniteFromBack(iterBest_, r);
        const double minPre      = minFiniteAll(iterBest_, r);

        const double finalPost   = lastFiniteFromBack(iterBestPost_, r);
        const double minPost     = minFiniteAll(iterBestPost_, r);

        const double totalMs     = lastFiniteFromBack(elapsedMs_, r);

        f << r
          << sep << seeds_.at(r)
          << sep << finalGlobal
          << sep << minGlobal
          << sep << finalPre
          << sep << minPre
          << sep << finalPost
          << sep << minPost
          << sep << totalMs
          << sep << stagnationEvents_.at(r)
          << sep << lkCalls_.at(r)
          << sep << lkSuccess_.at(r)
          << sep << toptCalls_.at(r)
          << sep << toptSuccess_.at(r)
          << sep << kickCalls_.at(r)
          << sep << kickSuccess_.at(r)
          << "\n";
    }
}

void RunRecorder::writeDynamicParamsCsv(const std::string& filename, char sep) const
{
    std::ofstream f(filename);
    f.imbue(std::locale::classic());

    // Header: iter, run0_q0, run0_LSW, run0_smoothing, run1_q0, run1_LSW, run1_smoothing, ...
    f << "iter";
    for (int r = 0; r < runs_; ++r) {
        f << sep << "run" << r << "_q0"
          << sep << "run" << r << "_LSW"
          << sep << "run" << r << "_smoothing";
    }
    f << "\n";

    // Rows: one per iteration
    for (int it = 0; it < iters_; ++it) {
        f << it;
        for (int r = 0; r < runs_; ++r) {
            const size_t idx = index(r, it);
            
            // q0
            const double q0 = q0_[idx];
            if (std::isfinite(q0)) {
                f << sep << q0;
            } else {
                f << sep << "";
            }
            
            // LSW
            const int lsw = LSW_[idx];
            if (lsw >= 0) {
                f << sep << lsw;
            } else {
                f << sep << "";
            }
            
            // smoothing (0 or 1)
            const bool smoothing = smoothing_[idx];
            f << sep << (smoothing ? 1 : 0);
        }
        f << "\n";
    }
}

void RunRecorder::writeSolutionsCsv(const std::string& filename, 
                                     const std::vector<RunResult>& results, 
                                     char sep) const
{
    std::ofstream f(filename);
    f.imbue(std::locale::classic());

    // Extract problem name from filename (without path and extension)
    std::string problemName = params_.filename;
    // Remove path
    size_t lastSlash = problemName.find_last_of("/\\");
    if (lastSlash != std::string::npos)
        problemName = problemName.substr(lastSlash + 1);
    // Remove extension
    size_t lastDot = problemName.find_last_of(".");
    if (lastDot != std::string::npos)
        problemName = problemName.substr(0, lastDot);

    // Write problem name at the beginning
    f << "problem" << sep << problemName << "\n";

    // Write key parameters (similar to writeParamsHeader but simplified)
    f << "useTBAS" << sep << (params_.useTBAS ? 1 : 0) << "\n";

    // Header (only for first row with run/seed/cost/time)
    f << "run" << sep << "seed" << sep << "cost" << sep << "time_ms" << "\n";

    for (const auto& result : results)
    {
        // Row 1: run, seed, cost, time_ms
        f << result.runIndex << sep << result.seed << sep;
        f << std::setprecision(15) << result.best.cost << sep;
        
        // Get elapsed time for this run (last iteration)
        const double totalMs = lastFiniteFromBack(elapsedMs_, result.runIndex);
        f << totalMs << "\n";

        // Row 2: nodes (comma-separated)
        const auto& nodes = result.best.sol.node;
        if (!nodes.empty())
        {
            f << nodes[0];
            for (size_t i = 1; i < nodes.size(); ++i)
                f << "," << nodes[i];
        }
        f << "\n";

        // Row 3: cars (comma-separated)
        const auto& cars = result.best.sol.car;
        if (!cars.empty())
        {
            f << cars[0];
            for (size_t i = 1; i < cars.size(); ++i)
                f << "," << cars[i];
        }
        else
        {
            // TSP case: no cars or all zeros
            if (!nodes.empty())
            {
                f << "0";
                for (size_t i = 1; i < nodes.size(); ++i)
                    f << ",0";
            }
        }
        f << "\n";
    }

    // Add best run at the end
    if (!results.empty())
    {
        auto bestIt = std::min_element(results.begin(), results.end(),
            [](const RunResult& a, const RunResult& b) 
            {
                return a.best.cost < b.best.cost;
            });

        if (bestIt != results.end())
        {
            // Row 1: best run marker, seed, cost, time_ms
            f << "best" << sep << bestIt->seed << sep;
            f << std::setprecision(15) << bestIt->best.cost << sep;
            const double totalMs = lastFiniteFromBack(elapsedMs_, bestIt->runIndex);
            f << totalMs << "\n";

            // Row 2: nodes (comma-separated)
            const auto& nodes = bestIt->best.sol.node;
            if (!nodes.empty())
            {
                f << nodes[0];
                for (size_t i = 1; i < nodes.size(); ++i)
                    f << "," << nodes[i];
            }
            f << "\n";

            // Row 3: cars (comma-separated)
            const auto& cars = bestIt->best.sol.car;
            if (!cars.empty())
            {
                f << cars[0];
                for (size_t i = 1; i < cars.size(); ++i)
                    f << "," << cars[i];
            }
            else
            {
                // TSP case: no cars or all zeros
                if (!nodes.empty())
                {
                    f << "0";
                    for (size_t i = 1; i < nodes.size(); ++i)
                        f << ",0";
                }
            }
            f << "\n";
        }
    }
}

double RunRecorder::getElapsedTimeForRun(int runIdx) const
{
    return lastFiniteFromBack(elapsedMs_, runIdx);
}

double RunRecorder::getMedianGlobalBestForIteration(int it) const
{
    if (it < 0 || it >= iters_) return StatsRow::nan();
    
    std::vector<double> values;
    values.reserve(static_cast<size_t>(runs_));
    
    for (int r = 0; r < runs_; ++r)
    {
        const size_t idx = index(r, it);
        const double val = globalBest_[idx];
        if (std::isfinite(val))
        {
            values.push_back(val);
        }
    }
    
    if (values.empty()) return StatsRow::nan();
    
    const StatsRow stats = calcStats(values);
    return stats.med;
}

double RunRecorder::getMeanGlobalBestForIteration(int it) const
{
    if (it < 0 || it >= iters_) return StatsRow::nan();
    
    std::vector<double> values;
    values.reserve(static_cast<size_t>(runs_));
    
    // Samo uključi runove koji su već izračunali svoju vrijednost (sporiji se ignoriraju)
    for (int r = 0; r < runs_; ++r)
    {
        const size_t idx = index(r, it);
        const double val = globalBest_[idx];
        if (std::isfinite(val))
        {
            values.push_back(val);
        }
    }
    
    if (values.empty()) return StatsRow::nan();
    
    // Izračunaj srednju vrijednost
    double sum = 0.0;
    for (double val : values)
    {
        sum += val;
    }
    return sum / static_cast<double>(values.size());
}

}