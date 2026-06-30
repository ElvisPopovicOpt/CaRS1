#pragma once
#include <vector>
#include <string>
#include <limits>

namespace aco_cli { struct ParamsData; }

namespace aco 
{
// Forward declarations
struct RunResult;
struct Solution;

struct StatsRow 
{
    double min = nan(), p10 = nan(), q1 = nan(), med = nan(), q3 = nan(), p90 = nan(), max = nan(), mean = nan();
    static double nan() { return std::numeric_limits<double>::quiet_NaN(); }
};

class RunRecorder
{
public:
    RunRecorder(const aco_cli::ParamsData& params)
        : params_(params),
          runs_(params.nRuns),
          iters_(params.iterations),
          seeds_(runs_, 0),
          iterBest_(static_cast<size_t>(runs_) * iters_, StatsRow::nan()),
          iterBestPost_(static_cast<size_t>(runs_) * iters_, StatsRow::nan()),
          globalBest_(static_cast<size_t>(runs_) * iters_, StatsRow::nan()),
          elapsedMs_(static_cast<size_t>(runs_) * iters_, StatsRow::nan()),
          weibullK_(static_cast<size_t>(runs_) * iters_, std::numeric_limits<double>::quiet_NaN()),
          q0_(static_cast<size_t>(runs_) * iters_, StatsRow::nan()),
          LSW_(static_cast<size_t>(runs_) * iters_, -1),
          smoothing_(static_cast<size_t>(runs_) * iters_, false),
          // Per-run counters
          stagnationEvents_(runs_, 0),
          lkCalls_(runs_, 0),
          lkSuccess_(runs_, 0),
          toptCalls_(runs_, 0),
          toptSuccess_(runs_, 0),
          kickCalls_(runs_, 0),
          kickSuccess_(runs_, 0),
          cacheEvalHitRate_(runs_, StatsRow::nan()),
          cacheReasHitRate_(runs_, StatsRow::nan()),
          improvementIters_(runs_),
          costStd_(static_cast<size_t>(runs_) * iters_, std::numeric_limits<double>::quiet_NaN()),
          costRange_(static_cast<size_t>(runs_) * iters_, std::numeric_limits<double>::quiet_NaN()),
          costMean_(static_cast<size_t>(runs_) * iters_, std::numeric_limits<double>::quiet_NaN()),
          costCV_(static_cast<size_t>(runs_) * iters_, std::numeric_limits<double>::quiet_NaN()),
          improvedThisIter_(static_cast<size_t>(runs_) * iters_, false),
          antCostMin_(static_cast<size_t>(runs_) * iters_, std::numeric_limits<double>::quiet_NaN()),
          antCostQ1_(static_cast<size_t>(runs_) * iters_, std::numeric_limits<double>::quiet_NaN()),
          antCostMedian_(static_cast<size_t>(runs_) * iters_, std::numeric_limits<double>::quiet_NaN()),
          antCostQ3_(static_cast<size_t>(runs_) * iters_, std::numeric_limits<double>::quiet_NaN()),
          antCostMax_(static_cast<size_t>(runs_) * iters_, std::numeric_limits<double>::quiet_NaN()),
          antCostMean_(static_cast<size_t>(runs_) * iters_, std::numeric_limits<double>::quiet_NaN()),
          antCostStd_(static_cast<size_t>(runs_) * iters_, std::numeric_limits<double>::quiet_NaN()),
          antCostN_(static_cast<size_t>(runs_) * iters_, 0)
    {}

    void setSeed(int runIdx, uint64_t seed) { seeds_.at(runIdx) = seed; }

    // Per-run counters
    void incStagnationEvent(int runIdx) { ++stagnationEvents_.at(runIdx); }

    void incLkCall(int runIdx)          { ++lkCalls_.at(runIdx); }
    void incLkSuccess(int runIdx)       { ++lkSuccess_.at(runIdx); }

    void incThreeOptCall(int runIdx)          { ++toptCalls_.at(runIdx); }
    void incThreeOptSuccess(int runIdx)       { ++toptSuccess_.at(runIdx); }

    void incKickCall(int runIdx)        { ++kickCalls_.at(runIdx); }
    void incKickSuccess(int runIdx)     { ++kickSuccess_.at(runIdx); }

    // Cache statistics (recorded at end of run)
    void setCacheStats(int runIdx, double evalHitRate, double reasHitRate);

    void recordIteration(int runIdx, int it, 
                          double iterBestCost, 
                          double globalBestCost, 
                          double elapsedMs);
    // Records both pre- and post-local-search iteration-best costs
    void recordIteration2(int runIdx, int it,
                          double iterBestPreCost,
                          double iterBestPostCost,
                          double globalBestCost,
                          double elapsedMs,
                          double weibullK = std::numeric_limits<double>::quiet_NaN());
    // Record q0, LSW and smoothing per iteration
    void recordIterationParams(int runIdx, int it, double q0, int LSW, bool smoothingActive = false);
    void writeIterationBestCsv(const std::string& filename, char sep=';') const;
    void writeGlobalBestCsv(const std::string& filename, char sep=';') const;
    void writeRunSummaryCsv(const std::string& filename, char sep=';') const;
    void writeDynamicParamsCsv(const std::string& filename, char sep=';') const;
    // Weibull shape k per iteration/run, cross-sectional from ant costs
    void writeWeibullCsv(const std::string& filename, char sep=';') const;

    // Waiting-time Weibull: records the iteration of a global-best improvement
    void recordGlobalBestImprovement(int runIdx, int iteration);
    // Waiting-time Weibull: MLE fit over gaps between improvements, per run + aggregated
    void writeWeibullWaitingTimeCsv(const std::string& filename, char sep=';') const;

    // Diversity stats per iteration (for the H5/H6 hypotheses)
    void recordDiversity(int runIdx, int it, double costStd, double costRange,
                         double costMean, double costCV, bool improved);
    void writeWeibullHypothesesCsv(const std::string& filename, char sep=';') const;

    // Extra outputs for external Weibull analysis
    void recordIntensifierDeactivationEvent(int runIdx, int it, int activationCount, double globalBest);
    void recordStagnationResetGlobalImprovementEvent(int runIdx, int it, int noImproveBefore, double globalBest);
    void writeIntensifierEventsCsv(const std::string& filename, char sep=';') const;

    // Ant population per (run,iter): cost quantiles (min/Q1/median/Q3/max) + mean/std
    void recordAntQuantiles(int runIdx, int it, double minCost, double q1Cost, double medianCost,
                            double q3Cost, double maxCost, double meanCost, double stdCost, int sampleN);
    void writeAntsQuantilesCsv(const std::string& filename, char sep=';') const;
    
    // Write solutions.csv with best solution per run
    // Format: run, nodes (comma-separated), cars (comma-separated), cost, time_ms
    // At the end, adds the best run overall (marked as "best")
    void writeSolutionsCsv(const std::string& filename, 
                          const std::vector<RunResult>& results, 
                          char sep=';') const;

    // Get elapsed time for a specific run (in milliseconds)
    double getElapsedTimeForRun(int runIdx) const;
    
    // Get median global best cost for a specific iteration across all runs
    // Returns NaN if no valid data available
    double getMedianGlobalBestForIteration(int it) const;
    
    // Get mean (average) global best cost for a specific iteration across all runs
    // Only includes runs that have already computed their value (faster runs)
    // Returns NaN if no valid data available
    double getMeanGlobalBestForIteration(int it) const;

private:
    aco_cli::ParamsData params_;
    int runs_, iters_;
    std::vector<uint64_t> seeds_;
    std::vector<double> iterBest_;
    std::vector<double> iterBestPost_; 
    std::vector<double> globalBest_;
    std::vector<double> elapsedMs_;
    std::vector<double> weibullK_;  // Weibull shape k per iteration (from ant costs)
    std::vector<double> q0_;  // q0 per iteration
    std::vector<int> LSW_;     // LSW per iteration
    std::vector<bool> smoothing_;  // smoothing active per iteration (true if smoothing was applied)

    // Per-run counters
    std::vector<int> stagnationEvents_;
    std::vector<int> lkCalls_;
    std::vector<int> lkSuccess_;
    std::vector<int> toptCalls_;
    std::vector<int> toptSuccess_;
    std::vector<int> kickCalls_;
    std::vector<int> kickSuccess_;

    // Cache statistics (per run)
    std::vector<double> cacheEvalHitRate_;  // eval hit rate (%)
    std::vector<double> cacheReasHitRate_;  // reassign hit rate (%)

    // Waiting-time Weibull: iterations where the global best improved (per run)
    std::vector<std::vector<int>> improvementIters_;

    // Diversity stats per iteration (for H5/H6)
    std::vector<double> costStd_;
    std::vector<double> costRange_;
    std::vector<double> costMean_;
    std::vector<double> costCV_;
    std::vector<bool>   improvedThisIter_;

    // Ant population cost quantiles per iteration/run
    std::vector<double> antCostMin_;
    std::vector<double> antCostQ1_;
    std::vector<double> antCostMedian_;
    std::vector<double> antCostQ3_;
    std::vector<double> antCostMax_;
    std::vector<double> antCostMean_;
    std::vector<double> antCostStd_;
    std::vector<int>    antCostN_;

    struct IntensifierEvent
    {
        int run = -1;
        int iter = -1;
        std::string eventType;
        int activationCount = -1;
        int noImproveBefore = -1;
        double globalBest = std::numeric_limits<double>::quiet_NaN();
    };
    std::vector<IntensifierEvent> intensifierEvents_;

    size_t index(int runIdx, int it) const;
    static StatsRow calcStats(std::vector<double> v);
    void writeParamsHeader(std::ofstream& f, char sep) const;
    void writeCsv(const std::string& filename, char sep, bool useGlobal) const;
    static double lastFiniteInRun(const std::vector<double>& v, int runIdx, int iters_, int runs_ignored_dummy = 0);
    static double minFiniteInRun (const std::vector<double>& v, int runIdx, int iters_);
    private:
    double lastFiniteFromBack(const std::vector<double>& arr, int runIdx) const;
    double minFiniteAll(const std::vector<double>& arr, int runIdx) const;
};

} // namespace aco
