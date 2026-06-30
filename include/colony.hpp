#pragma once
#include <memory>
#include <vector>
#include <cstdint>
#include <array>
#include <algorithm>
#include <functional>
#include <rng.hpp>
#include <interfaces.hpp>

// forward declarations
namespace cars_tsplib { struct Instance; }
namespace aco_cli { struct ParamsData; }
namespace aco { class RunRecorder; class Logger; class IntensifierBinomialLogger; class SurrogateCorrelationLogger; }
namespace localSearch { class LkLiteLocalSearch; }

namespace aco 
{

struct RunResult 
{
    int runIndex = -1;
    uint64_t seed = 0;

    EvaluatedSolution best;                 // global best of this run
    std::vector<double> iterBestCost;       // iteration-best cost per iteration
};

class Colony 
{
public:
    Colony(std::shared_ptr<const cars_tsplib::Instance> inst,
           const aco_cli::ParamsData& params,
           int runIndex,
           uint64_t seed,
           std::shared_ptr<const IFixedTourCarAssigner> dp,
           std::shared_ptr<const IAntPolicy> antPolicy,
           std::shared_ptr<IPheromoneModel> pherModel,
           std::shared_ptr<const ICostModel> costModel,
           std::shared_ptr<const ILocalSearch> localSearch,
           std::shared_ptr<const ILocalSearch> lkLite  = nullptr,
           std::shared_ptr<const ILocalSearch> threeOptLite = nullptr,
           std::shared_ptr<const IIntensifier> intensifier = nullptr);  // nullptr + intensifierEnabled => auto-create; -intf 0 keeps it disabled

    RunResult run();
    void setRecorder(std::shared_ptr<aco::RunRecorder> r) { recorder_ = std::move(r); }
    void setLogger(std::shared_ptr<aco::Logger> logger) { logger_ = std::move(logger); }
    void setResearchLogger(std::shared_ptr<aco::IntensifierBinomialLogger> r) { researchLogger_ = std::move(r); }
    void setResearchSurrogateLogger(std::shared_ptr<aco::SurrogateCorrelationLogger> r) { researchSurrogateLogger_ = std::move(r); }

    // for console/log output
    void setVerbose(bool v, std::mutex* mx) { verbose_ = v; printMx_ = mx; }
    void updateVerbose(bool v) { verbose_ = v; } // update verbose flag while running

    // Callback for reporting current iteration (progress tracking)
    using ProgressCallback = std::function<void(int iteration)>;
    void setProgressCallback(ProgressCallback cb) { progressCallback_ = std::move(cb); }

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    const aco_cli::ParamsData& params_;
    int runIndex_ = 0;
    uint64_t seed_ = 0;
    Rng rng_;

    std::shared_ptr<const aco::IFixedTourCarAssigner> dp_;

    std::shared_ptr<const IAntPolicy> antPolicy_;
    std::shared_ptr<IPheromoneModel> pherModel_;
    std::shared_ptr<const ICostModel> costModel_;
    std::shared_ptr<const ILocalSearch> localSearch_;
    std::shared_ptr<const aco::ILocalSearch> lkLite_;
    std::shared_ptr<const aco::ILocalSearch> threeOptLite_;
    std::shared_ptr<const IIntensifier> intensifier_;
    std::shared_ptr<aco::RunRecorder> recorder_;
    std::shared_ptr<aco::Logger> logger_; // logger for file and console output
    std::shared_ptr<aco::IntensifierBinomialLogger> researchLogger_; // research logging (binomial analysis, ablations)
    std::shared_ptr<aco::SurrogateCorrelationLogger> researchSurrogateLogger_; // surrogate vs DP cost correlation logging

    // for console/log output
    bool verbose_ = false;
    std::mutex* printMx_ = nullptr;
    ProgressCallback progressCallback_; // progress-reporting callback

    // Precomputed min travel cost cache (for surrogate evaluation)
    std::vector<std::vector<double>> minTravelCache_;
    double avgReturnCost_ = 0.0;  // average return cost (for surrogate's return term)
    void precomputeMinTravel_();
};

} // namespace aco
