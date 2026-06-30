#pragma once
#include <memory>
#include <vector>
#include <cstdint>
#include <array>
#include <algorithm>
#include <functional>
#include <rng.hpp>
#include <interfaces.hpp>

//forward declaration
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

    EvaluatedSolution best;                 // global best ovog run-a
    std::vector<double> iterBestCost;       // po iteraciji (iteration-best cost)
    // kasnije: iterGlobalBestCost, vremena, itd.
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
           std::shared_ptr<const IIntensifier> intensifier = nullptr);  // nullptr + intensifierEnabled => auto-create; -intf 0 => ostaje isključen

    RunResult run();
    void setRecorder(std::shared_ptr<aco::RunRecorder> r) { recorder_ = std::move(r); }
    void setLogger(std::shared_ptr<aco::Logger> logger) { logger_ = std::move(logger); }
    void setResearchLogger(std::shared_ptr<aco::IntensifierBinomialLogger> r) { researchLogger_ = std::move(r); }
    void setResearchSurrogateLogger(std::shared_ptr<aco::SurrogateCorrelationLogger> r) { researchSurrogateLogger_ = std::move(r); }

    // za ispis
    void setVerbose(bool v, std::mutex* mx) { verbose_ = v; printMx_ = mx; }
    void updateVerbose(bool v) { verbose_ = v; } // za promjenu verbose flag-a tijekom izvršavanja
    
    // Callback za ažuriranje trenutne iteracije (za praćenje napretka)
    using ProgressCallback = std::function<void(int iteration)>;
    void setProgressCallback(ProgressCallback cb) { progressCallback_ = std::move(cb); }

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    const aco_cli::ParamsData& params_;
    int runIndex_ = 0;
    uint64_t seed_ = 0;
    Rng rng_;

    //sad i kolonija koristi DP
    std::shared_ptr<const aco::IFixedTourCarAssigner> dp_;


    std::shared_ptr<const IAntPolicy> antPolicy_;
    std::shared_ptr<IPheromoneModel> pherModel_;
    std::shared_ptr<const ICostModel> costModel_;
    std::shared_ptr<const ILocalSearch> localSearch_;
    std::shared_ptr<const aco::ILocalSearch> lkLite_;
    std::shared_ptr<const aco::ILocalSearch> threeOptLite_;
    std::shared_ptr<const IIntensifier> intensifier_;
    std::shared_ptr<aco::RunRecorder> recorder_;
    std::shared_ptr<aco::Logger> logger_; // logger za zapisivanje u fajl i konzolu
    std::shared_ptr<aco::IntensifierBinomialLogger> researchLogger_; // istraživačko logiranje (binomna analiza, ablacije)
    std::shared_ptr<aco::SurrogateCorrelationLogger> researchSurrogateLogger_; // surogat vs DP cost za korelaciju

    // za ispis
    bool verbose_ = false;
    std::mutex* printMx_ = nullptr;
    ProgressCallback progressCallback_; // callback za ažuriranje napretka
    
    // Precomputed min travel cost cache (za surrogate evaluaciju)
    std::vector<std::vector<double>> minTravelCache_;
    double avgReturnCost_ = 0.0;  // prosječni return cost (za surogat s return članom)
    void precomputeMinTravel_();
};

} // namespace aco
