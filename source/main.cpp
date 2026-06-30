#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <memory>
#include <params_cli.hpp>
#include <parser.hpp>
#include <runner.hpp>
#include <colony.hpp>
#include <logger.hpp>
#include <solutionIO.hpp>
#include <runRecorder.hpp>
#include <candidate_list.hpp>
#include <cost_model.hpp>
#include <antPolicy.hpp>
#include <pheromone_model.hpp>
#include <FixedTourCarAssignerDP.hpp>
#include <fixedTourCarAssigner_adapter.hpp>
#include "LkLiteLocalSearch.hpp"
#include "ThreeOptLiteLocalSearch.hpp"
#include "LinKernighanLocalSearch.hpp"
#include "ThreeOptLocalSearch.hpp"
#include "TwoOptLocalSearch.hpp"
#include "ReallocationLocalSearch.hpp"
#include "ReinsertionLocalSearch.hpp"
#include "CachedFixedTourCarAssigner.hpp"
#include <intensifier.hpp>
#include <research/intensifier_binomial_logger.hpp>
#include <research/surrogate_correlation_logger.hpp>

// Helper structure for LS configuration
struct LSConfig {
    int twoOpt1Passes;
    int relocPasses;
    int twoOpt2Passes;
    int dpRetries;
};

static LSConfig getLSConfig(const aco_cli::ParamsData& params) {
    // Mode values: [TwoOpt #1, Relocation, TwoOpt #2, DP Retries]
    static const LSConfig MODES[] = {
        {6, 5, 3, 1},    // Low (1)
        {15, 10, 5, 2},  // Medium (2)
        {30, 20, 10, 3}, // High (3) - default
        {60, 50, 30, 5}  // Extreme (4)
    };

    const int mode = std::clamp(params.lsQualityMode, 1, 4) - 1;
    LSConfig config = MODES[mode];

    // Apply overrides if set
    if (params.lsTwoOpt1Passes > 0) config.twoOpt1Passes = params.lsTwoOpt1Passes;
    if (params.lsRelocPasses > 0) config.relocPasses = params.lsRelocPasses;
    if (params.lsTwoOpt2Passes > 0) config.twoOpt2Passes = params.lsTwoOpt2Passes;
    if (params.lsDpRetries > 0) config.dpRetries = params.lsDpRetries;
    
    return config;
}

void runAco(const aco_cli::ParamsData& params, std::shared_ptr<cars_tsplib::Instance const> inst, std::shared_ptr<aco::Logger> logger)
{
    aco::Runner runner;
    const uint64_t baseSeed = params.seed;

    // Create a shared RunRecorder for all runs
    auto recorder = std::make_shared<aco::RunRecorder>(params);

    // Research logging (intensifier binomial analysis, ablations) - optional
    auto researchLogger = params.researchLog ? std::make_shared<aco::IntensifierBinomialLogger>() : nullptr;
    // Research logging: surrogate vs DP cost correlation - optional, slows down the run
    auto researchSurrogateLogger = params.researchSurrogateLog ? std::make_shared<aco::SurrogateCorrelationLogger>() : nullptr;

    auto results = runner.runAll(params, baseSeed,
    [&](int runIdx, uint64_t seed) -> std::unique_ptr<aco::Colony>
    {
        const int N = inst->n();

        // Candidate list size - use the user's choice, no automatic adjustment
        int candidateListSize = params.favorites;
        if (candidateListSize <= 0) {  // auto mode (favorites <= 0)
            candidateListSize = 40;  // default for all instances
        }
        
        auto cl = std::make_shared<aco::CandidateListCache>(inst, candidateListSize);
        auto cost = std::make_shared<aco::CostModelCars>(inst);

        auto antPolicy = std::make_shared<aco::AntPolicyCandidateListRoulette>(
            cl,
            params.epsilonEtaCars,
            params.returnFactorCars,
            params.sampleCars,
            params.enforceNoRerent,
            params.policyQ0);

        aco::ExplorationOptions explOpt;
        explOpt.enabled = params.exploreEnabledCars;
        explOpt.iters = params.exploreItersCars;
        explOpt.rhoMultiplier = params.exploreRhoMultiplierCars;

        // Choose pheromone model: MMAS (default) or TBAS
        std::shared_ptr<aco::IPheromoneModel> pher;
        if (params.useTBAS) {
            pher = std::make_shared<aco::PheromoneModelTBASMove>(inst, params, explOpt);
        } else {
            pher = std::make_shared<aco::PheromoneModelMMASMove>(inst, params, explOpt);
        }

        // Get LS config from params (user's choice, no automatic adjustment)
        auto lsConfig = getLSConfig(params);

        localSearch::TwoOptOptions opt;
        opt.maxPasses = lsConfig.twoOpt1Passes;
        opt.firstImprovement = true;
        opt.useDontLookBits = true;
        opt.dontLookWakeRadius = 1;
        opt.reassignCarsAtStart = false;
        opt.dpOpt.validateTour = false;
        opt.cand = cl;

        const int C = inst->cars();
        const int minL = (C > 0) ? ((N + C - 1) / C) : N;
        opt.dpOpt.maxSegmentLen = (params.dpMaxSegmentLenOffset < 0) ? -1 : (minL + params.dpMaxSegmentLenOffset);

        opt.dpOpt.fallbackToUnlimitedOnInf = true;
        opt.dpOpt.fallbackMaxRetries = lsConfig.dpRetries;
        opt.dpOpt.fallbackGrowthFactor = 2;

        auto dpImpl = std::make_shared<localSearch::FixedTourCarAssignerDP>(inst, opt.dpOpt);
        auto dpAco = std::make_shared<aco::FixedTourCarAssignerDPAdapter>(dpImpl);

        aco::DPCacheOptions cacheOpt;
        cacheOpt.capacity = params.dpCacheCapacity;
        cacheOpt.storeCars = true;
        cacheOpt.autoDetectSymmetryWhenUnknown = true;
        cacheOpt.symmetryDetectSamples = 2000;
        cacheOpt.threadSafe = true;
        auto dpCached = std::make_shared<aco::CachedFixedTourCarAssigner>(inst, dpAco, cacheOpt);

        auto twoOpt1 = std::make_shared<localSearch::TwoOptLocalSearch>(inst, opt, dpImpl);

        localSearch::RelocationOptions ropt;
        ropt.cand = cl;
        ropt.reassignCarsAtStart = false;
        ropt.maxPasses = lsConfig.relocPasses;
        ropt.firstImprovement = true;
        ropt.maxMoveEvaluations = 3000;
        ropt.minSurrogateGain = 0.0;

        auto reloc = std::make_shared<localSearch::RelocationLocalSearch>(inst, ropt, dpImpl);

        auto opt2 = opt;
        opt2.maxPasses = lsConfig.twoOpt2Passes;
        opt2.reassignCarsAtStart = false;
        auto twoOpt2 = std::make_shared<localSearch::TwoOptLocalSearch>(inst, opt2, dpImpl);

        std::vector<std::shared_ptr<const aco::ILocalSearch>> lsSeq{ twoOpt1, reloc };
        if (params.useReinsertion) {
            localSearch::ReinsertionOptions riopt;
            riopt.maxPasses = params.reinsertionMaxPasses;
            riopt.firstImprovement = true;
            auto reinsert = std::make_shared<localSearch::ReinsertionLocalSearch>(inst, dpImpl, cl, riopt);
            lsSeq.push_back(reinsert);
        }
        lsSeq.push_back(twoOpt2);
        auto ls = std::make_shared<localSearch::ChainedLocalSearch>(std::move(lsSeq));

        std::shared_ptr<const aco::ILocalSearch> lkPolish = nullptr;
        std::shared_ptr<const aco::ILocalSearch> threeOptPolish = nullptr;

        if (params.useFullPolish)
        {
            localSearch::LinKernighanOptions lkopt;
            lkopt.cand = cl;
            lkopt.attempts = 3;
            lkopt.maxDepth = 20;
            lkopt.maxSequences = 1000;
            lkopt.eps = 1e-12;
            lkopt.polish = ls;
            lkPolish = std::make_shared<localSearch::LinKernighanLocalSearch>(inst, lkopt, dpImpl);

            localSearch::ThreeOptOptions topt;
            topt.cand = cl;
            topt.attempts = 3;
            topt.maxPasses = 5;
            topt.eps = 1e-12;
            topt.polish = ls;
            threeOptPolish = std::make_shared<localSearch::ThreeOptLocalSearch>(inst, topt, dpImpl);
        }
        else
        {
            localSearch::LkLiteOptions lkopt;
            lkopt.cand = cl;
            lkopt.attempts = 10;
            lkopt.chainLen = 12;
            lkopt.triesPerStep = 40;
            lkopt.topKPerI = 8;
            lkopt.requireNegativeSurrogate = true;
            lkopt.eps = 1e-12;
            lkopt.polish = ls;
            lkPolish = std::make_shared<localSearch::LkLiteLocalSearch>(inst, lkopt, dpImpl);

            localSearch::ThreeOptLiteOptions topt;
            topt.cand = cl;
            topt.attempts = 10;
            topt.chainLen = 8;
            topt.triesPerStep = 30;
            topt.topKPerI = 8;
            topt.requireNegativeSurrogate = true;
            topt.eps = 1e-12;
            topt.polish = ls;
            threeOptPolish = std::make_shared<localSearch::ThreeOptLiteLocalSearch>(inst, topt, dpImpl);
        }

        // Intensifier: applies an extra LS pass to the global-best solution for deeper
        // exploitation. Activates when the best changes, or for the first N iterations
        // after stagnation (maxStagnationIterations, from params, default 3).
        // resetPheromonesOnDeactivationCount controls pheromone resets on deactivation:
        // 0 or negative = never reset, 1 = after the 1st deactivation, 2 = after the 1st
        // and 2nd, etc. Passed as the last (optional) Colony constructor argument; if
        // params.intensifierEnabled is false, nullptr is passed (no extra LS on global best).
        std::shared_ptr<aco::IIntensifier> intensifier;
        if (params.intensifierEnabled)
            intensifier = std::make_shared<aco::Intensifier>(ls, dpCached, inst,
                                                            params.intensifierMaxStagnationIterations,
                                                            params.intensifierResetPheromonesCount);

        auto colony = std::make_unique<aco::Colony>(
            inst, params, runIdx, seed,
            dpCached,
            antPolicy,
            pher,
            cost,
            ls,
            lkPolish,
            threeOptPolish,
            intensifier
        );

        colony->setRecorder(recorder);
        if (logger) colony->setLogger(logger);
        if (researchLogger) colony->setResearchLogger(researchLogger);
        if (researchSurrogateLogger) colony->setResearchSurrogateLogger(researchSurrogateLogger);

        return colony;
    });

    // Final cost evaluation for all runs: re-evaluate cost for the final solutions to
    // ensure consistency before printing. Only re-evaluates the existing car assignment,
    // does not change it.
    {
        // CostModelCars only evaluates cost for an existing Solution (nodes + cars)
        auto costModel = std::make_shared<aco::CostModelCars>(inst);

        // Re-evaluate cost for each run's global best
        const int N = inst->n();
        for (auto& result : results)
        {
            if (std::isfinite(result.best.cost) && !result.best.sol.node.empty())
            {
                // Validate size before evaluating (guard against node.size != N)
                if ((int)result.best.sol.node.size() != N)
                {
                    std::cerr << "[WARN] Run " << result.runIndex << " final eval: node.size="
                              << result.best.sol.node.size() << " != N=" << N << " - skipping eval\n";
                    continue;
                }
                if ((int)result.best.sol.car.size() != N)
                    result.best.sol.car.assign((size_t)N, 0);
                
                try {
                    const double finalCost = costModel->evaluate(result.best.sol);
                    if (std::isfinite(finalCost))
                        result.best.cost = finalCost;
                } catch (const std::exception& e) {
                    std::cerr << "[WARN] Run " << result.runIndex << " final eval failed: " << e.what() << "\n";
                }
            }
        }
    }
    
    // Find the best run before printing (after final evaluation)
    auto bestIt = results.end();
    if (!results.empty()) {
        bestIt = std::min_element(results.begin(), results.end(),
            [](const aco::RunResult& a, const aco::RunResult& b) {
                return a.best.cost < b.best.cost;
            });
    }
    
    // Output results
    if (logger) {
        logger->log("\n=== Run Results ===");
    } else {
        std::cout << "\n=== Run Results ===\n";
    }
    
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& result = results[i];
        const double elapsedMs = recorder->getElapsedTimeForRun(result.runIndex);
        const double elapsedSec = elapsedMs / 1000.0;
        
        std::ostringstream runMsg;
        runMsg << "Run " << result.runIndex << " (seed: " << result.seed << "): cost = " 
               << std::fixed << std::setprecision(3) << result.best.cost
               << ", time = " << std::setprecision(2) << elapsedSec << "s";
        
        if (logger) {
            logger->log(runMsg.str());
        } else {
            std::cout << runMsg.str() << "\n";
        }
        
        // Only print the solution if this isn't the best run (it's printed separately below)
        if (bestIt == results.end() || &result != &(*bestIt)) {
            if (std::isfinite(result.best.cost) && !result.best.sol.node.empty()) {
                std::ostringstream solMsg;
                solMsg << "  Solution nodes: [";
                const auto& nodes = result.best.sol.node;
                const size_t maxShow = std::min(nodes.size(), size_t(20));
                for (size_t j = 0; j < maxShow; ++j) {
                    if (j > 0) solMsg << ", ";
                    solMsg << nodes[j];
                }
                if (nodes.size() > 20) solMsg << ", ... (" << nodes.size() << " nodes)";
                solMsg << "]";
                
                if (logger) {
                    logger->log(solMsg.str());
                } else {
                    std::cout << solMsg.str() << "\n";
                }
                
                // Print cars if present and not all zero
                if (!result.best.sol.car.empty()) {
                    bool hasNonZeroCars = false;
                    for (int c : result.best.sol.car) {
                        if (c != 0) {
                            hasNonZeroCars = true;
                            break;
                        }
                    }
                    if (hasNonZeroCars) {
                        std::ostringstream carMsg;
                        carMsg << "  Solution cars: [";
                        const auto& cars = result.best.sol.car;
                        const size_t maxShowCars = std::min(cars.size(), size_t(20));
                        for (size_t j = 0; j < maxShowCars; ++j) {
                            if (j > 0) carMsg << ", ";
                            carMsg << cars[j];
                        }
                        if (cars.size() > 20) carMsg << ", ... (" << cars.size() << " cars)";
                        carMsg << "]";
                        
                        if (logger) {
                            logger->log(carMsg.str());
                        } else {
                            std::cout << carMsg.str() << "\n";
                        }
                    }
                }
            }
        }
    }
    
    // Print the best run separately, if any
    if (bestIt != results.end()) {
        if (logger) {
            logger->log("\n=== Best Run ===");
        } else {
            std::cout << "\n=== Best Run ===\n";
        }
        
        const double elapsedMs = recorder->getElapsedTimeForRun(bestIt->runIndex);
        const double elapsedSec = elapsedMs / 1000.0;
        
        std::ostringstream bestMsg;
        bestMsg << "Run " << bestIt->runIndex << " (seed: " << bestIt->seed << "): cost = " 
                << std::fixed << std::setprecision(3) << bestIt->best.cost
                << ", time = " << std::setprecision(2) << elapsedSec << "s";
        
        if (logger) {
            logger->log(bestMsg.str());
        } else {
            std::cout << bestMsg.str() << "\n";
        }
        
        // Print the best solution (nodes and cars only once)
        if (std::isfinite(bestIt->best.cost) && !bestIt->best.sol.node.empty()) {
            std::ostringstream bestSolMsg;
            bestSolMsg << "  Best Solution nodes: [";
            const auto& nodes = bestIt->best.sol.node;
            const size_t maxShow = std::min(nodes.size(), size_t(20));
            for (size_t j = 0; j < maxShow; ++j) {
                if (j > 0) bestSolMsg << ", ";
                bestSolMsg << nodes[j];
            }
            if (nodes.size() > 20) bestSolMsg << ", ... (" << nodes.size() << " nodes)";
            bestSolMsg << "]";
            
            if (logger) {
                logger->log(bestSolMsg.str());
            } else {
                std::cout << bestSolMsg.str() << "\n";
            }
            
            // Print cars if present and not all zero
            if (!bestIt->best.sol.car.empty()) {
                bool hasNonZeroCars = false;
                for (int c : bestIt->best.sol.car) {
                    if (c != 0) {
                        hasNonZeroCars = true;
                        break;
                    }
                }
                if (hasNonZeroCars) {
                    std::ostringstream bestCarMsg;
                    bestCarMsg << "  Best Solution cars: [";
                    const auto& cars = bestIt->best.sol.car;
                    const size_t maxShowCars = std::min(cars.size(), size_t(20));
                    for (size_t j = 0; j < maxShowCars; ++j) {
                        if (j > 0) bestCarMsg << ", ";
                        bestCarMsg << cars[j];
                    }
                    if (cars.size() > 20) bestCarMsg << ", ... (" << cars.size() << " cars)";
                    bestCarMsg << "]";
                    
                    if (logger) {
                        logger->log(bestCarMsg.str());
                    } else {
                        std::cout << bestCarMsg.str() << "\n";
                    }
                }
            }
        }
    }
    
    if (logger) {
        logger->log("");
    } else {
        std::cout << "\n";
    }
    
    // Write CSV files
    if (logger) {
        logger->log("Writing CSV files...");
    } else {
        std::cout << "Writing CSV files...\n";
    }
    
    try {
        recorder->writeIterationBestCsv("outputData/iteration_best.csv");
        recorder->writeGlobalBestCsv("outputData/global_best.csv");
        recorder->writeWeibullCsv("outputData/weibulls.csv");
        recorder->writeWeibullWaitingTimeCsv("outputData/weibull_waiting_times.csv");
        recorder->writeWeibullHypothesesCsv("outputData/weibull_hypotheses.csv");
        recorder->writeIntensifierEventsCsv("outputData/intensifier_events.csv");
        recorder->writeAntsQuantilesCsv("outputData/ants_quantiles.csv");
        recorder->writeRunSummaryCsv("outputData/runs_summary.csv");
        recorder->writeDynamicParamsCsv("outputData/dynamic_params.csv");
        recorder->writeSolutionsCsv("outputData/solutions.csv", results);

        if (researchLogger && researchLogger->numRunsWithData() > 0)
            researchLogger->writeCsv("outputData/intensifier_binomial.csv");
        if (researchSurrogateLogger && researchSurrogateLogger->numRecords() > 0)
            researchSurrogateLogger->writeCsv("outputData/surrogate_correlation.csv");

        if (logger) {
            logger->log("CSV files written successfully.");
        } else {
            std::cout << "CSV files written successfully.\n";
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error writing CSV files: " << e.what() << "\n";
        if (logger) {
            logger->log("Error writing CSV files: " + std::string(e.what()));
        }
    }
}

int main(int argc, char* argv[])
{
    try {
        aco_cli::Parser parser;
        auto params = parser.parse(argc, argv);

        cars_tsplib::Parser parserInst;
        auto inst = std::make_shared<cars_tsplib::Instance const>(parserInst.parseFileSpec(params.filename));
        
        // Build the log file path from the problem name
        std::string logFilePath = "";
        if (!params.filename.empty()) {
            // Extract the problem name (strip path and extension)
            std::string problemName = params.filename;
            size_t lastSlash = problemName.find_last_of("/\\");
            if (lastSlash != std::string::npos) {
                problemName = problemName.substr(lastSlash + 1);
            }
            size_t lastDot = problemName.find_last_of(".");
            if (lastDot != std::string::npos) {
                problemName = problemName.substr(0, lastDot);
            }
            logFilePath = "logs/" + problemName + ".log";
        }
        
        auto logger = std::make_shared<aco::Logger>(logFilePath);

        runAco(params, inst, logger);
    }
    catch (const aco_cli::ParseError& e) {
        if (std::string(e.what()) == "HELP_REQUESTED") {
            aco_cli::Parser parser;
            parser.printUsage(std::cout, argv[0]);
            return 0;
        }
        std::cerr << "Parse error: " << e.what() << "\n";
        return 1;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
