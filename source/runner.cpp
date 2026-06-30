#include <algorithm>
#include <iostream>
#include <iomanip>
#include <params_cli.hpp>
#include <runner.hpp>
#include <colony.hpp>

namespace aco 
{

Runner::Runner(RunnerConfig cfg) : cfg_(cfg) {}

std::vector<RunResult> Runner::runAll(const aco_cli::ParamsData& params, uint64_t baseSeed, ColonyFactory factory) 
{
    const int nRuns = params.nRuns;
    std::vector<RunResult> results;
    results.resize(static_cast<size_t>(nRuns));

    // Tracks which run currently owns verbose logging
    std::atomic<int> verboseOwner{-1};
    std::mutex printMx;

    // Progress tracking per run (current iteration)
    std::vector<std::atomic<int>> currentIteration(static_cast<size_t>(nRuns));
    for (auto& iter : currentIteration) iter.store(-1, std::memory_order_relaxed);

    // Run status: -1 = not started, 0 = in progress, 1 = finished
    std::vector<std::atomic<int>> runStatus(static_cast<size_t>(nRuns));
    for (auto& status : runStatus) status.store(-1, std::memory_order_relaxed);

    std::atomic<int> next{0};

    int nThreads = cfg_.threads;
    if (nThreads <= 0) {
        nThreads = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
    }
    nThreads = std::min(nThreads, nRuns);

    // Selects the next run to make verbose (the one with the most iterations remaining, i.e. will run longest)
    auto selectNextVerbose = [&]() -> int {
        int bestIdx = -1;
        int maxIter = -1;

        for (int i = 0; i < nRuns; ++i) {
            int status = runStatus[i].load(std::memory_order_acquire);
            if (status == 0) { // in progress
                int iter = currentIteration[i].load(std::memory_order_acquire);
                if (iter >= 0) {
                    int remaining = params.iterations - iter;
                    if (remaining > maxIter) {
                        maxIter = remaining;
                        bestIdx = i;
                    }
                }
            }
        }
        return bestIdx;
    };

    // Stores colony references so the verbose flag can be updated while running.
    // Uses weak_ptr since ownership stays with the worker (colony is destroyed after run()).
    std::vector<std::weak_ptr<Colony>> colonies(static_cast<size_t>(nRuns));
    std::mutex coloniesMx; // guards thread-safe access to the colonies vector

    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(nThreads));

    for (int t = 0; t < nThreads; ++t) {
        workers.emplace_back([&, t]() {
            (void)t;
            while (true) {
                int idx = next.fetch_add(1, std::memory_order_relaxed);
                if (idx >= nRuns) break;

                const uint64_t colonySeed = aco::mixSeedRun(baseSeed, (uint64_t)idx);

                // Mark the run as started
                runStatus[idx].store(0, std::memory_order_release);

                // Try to become the "verbose" run if none is currently set
                bool iAmVerbose = false;
                {
                    int expected = -1;
                    if (verboseOwner.compare_exchange_strong(expected, idx, std::memory_order_acq_rel)) {
                        iAmVerbose = true;
                    }
                }

                auto colony = factory(idx, colonySeed);

                // Keep a shared_ptr reference to the colony for later verbose-flag updates
                std::shared_ptr<Colony> colonyShared = std::shared_ptr<Colony>(colony.release());
                {
                    std::lock_guard<std::mutex> lock(coloniesMx);
                    colonies[static_cast<size_t>(idx)] = colonyShared;
                }

                // Set the colony's logger/verbose flag
                colonyShared->setVerbose(iAmVerbose, &printMx);

                // Callback to update progress tracking
                colonyShared->setProgressCallback([&, idx](int iteration) {
                    currentIteration[idx].store(iteration, std::memory_order_release);

                    // Try to become the verbose run, but only if there isn't one already
                    int currentVerbose = verboseOwner.load(std::memory_order_acquire);
                    if (currentVerbose == -1) {
                        int expected = -1;
                        if (verboseOwner.compare_exchange_strong(expected, idx, std::memory_order_acq_rel)) {
                            {
                                std::lock_guard<std::mutex> lock(coloniesMx);
                                auto idxColony = colonies[static_cast<size_t>(idx)].lock();
                                if (idxColony) {
                                    idxColony->updateVerbose(true);
                                }
                            }
                            {
                                std::lock_guard<std::mutex> lock(printMx);
                                std::cout << "\n[Switching verbose] Switching to run " << idx 
                                          << " (most iterations remaining: " << (params.iterations - iteration) << ").\n";
                            }
                        }
                    }
                    // The verbose run is only reassigned when a run finishes, not mid-run
                });

                results[static_cast<size_t>(idx)] = colonyShared->run();

                // Mark the run as finished
                runStatus[idx].store(1, std::memory_order_release);

                // Hand off verbose logging to another run if this one owned it
                if (iAmVerbose)
                {
                    int expected = idx;
                    verboseOwner.compare_exchange_strong(expected, -1, std::memory_order_acq_rel);

                    int nextVerbose = selectNextVerbose();
                    if (nextVerbose >= 0) {
                        int expected2 = -1;
                        if (verboseOwner.compare_exchange_strong(expected2, nextVerbose, std::memory_order_acq_rel)) {
                            {
                                std::lock_guard<std::mutex> lock(coloniesMx);
                                auto nextVerboseColony = colonies[static_cast<size_t>(nextVerbose)].lock();
                                if (nextVerboseColony) {
                                    nextVerboseColony->updateVerbose(true);
                                }
                            }
                            {
                                std::lock_guard<std::mutex> lock(printMx);
                                int nextIter = currentIteration[nextVerbose].load(std::memory_order_acquire);
                                int remaining = (nextIter >= 0) ? (params.iterations - nextIter) : params.iterations;
                                std::cout << "\n[Switching verbose] Run " << idx << " completed. Switching to run " 
                                          << nextVerbose << " (most iterations remaining: " << remaining << ").\n";
                            }
                        }
                    }
                }
            }
        });
    }

    for (auto& th : workers) th.join();
    return results;
}

} // namespace aco
