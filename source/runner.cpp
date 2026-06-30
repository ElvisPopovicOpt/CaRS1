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

    // za ispis 
    std::atomic<int> verboseOwner{-1};
    std::mutex printMx;

    // Praćenje napretka svakog runa (trenutna iteracija)
    std::vector<std::atomic<int>> currentIteration(static_cast<size_t>(nRuns));
    for (auto& iter : currentIteration) iter.store(-1, std::memory_order_relaxed);
    
    // Status runova: -1 = nije počeo, 0 = u tijeku, 1 = završen
    std::vector<std::atomic<int>> runStatus(static_cast<size_t>(nRuns));
    for (auto& status : runStatus) status.store(-1, std::memory_order_relaxed);

    std::atomic<int> next{0};

    int nThreads = cfg_.threads;
    if (nThreads <= 0) {
        nThreads = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
    }
    nThreads = std::min(nThreads, nRuns);

    // Funkcija za odabir sljedećeg verbose runa (najviše iteracija ostalo - najdulje će raditi)
    auto selectNextVerbose = [&]() -> int {
        int bestIdx = -1;
        int maxIter = -1;
        
        for (int i = 0; i < nRuns; ++i) {
            int status = runStatus[i].load(std::memory_order_acquire);
            if (status == 0) { // u tijeku
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

    // Pohrana colony referenci za ažuriranje verbose flag-a tijekom izvršavanja
    // Koristimo weak_ptr jer ne želimo zadržati ownership (colony se uništi nakon run())
    std::vector<std::weak_ptr<Colony>> colonies(static_cast<size_t>(nRuns));
    std::mutex coloniesMx; // za thread-safe pristup colonies vektoru

    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(nThreads));

    for (int t = 0; t < nThreads; ++t) {
        workers.emplace_back([&, t]() {
            (void)t;
            while (true) {
                int idx = next.fetch_add(1, std::memory_order_relaxed);
                if (idx >= nRuns) break;

                const uint64_t colonySeed = aco::mixSeedRun(baseSeed, (uint64_t)idx);

                // Označi da je run počeo
                runStatus[idx].store(0, std::memory_order_release);

                // za ispis
                // pokušaj postati "verbose" run ako ga trenutno nema
                bool iAmVerbose = false;
                {
                    int expected = -1;
                    if (verboseOwner.compare_exchange_strong(expected, idx, std::memory_order_acq_rel)) {
                        iAmVerbose = true;
                    }
                }

                auto colony = factory(idx, colonySeed);
                
                // Pohrani shared_ptr referencu na colony za kasnije ažuriranje verbose flag-a
                std::shared_ptr<Colony> colonyShared = std::shared_ptr<Colony>(colony.release());
                {
                    std::lock_guard<std::mutex> lock(coloniesMx);
                    colonies[static_cast<size_t>(idx)] = colonyShared;
                }
                
                // ispis: postavi colonyju logger/flag
                colonyShared->setVerbose(iAmVerbose, &printMx);
                
                // Callback za ažuriranje napretka
                colonyShared->setProgressCallback([&, idx](int iteration) {
                    currentIteration[idx].store(iteration, std::memory_order_release);
                    
                    // Provjeri treba li postati verbose run (samo ako nema verbose runa)
                    int currentVerbose = verboseOwner.load(std::memory_order_acquire);
                    if (currentVerbose == -1) {
                        // Nema verbose runa, pokušaj postati verbose
                        int expected = -1;
                        if (verboseOwner.compare_exchange_strong(expected, idx, std::memory_order_acq_rel)) {
                            // Postani verbose
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
                    // Ne prebacujemo verbose run tijekom izvršavanja - samo kada završi
                });

                results[static_cast<size_t>(idx)] = colonyShared->run();
                
                // Označi da je run završen
                runStatus[idx].store(1, std::memory_order_release);
                
                // za ispis nakon zavrsetka
                if (iAmVerbose) 
                {
                    int expected = idx;
                    verboseOwner.compare_exchange_strong(expected, -1, std::memory_order_acq_rel);
                    
                    // Odaberi sljedeći run za ispis (najmanje iteracija)
                    int nextVerbose = selectNextVerbose();
                    if (nextVerbose >= 0) {
                        int expected2 = -1;
                        if (verboseOwner.compare_exchange_strong(expected2, nextVerbose, std::memory_order_acq_rel)) {
                            // Ažuriraj verbose flag za novi run
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
