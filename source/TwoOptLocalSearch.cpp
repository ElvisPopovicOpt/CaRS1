#include <algorithm>
#include <cmath>
#include <chrono>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <parser.hpp>
#include <interfaces.hpp>

#include "TwoOptLocalSearch.hpp"
#include "parser.hpp"          // cars_tsplib::Instance definition
#include "candidate_list.hpp"  // aco::CandidateListCache

#include "solution.hpp"                // aco::Solution
#include "FixedTourCarAssignerDP.hpp"  // FixedTourCarAssignerDP

namespace localSearch 
{

// --------------------------
// Utilities (cpp-only)
// --------------------------

static std::uint64_t nowNs()
{
    using clock = std::chrono::steady_clock;
    return (std::uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(
               clock::now().time_since_epoch())
        .count();
}

static inline bool isValid2OptMove(int i, int k, int N)
{
    // fixed start: i must begin at 1
    if (i <= 0) return false;
    if (k <= i) return false;
    if (k >= N) return false;
    if (k == i + 1) return false; // trivial move (reverse of length 2)
    return true;
}

static inline void applyReverseAndUpdatePos(std::vector<int>& tour,
                                           std::vector<int>& pos,
                                           int i,
                                           int k)
{
    std::reverse(tour.begin() + i, tour.begin() + (k + 1));
    for (int t = i; t <= k; ++t) {
        pos[(std::size_t)tour[(std::size_t)t]] = t;
    }
}

// Virtual reverse view: doesn't mutate the base tour
struct TwoOptView
{
    const std::vector<int>& base;
    int i, k; // reversed range [i..k]

    int operator()(int idx) const
    {
        if (idx < i || idx > k) return base[(std::size_t)idx];
        return base[(std::size_t)(i + (k - idx))];
    }
};

// Symmetric TSP 2-opt delta (only 2 edges change)
static inline double tsp2OptDeltaSym(const cars_tsplib::Instance& inst,
                                     const std::vector<int>& tour,
                                     int i,
                                     int k)
{
    const int N = (int)tour.size();
    const int a = tour[(std::size_t)(i - 1)];
    const int b = tour[(std::size_t)i];
    const int c = tour[(std::size_t)k];
    const int d = tour[(std::size_t)((k + 1) % N)];

    return inst.travelCost(0, a, c) + inst.travelCost(0, b, d) - inst.travelCost(0, a, b) -
           inst.travelCost(0, c, d);
}

// --------------------------
// TwoOptLocalSearch
// --------------------------

TwoOptLocalSearch::TwoOptLocalSearch(std::shared_ptr<const cars_tsplib::Instance> inst,
                                     TwoOptOptions opt,
                                     std::shared_ptr<FixedTourCarAssignerDP> dp)
    : inst_(std::move(inst))
    , opt_(std::move(opt))
    , dp_(std::move(dp))
{
    if (!dp_ && inst_ && inst_->cars() > 1) {
        dp_ = std::make_shared<FixedTourCarAssignerDP>(inst_, opt_.dpOpt);
    }
}

double TwoOptLocalSearch::singleCarCycleCost_(const std::vector<int>& nodes) const
{
    const int N = (int)nodes.size();
    double c = 0.0;
    for (int i = 0; i < N; ++i) {
        const int u = nodes[(std::size_t)i];
        const int v = nodes[(std::size_t)((i + 1) % N)];
        c += inst_->travelCost(0, u, v);
    }
    return c;
}

bool TwoOptLocalSearch::timeExceeded_(std::uint64_t startNs) const
{
    if (opt_.timeLimitMs <= 0.0) return false;
    const double elapsedMs = (double)(nowNs() - startNs) / 1e6;
    return elapsedMs >= opt_.timeLimitMs;
}

bool TwoOptLocalSearch::likelySymmetricTsp_() const
{
    if (!inst_) return true;
    const int N = inst_->n();
    if (N <= 2) return true;

    // Heuristic symmetry check for car=0 (TSP fallback branch), using a deterministic sample (no RNG).
    const int M = std::min(N, 25);
    for (int i = 0; i < M; ++i) {
        for (int j = i + 1; j < M; ++j) {
            const double a = inst_->travelCost(0, i, j);
            const double b = inst_->travelCost(0, j, i);
            if (std::isfinite(a) && std::isfinite(b)) {
                if (std::fabs(a - b) > 1e-12) return false;
            }
        }
    }
    return true;
}

void TwoOptLocalSearch::improve(aco::Solution& s, double& cost) const
{
    const int N = inst_->n();
    const bool useCarDP = (dp_ != nullptr);

    if ((int)s.node.size() != N)
        throw std::runtime_error("TwoOptLocalSearch: Solution node size mismatch with instance N.");
    if ((int)s.car.size() != N) s.car.assign((std::size_t)N, 0);
    if (N < 4) return;
    if (s.node[0] != 0) throw std::runtime_error("TwoOptLocalSearch: node[0] must be 0.");

    // baseline
    if (opt_.reassignCarsAtStart) {
        if (useCarDP) {
            cost = dp_->reassignCars(s.node, s.car);
        } else {
            std::fill(s.car.begin(), s.car.end(), 0);
            cost = singleCarCycleCost_(s.node);
        }
    } else {
        // Not reassigning at start; for the TSP fallback at least clear car. Cost is assumed valid from the caller.
        if (!useCarDP) std::fill(s.car.begin(), s.car.end(), 0);
    }

    if (!std::isfinite(cost))
        throw std::runtime_error("TwoOptLocalSearch: initial cost is not finite.");

    // pos[node] = index
    std::vector<int> pos((std::size_t)N, -1);
    for (int idx = 0; idx < N; ++idx) pos[(std::size_t)s.node[(std::size_t)idx]] = idx;

    std::int64_t evals = 0;
    const std::uint64_t t0 = nowNs();

    const bool useCand = (opt_.cand != nullptr) && (opt_.cand->N() == N) && (opt_.cand->K() > 0);

    // DLB keyed by node id (stable across reverses)
    std::vector<std::uint8_t> dontLookNode;
    if (opt_.useDontLookBits) dontLookNode.assign((std::size_t)N, 0);

    // TSP fallback: use the symmetric delta only if the heuristic says the instance is symmetric.
    const bool useSymTspDelta = (!useCarDP) ? likelySymmetricTsp_() : false;

    bool anyAccepted = false;

    if (opt_.firstImprovement) {
        improveFirst_(s, cost, useCarDP, useCand, useSymTspDelta, pos, dontLookNode, evals, t0,
                      anyAccepted);
    } else {
        improveBest_(s, cost, useCarDP, useCand, useSymTspDelta, pos, dontLookNode, evals, t0,
                     anyAccepted);
    }

    // At the end, reconcile cars only once (if anything changed and DP is in use)
    if (useCarDP && anyAccepted) {
        cost = dp_->reassignCars(s.node, s.car);
        if (!std::isfinite(cost))
            throw std::runtime_error(
                "TwoOptLocalSearch: final reassignCars produced non-finite cost.");
    } else if (!useCarDP) {
        std::fill(s.car.begin(), s.car.end(), 0);
    }
}

// --------------------------
// First-improvement
// --------------------------

void TwoOptLocalSearch::improveFirst_(aco::Solution& s,
                                      double& cost,
                                      bool useCarDP,
                                      bool useCand,
                                      bool useSymTspDelta,
                                      std::vector<int>& pos,
                                      std::vector<std::uint8_t>& dontLookNode,
                                      std::int64_t& evals,
                                      std::uint64_t t0,
                                      bool& anyAccepted) const
{
    const int N = (int)s.node.size();
    const bool useDLB = opt_.useDontLookBits;

    auto stopRequested = [&]() -> bool {
        if (opt_.maxMoveEvaluations > 0 && evals >= opt_.maxMoveEvaluations) return true;
        if (timeExceeded_(t0)) return true;
        return false;
    };

    auto wakePos = [&](int posIdx) {
        if (!useDLB) return;
        posIdx %= N;
        if (posIdx < 0) posIdx += N;
        const int v = s.node[(std::size_t)posIdx];
        dontLookNode[(std::size_t)v] = 0;
    };

    auto wakeAroundCuts = [&](int i, int k) {
        if (!useDLB) return;
        const int r = std::max(0, opt_.dontLookWakeRadius);
        for (int d = -r; d <= r; ++d) {
            wakePos((i - 1) + d);
            wakePos(i + d);
            wakePos(k + d);
            wakePos((k + 1) + d);
        }
    };

    // Stamp-based dedup for k (k are positions 0..N-1)
    std::vector<std::uint32_t> seenStamp((std::size_t)N, 0u);
    std::uint32_t stamp = 1u;

    // Early termination: stop after maxNoImprove consecutive passes with no improvement
    int noImproveCount = 0;
    const int maxNoImprove = std::max(10, N / 4); // scales with instance size

    for (int pass = 0; pass < opt_.maxPasses; ++pass) {
        bool improvedThisPass = false;
        if (stopRequested()) break;

        for (int i = 1; i <= N - 2; ++i) {
            if (stopRequested()) return;

            const int a = s.node[(std::size_t)(i - 1)];
            const int b = s.node[(std::size_t)i];

            // never let DLB block the depot (0)
            if (useDLB && a != 0 && b != 0 && dontLookNode[(std::size_t)a] &&
                dontLookNode[(std::size_t)b]) {
                continue;
            }

            bool acceptedForThisI = false;

            std::vector<int> ks;

            auto pushK = [&](int k) {
                if (!isValid2OptMove(i, k, N)) return;
                if (seenStamp[(std::size_t)k] == stamp) return;
                seenStamp[(std::size_t)k] = stamp;
                ks.push_back(k);
            };

            ++stamp;
            if (stamp == 0u) {
                std::fill(seenStamp.begin(), seenStamp.end(), 0u);
                stamp = 1u;
            }

            if (useCand) {
                auto scanAnchor = [&](int anchorNode) {
                    const auto& cand = opt_.cand->candidates(anchorNode);
                    for (int j : cand) {
                        const int k = pos[(std::size_t)j];
                        pushK(k);
                    }
                };
                scanAnchor(a);
                scanAnchor(b);
            } else {
                for (int k = i + 1; k <= N - 1; ++k) pushK(k);
            }

            // try ks in order (FIRST improvement)
            for (int k : ks) {
                if (stopRequested()) return;

                double candCost = 0.0;

                if (useCarDP) {
                    TwoOptView view{s.node, i, k};
                    candCost = dp_->evaluateCostViewScratch([&](int idx) { return view(idx); }, scratch_);
                } else {
                    if (useSymTspDelta) {
                        candCost = cost + tsp2OptDeltaSym(*inst_, s.node, i, k);
                    } else {
                        // ATSP-safe: full O(N) evaluation over the virtual tour
                        TwoOptView view{s.node, i, k};
                        candCost = 0.0;
                        for (int t = 0; t < N; ++t) {
                            const int u = view(t);
                            const int v = view((t + 1) % N);
                            candCost += inst_->travelCost(0, u, v);
                        }
                    }
                }

                ++evals;
                if (!std::isfinite(candCost)) continue;

                if (candCost + opt_.eps < cost) {
                    // ACCEPT
                    applyReverseAndUpdatePos(s.node, pos, i, k);
                    wakeAroundCuts(i, k);

                    cost = candCost;

                    if (!useCarDP) {
                        std::fill(s.car.begin(), s.car.end(), 0);
                    } else {
                        anyAccepted = true; // cars are stale; reassigned once at the end of improve()
                    }

                    acceptedForThisI = true;
                    improvedThisPass = true;
                    break; // first improvement
                }
            }

            if (!acceptedForThisI) {
                if (useDLB) {
                    if (a != 0) dontLookNode[(std::size_t)a] = 1;
                    if (b != 0) dontLookNode[(std::size_t)b] = 1;
                }
            } else {
                // first improvement: restart the pass as soon as a move is accepted
                break;
            }
        }

        // Early termination: stop after maxNoImprove consecutive passes with no improvement
        if (improvedThisPass) {
            noImproveCount = 0;
        } else {
            noImproveCount++;
            if (noImproveCount >= maxNoImprove) {
                break;
            }
        }
    }
}

// --------------------------
// Best-improvement
// --------------------------

void TwoOptLocalSearch::improveBest_(aco::Solution& s,
                                     double& cost,
                                     bool useCarDP,
                                     bool useCand,
                                     bool useSymTspDelta,
                                     std::vector<int>& pos,
                                     std::vector<std::uint8_t>& dontLookNode,
                                     std::int64_t& evals,
                                     std::uint64_t t0,
                                     bool& anyAccepted) const
{
    const int N = (int)s.node.size();
    const bool useDLB = opt_.useDontLookBits;

    auto stopRequested = [&]() -> bool {
        if (opt_.maxMoveEvaluations > 0 && evals >= opt_.maxMoveEvaluations) return true;
        if (timeExceeded_(t0)) return true;
        return false;
    };

    auto wakePos = [&](int posIdx) {
        if (!useDLB) return;
        posIdx %= N;
        if (posIdx < 0) posIdx += N;
        const int v = s.node[(std::size_t)posIdx];
        dontLookNode[(std::size_t)v] = 0;
    };

    auto wakeAroundCuts = [&](int i, int k) {
        if (!useDLB) return;
        const int r = std::max(0, opt_.dontLookWakeRadius);
        for (int d = -r; d <= r; ++d) {
            wakePos((i - 1) + d);
            wakePos(i + d);
            wakePos(k + d);
            wakePos((k + 1) + d);
        }
    };

    // Stamp-based dedup
    std::vector<std::uint32_t> seenStamp((std::size_t)N, 0u);
    std::uint32_t stamp = 1u;

    // Early termination: stop after maxNoImprove consecutive passes with no improvement
    int noImproveCount = 0;
    const int maxNoImprove = std::max(10, N / 4); // scales with instance size

    for (int pass = 0; pass < opt_.maxPasses; ++pass) {
        if (stopRequested()) break;

        double bestCost = cost;
        int bestI = -1, bestK = -1;

        for (int i = 1; i <= N - 2; ++i) {
            if (stopRequested()) return;

            const int a = s.node[(std::size_t)(i - 1)];
            const int b = s.node[(std::size_t)i];

            // never let DLB block the depot (0)
            if (useDLB && a != 0 && b != 0 && dontLookNode[(std::size_t)a] &&
                dontLookNode[(std::size_t)b]) {
                continue;
            }

            // Standard DLB criterion for best-improvement mode: only sleep (a,b) if no move improves the current cost.
            bool hasImprovingMoveForThisI = false;

            std::vector<int> ks;

            auto pushK = [&](int k) {
                if (!isValid2OptMove(i, k, N)) return;
                if (seenStamp[(std::size_t)k] == stamp) return;
                seenStamp[(std::size_t)k] = stamp;
                ks.push_back(k);
            };

            ++stamp;
            if (stamp == 0u) {
                std::fill(seenStamp.begin(), seenStamp.end(), 0u);
                stamp = 1u;
            }

            if (useCand) {
                auto scanAnchor = [&](int anchorNode) {
                    const auto& cand = opt_.cand->candidates(anchorNode);
                    for (int j : cand) {
                        const int k = pos[(std::size_t)j];
                        pushK(k);
                    }
                };
                scanAnchor(a);
                scanAnchor(b);
            } else {
                for (int k = i + 1; k <= N - 1; ++k) pushK(k);
            }

            for (int k : ks) {
                if (stopRequested()) return;

                double candCost = 0.0;

                if (useCarDP) {
                    TwoOptView view{s.node, i, k};
                    candCost = dp_->evaluateCostViewScratch([&](int idx) { return view(idx); }, scratch_);
                } else {
                    if (useSymTspDelta) {
                        candCost = cost + tsp2OptDeltaSym(*inst_, s.node, i, k);
                    } else {
                        // ATSP-safe: full O(N) evaluation over the virtual tour
                        TwoOptView view{s.node, i, k};
                        candCost = 0.0;
                        for (int t = 0; t < N; ++t) {
                            const int u = view(t);
                            const int v = view((t + 1) % N);
                            candCost += inst_->travelCost(0, u, v);
                        }
                    }
                }

                ++evals;
                if (!std::isfinite(candCost)) continue;

                if (candCost + opt_.eps < cost) {
                    hasImprovingMoveForThisI = true;
                }
                if (candCost + opt_.eps < bestCost) {
                    bestCost = candCost;
                    bestI = i;
                    bestK = k;
                }
            }

            if (!hasImprovingMoveForThisI) {
                if (useDLB) {
                    if (a != 0) dontLookNode[(std::size_t)a] = 1;
                    if (b != 0) dontLookNode[(std::size_t)b] = 1;
                }
            }
        }

        // Early termination: stop after maxNoImprove consecutive passes with no improvement
        if (bestI < 0) {
            noImproveCount++;
            if (noImproveCount >= maxNoImprove) {
                break;
            }
        } else {
            noImproveCount = 0;
        }

        if (bestI < 0) break; // no improving move found this pass

        // ACCEPT best
        applyReverseAndUpdatePos(s.node, pos, bestI, bestK);
        wakeAroundCuts(bestI, bestK);

        cost = bestCost;

        if (!useCarDP) {
            std::fill(s.car.begin(), s.car.end(), 0);
        } else {
            anyAccepted = true; // cars are stale; reassigned once at the end of improve()
        }
    }
}

} // namespace localSearch

