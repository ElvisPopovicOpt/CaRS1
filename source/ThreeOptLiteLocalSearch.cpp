#include "ThreeOptLiteLocalSearch.hpp"

#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "parser.hpp"
#include "candidate_list.hpp"
#include "FixedTourCarAssignerDP.hpp"

namespace localSearch
{

ThreeOptLiteLocalSearch::ThreeOptLiteLocalSearch(std::shared_ptr<const cars_tsplib::Instance> inst,
                                                 ThreeOptLiteOptions opt,
                                                 std::shared_ptr<FixedTourCarAssignerDP> dp)
    : inst_(std::move(inst)), opt_(std::move(opt)), dp_(std::move(dp)), rng_(0)
{
    if (!inst_) throw std::runtime_error("ThreeOptLiteLocalSearch: inst is null");
    if (!dp_)   throw std::runtime_error("ThreeOptLiteLocalSearch: dp is null");
    if (!opt_.cand) throw std::runtime_error("ThreeOptLiteLocalSearch: cand is null");

    const int N = inst_->n();
    const int C = inst_->cars();
    if (N <= 0) throw std::runtime_error("ThreeOptLiteLocalSearch: inst->n() <= 0");
    if (C <= 0) throw std::runtime_error("ThreeOptLiteLocalSearch: inst->cars() <= 0");

    // Precompute minTravel (O(C*N^2)) so each later call is O(1).
    minTravel_.assign((std::size_t)N * (std::size_t)N, std::numeric_limits<double>::infinity());

    for (int u = 0; u < N; ++u)
    {
        for (int v = 0; v < N; ++v)
        {
            if (u == v) continue;
            double best = std::numeric_limits<double>::infinity();
            for (int c = 0; c < C; ++c)
                best = std::min(best, inst_->travelCost(c, u, v));
            minTravel_[(std::size_t)u * (std::size_t)N + (std::size_t)v] = best;
        }
    }
}

void ThreeOptLiteLocalSearch::resetSeed(uint64_t seed) const
{
    rng_.reseed(seed);
}

double ThreeOptLiteLocalSearch::minTravelFast_(int u, int v) const
{
    const int N = inst_->n();
    return minTravel_[(std::size_t)u * (std::size_t)N + (std::size_t)v];
}

bool ThreeOptLiteLocalSearch::isValidCut_(int i, int j, int k, int N)
{
    // Keep node[0]==0 fixed: i must be >=1, and we never cut around the 0-edges.
    if (N < 6) return false;
    if (i < 1) return false;
    if (!(i < j && j < k)) return false;
    if (k >= N) return false;

    // Ensure segments S1=[i..j-1] and S2=[j..k-1] are non-empty.
    if (j == i) return false;
    if (k == j) return false;

    // k<N already keeps us from cutting the closing edge (N-1 -> 0),
    // and i>=1 already keeps us from cutting the (0->1) edge.
    return true;
}

double ThreeOptLiteLocalSearch::surrogate3OptDeltaSwap_(const std::vector<int>& tour,
                                                       int i, int j, int k,
                                                       bool revS1, bool revS2) const
{
    // Cutting edges A->B, C->D, E->F:
    // A=tour[i-1], B=tour[i]; C=tour[j-1], D=tour[j]; E=tour[k-1], F=tour[k]
    const int A = tour[(std::size_t)i - 1];
    const int B = tour[(std::size_t)i];
    const int C = tour[(std::size_t)j - 1];
    const int D = tour[(std::size_t)j];
    const int E = tour[(std::size_t)k - 1];
    const int F = tour[(std::size_t)k];

    const double rem = minTravelFast_(A, B) + minTravelFast_(C, D) + minTravelFast_(E, F);

    // start/end for S1=[i..j-1]
    const int start1 = revS1 ? C : B;
    const int end1   = revS1 ? B : C;

    // start/end for S2=[j..k-1]
    const int start2 = revS2 ? E : D;
    const int end2   = revS2 ? D : E;

    // After swap: S0 + S2 + S1 + S3
    const double add =
        minTravelFast_(A, start2) +
        minTravelFast_(end2, start1) +
        minTravelFast_(end1, F);

    return add - rem;
}

void ThreeOptLiteLocalSearch::applySwap_(std::vector<int>& tour,
                                        std::vector<int>& pos,
                                        int i, int j, int k,
                                        bool revS1, bool revS2) const
{
    // tour = S0 + S2 + S1 + S3, with reversal options for S1 and S2.
    const int N = (int)tour.size();
    std::vector<int> tmp;
    tmp.resize((std::size_t)N);

    int p = 0;

    // S0: [0..i-1]
    for (int a = 0; a < i; ++a) tmp[(std::size_t)p++] = tour[(std::size_t)a];

    // S2: [j..k-1]
    if (!revS2)
    {
        for (int a = j; a < k; ++a) tmp[(std::size_t)p++] = tour[(std::size_t)a];
    }
    else
    {
        for (int a = k - 1; a >= j; --a) tmp[(std::size_t)p++] = tour[(std::size_t)a];
    }

    // S1: [i..j-1]
    if (!revS1)
    {
        for (int a = i; a < j; ++a) tmp[(std::size_t)p++] = tour[(std::size_t)a];
    }
    else
    {
        for (int a = j - 1; a >= i; --a) tmp[(std::size_t)p++] = tour[(std::size_t)a];
    }

    // S3: [k..N-1]
    for (int a = k; a < N; ++a) tmp[(std::size_t)p++] = tour[(std::size_t)a];

    tour.swap(tmp);

    // rebuild pos (simple and fast enough for N<=300)
    std::fill(pos.begin(), pos.end(), -1);
    for (int idx = 0; idx < N; ++idx)
        pos[(std::size_t)tour[(std::size_t)idx]] = idx;
}

void ThreeOptLiteLocalSearch::improve(aco::Solution& s, double& cost) const
{
    const int N = inst_->n();
    if (N < 6) return;
    if ((int)s.node.size() != N) throw std::runtime_error("ThreeOptLiteLocalSearch: node size mismatch");
    if (s.node[0] != 0) throw std::runtime_error("ThreeOptLiteLocalSearch: node[0] must be 0");
    if (!std::isfinite(cost)) return;

    if ((int)s.car.size() != N) s.car.assign((std::size_t)N, 0);

    const double startCost = cost;

    aco::Solution bestSol = s;
    double bestCost = cost;

    for (int att = 0; att < opt_.attempts; ++att)
    {
        aco::Solution cand = s;

        std::vector<int> pos((std::size_t)N, -1);
        for (int idx = 0; idx < N; ++idx)
            pos[(std::size_t)cand.node[(std::size_t)idx]] = idx;

        // 3-opt-lite chain (surrogate cost only)
        for (int step = 0; step < opt_.chainLen; ++step)
        {
            int bestI = -1, bestJ = -1, bestK = -1;
            bool bestRev1 = false, bestRev2 = false;
            double bestDelta = std::numeric_limits<double>::infinity();

            for (int t = 0; t < opt_.triesPerStep; ++t)
            {
                const int i = rng_.uniformInt(1, N - 3);

                // Anchors as in LK-lite (try both)
                const int anchors1[2] = {
                    cand.node[(std::size_t)(i - 1)],
                    cand.node[(std::size_t)i]
                };

                for (int a1 = 0; a1 < 2; ++a1)
                {
                    const auto& neighJ = opt_.cand->candidates(anchors1[a1]);
                    int takenJ = 0;

                    for (int v : neighJ)
                    {
                        const int j = pos[(std::size_t)v];
                        // Need i<j<k<=N-1, and i must not be a 0-cut.
                        if (j <= i || j >= N - 1) continue; // j in (i..N-2)

                        // Second anchor around the cut at j
                        const int anchors2[2] = {
                            cand.node[(std::size_t)(j - 1)],
                            cand.node[(std::size_t)j]
                        };

                        for (int a2 = 0; a2 < 2; ++a2)
                        {
                            const auto& neighK = opt_.cand->candidates(anchors2[a2]);
                            int takenK = 0;

                            for (int w : neighK)
                            {
                                const int k = pos[(std::size_t)w];
                                if (k <= j || k >= N) continue; // k in (j..N-1)

                                if (!isValidCut_(i, j, k, N)) continue;

                                // 4 variants of revS1 x revS2
                                for (int mask = 0; mask < 4; ++mask)
                                {
                                    const bool revS1 = (mask & 1) != 0;
                                    const bool revS2 = (mask & 2) != 0;

                                    const double delta = surrogate3OptDeltaSwap_(cand.node, i, j, k, revS1, revS2);
                                    if (delta < bestDelta)
                                    {
                                        bestDelta = delta;
                                        bestI = i; bestJ = j; bestK = k;
                                        bestRev1 = revS1; bestRev2 = revS2;
                                    }
                                }

                                if (++takenK >= opt_.topKPerI) break;
                            }
                        }

                        if (++takenJ >= opt_.topKPerI) break;
                    }
                }
            }

            if (bestI < 0) break;
            if (opt_.requireNegativeSurrogate && !(bestDelta < -opt_.eps)) break;

            applySwap_(cand.node, pos, bestI, bestJ, bestK, bestRev1, bestRev2);
        }

        // DP verify: once per attempt (as in LK-lite)
        double candCost = dp_->reassignCars(cand.node, cand.car);
        if (!std::isfinite(candCost)) continue;

        // Only polish if already better
        if (candCost + opt_.eps < bestCost)
        {
            if (opt_.polish) {
                opt_.polish->improve(cand, candCost);
                // polish is an LS chain that may modify the tour
                if (dp_) {
                    candCost = dp_->reassignCars(cand.node, cand.car);
                    if (!std::isfinite(candCost)) continue;
                }
            }

            if (candCost + opt_.eps < bestCost)
            {
                bestCost = candCost;
                bestSol = std::move(cand);
            }
        }
    }

    if (bestCost + opt_.eps < startCost)
    {
        s = std::move(bestSol);
        cost = bestCost;
    }
}

} // namespace localSearch
