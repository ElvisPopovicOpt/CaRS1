#include "ThreeOptLocalSearch.hpp"

#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "parser.hpp"
#include "candidate_list.hpp"
#include "FixedTourCarAssignerDP.hpp"
#include "rng.hpp"

namespace localSearch
{

ThreeOptLocalSearch::ThreeOptLocalSearch(std::shared_ptr<const cars_tsplib::Instance> inst,
                                         ThreeOptOptions opt,
                                         std::shared_ptr<FixedTourCarAssignerDP> dp)
    : inst_(std::move(inst)), opt_(std::move(opt)), dp_(std::move(dp)), rng_(0)
{
    if (!inst_) throw std::runtime_error("ThreeOptLocalSearch: inst is null");
    if (!dp_)   throw std::runtime_error("ThreeOptLocalSearch: dp is null");
    if (!opt_.cand) throw std::runtime_error("ThreeOptLocalSearch: cand is null");

    const int N = inst_->n();
    const int C = inst_->cars();
    if (N <= 0) throw std::runtime_error("ThreeOptLocalSearch: inst->n() <= 0");
    if (C <= 0) throw std::runtime_error("ThreeOptLocalSearch: inst->cars() <= 0");

    // Precompute minTravel
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

void ThreeOptLocalSearch::resetSeed(uint64_t seed) const
{
    rng_.reseed(seed);
}

double ThreeOptLocalSearch::minTravelFast_(int u, int v) const
{
    const int N = inst_->n();
    return minTravel_[(std::size_t)u * (std::size_t)N + (std::size_t)v];
}

bool ThreeOptLocalSearch::isValidCut_(int i, int j, int k, int N)
{
    if (N < 6) return false;
    if (i < 1) return false;
    if (!(i < j && j < k)) return false;
    if (k >= N) return false;
    if (j == i) return false;
    if (k == j) return false;
    return true;
}

double ThreeOptLocalSearch::surrogate3OptDelta_(const std::vector<int>& tour,
                                                 int i, int j, int k,
                                                 int variant) const
{
    const int N = (int)tour.size();
    const int a = tour[(std::size_t)(i - 1)];
    const int b = tour[(std::size_t)i];
    const int c = tour[(std::size_t)(j - 1)];
    const int d = tour[(std::size_t)j];
    const int e = tour[(std::size_t)(k - 1)];
    const int f = tour[(std::size_t)k];

    const double rem = minTravelFast_(a, b) + minTravelFast_(c, d) + minTravelFast_(e, f);

    double add = 0.0;
    
    // 7 varijanti 3-opt move-a
    switch (variant)
    {
        case 0: // Original (no change) - ne koristimo
            add = rem;
            break;
        case 1: // Reverse [i, j)
            add = minTravelFast_(a, c) + minTravelFast_(b, d) + minTravelFast_(e, f);
            break;
        case 2: // Reverse [j, k)
            add = minTravelFast_(a, b) + minTravelFast_(c, e) + minTravelFast_(d, f);
            break;
        case 3: // Reverse [i, k)
            add = minTravelFast_(a, e) + minTravelFast_(b, f) + minTravelFast_(c, d);
            break;
        case 4: // Swap segments: S0 + S2 + S1 + S3 (no reverse)
            add = minTravelFast_(a, d) + minTravelFast_(c, b) + minTravelFast_(e, f);
            break;
        case 5: // Swap segments: S0 + rev(S2) + S1 + S3
            add = minTravelFast_(a, e) + minTravelFast_(d, b) + minTravelFast_(c, f);
            break;
        case 6: // Swap segments: S0 + S2 + rev(S1) + S3
            add = minTravelFast_(a, d) + minTravelFast_(e, c) + minTravelFast_(b, f);
            break;
        default:
            add = rem;
            break;
    }

    return add - rem;
}

void ThreeOptLocalSearch::apply3OptMove_(std::vector<int>& tour,
                                         int i, int j, int k,
                                         int variant) const
{
    const int N = (int)tour.size();
    std::vector<int> tmp;
    tmp.reserve((std::size_t)N);

    switch (variant)
    {
        case 1: // Reverse [i, j)
            tmp = tour;
            std::reverse(tmp.begin() + i, tmp.begin() + j);
            tour.swap(tmp);
            break;
        case 2: // Reverse [j, k)
            tmp = tour;
            std::reverse(tmp.begin() + j, tmp.begin() + k);
            tour.swap(tmp);
            break;
        case 3: // Reverse [i, k)
            tmp = tour;
            std::reverse(tmp.begin() + i, tmp.begin() + k);
            tour.swap(tmp);
            break;
        case 4: // Swap: S0 + S2 + S1 + S3
            tmp.clear();
            // S0: [0, i)
            for (int a = 0; a < i; ++a) tmp.push_back(tour[(std::size_t)a]);
            // S2: [j, k)
            for (int a = j; a < k; ++a) tmp.push_back(tour[(std::size_t)a]);
            // S1: [i, j)
            for (int a = i; a < j; ++a) tmp.push_back(tour[(std::size_t)a]);
            // S3: [k, N)
            for (int a = k; a < N; ++a) tmp.push_back(tour[(std::size_t)a]);
            tour.swap(tmp);
            break;
        case 5: // Swap: S0 + rev(S2) + S1 + S3
            tmp.clear();
            for (int a = 0; a < i; ++a) tmp.push_back(tour[(std::size_t)a]);
            for (int a = k - 1; a >= j; --a) tmp.push_back(tour[(std::size_t)a]);
            for (int a = i; a < j; ++a) tmp.push_back(tour[(std::size_t)a]);
            for (int a = k; a < N; ++a) tmp.push_back(tour[(std::size_t)a]);
            tour.swap(tmp);
            break;
        case 6: // Swap: S0 + S2 + rev(S1) + S3
            tmp.clear();
            for (int a = 0; a < i; ++a) tmp.push_back(tour[(std::size_t)a]);
            for (int a = j; a < k; ++a) tmp.push_back(tour[(std::size_t)a]);
            for (int a = j - 1; a >= i; --a) tmp.push_back(tour[(std::size_t)a]);
            for (int a = k; a < N; ++a) tmp.push_back(tour[(std::size_t)a]);
            tour.swap(tmp);
            break;
        default:
            break;
    }
}

void ThreeOptLocalSearch::improve(aco::Solution& s, double& cost) const
{
    const int N = inst_->n();
    if (N < 6) return;
    if ((int)s.node.size() != N) throw std::runtime_error("ThreeOptLocalSearch: node size mismatch");
    if (s.node[0] != 0) throw std::runtime_error("ThreeOptLocalSearch: node[0] must be 0");
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

        bool improved = true;
        for (int pass = 0; pass < opt_.maxPasses && improved; ++pass)
        {
            improved = false;
            double bestDelta = 0.0;
            int bestI = -1, bestJ = -1, bestK = -1, bestVariant = -1;

            // Testiraj sve moguće 3-opt poteze
            for (int i = 1; i < N - 4; ++i)
            {
                const int anchor1 = cand.node[(std::size_t)(i - 1)];
                const auto& neighJ = opt_.cand->candidates(anchor1);
                
                for (int vj : neighJ)
                {
                    const int j = pos[(std::size_t)vj];
                    if (j <= i || j >= N - 2) continue;

                    const int anchor2 = cand.node[(std::size_t)(j - 1)];
                    const auto& neighK = opt_.cand->candidates(anchor2);
                    
                    for (int vk : neighK)
                    {
                        const int k = pos[(std::size_t)vk];
                        if (k <= j || k >= N) continue;

                        if (!isValidCut_(i, j, k, N)) continue;

                        // Testiraj sve 7 varijanti (osim 0 - original)
                        for (int variant = 1; variant <= 6; ++variant)
                        {
                            const double delta = surrogate3OptDelta_(cand.node, i, j, k, variant);
                            if (delta < bestDelta)
                            {
                                bestDelta = delta;
                                bestI = i;
                                bestJ = j;
                                bestK = k;
                                bestVariant = variant;
                            }
                        }
                    }
                }
            }

            if (bestI >= 0 && bestDelta < -opt_.eps)
            {
                apply3OptMove_(cand.node, bestI, bestJ, bestK, bestVariant);
                
                // Rebuild pos
                std::fill(pos.begin(), pos.end(), -1);
                for (int idx = 0; idx < N; ++idx)
                    pos[(std::size_t)cand.node[(std::size_t)idx]] = idx;
                
                improved = true;
            }
        }

        // DP verify
        double candCost = dp_->reassignCars(cand.node, cand.car);
        if (!std::isfinite(candCost)) continue;

        // Polish ako je postavljen
        if (candCost + opt_.eps < bestCost)
        {
            if (opt_.polish) {
                opt_.polish->improve(cand, candCost);
                if (dp_) {
                    candCost = dp_->reassignCars(cand.node, cand.car);
                    if (!std::isfinite(candCost)) continue;
                }
            }

            if (candCost + opt_.eps < bestCost) {
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
