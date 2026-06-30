#include <iostream>
#include "LkLiteLocalSearch.hpp"
#include "parser.hpp"

namespace localSearch
{

LkLiteLocalSearch::LkLiteLocalSearch(std::shared_ptr<const cars_tsplib::Instance> inst,
                                     LkLiteOptions opt,
                                     std::shared_ptr<FixedTourCarAssignerDP> dp)
    : inst_(std::move(inst)), opt_(std::move(opt)), dp_(std::move(dp)), rng_(0)
{
    if (!inst_) throw std::runtime_error("LkLiteLocalSearch: inst is null");
    if (!dp_)   throw std::runtime_error("LkLiteLocalSearch: dp is null");
    if (!opt_.cand) throw std::runtime_error("LkLiteLocalSearch: cand is null");
}

void LkLiteLocalSearch::resetSeed(uint64_t seed) const
{
    rng_.reseed(seed);
}

bool LkLiteLocalSearch::isValid2OptMove(int i, int k, int N)
{
    if (i < 1) return false;        // ne diraj node<a href="" class="citation-link" target="_blank" style="vertical-align: super; font-size: 0.8em; margin-left: 3px;">[0]</a>==0
    if (k <= i) return false;
    if (k >= N) return false;
    if (k == i + 1) return false;   // no-op
    if (i == 0 && k == N - 1) return false;
    return true;
}

void LkLiteLocalSearch::applyReverseAndUpdatePos(std::vector<int>& tour,
                                                std::vector<int>& pos,
                                                int i, int k)
{
    std::reverse(tour.begin() + i, tour.begin() + (k + 1));
    for (int idx = i; idx <= k; ++idx) {
        pos[(size_t)tour[(size_t)idx]] = idx;
    }
}

double LkLiteLocalSearch::minTravel_(int i, int j) const
{
    const int C = inst_->cars();
    double best = std::numeric_limits<double>::infinity();
    for (int c = 0; c < C; ++c) {
        best = std::min(best, inst_->travelCost(c, i, j));
    }
    return best;
}

double LkLiteLocalSearch::surrogate2OptDelta_(const std::vector<int>& tour, int i, int k) const
{
    const int N = (int)tour.size();
    const int a = tour[(size_t)(i - 1)];
    const int b = tour[(size_t)i];
    const int c = tour[(size_t)k];
    const int d = tour[(size_t)((k + 1) % N)];

    const double rem = minTravel_(a, b) + minTravel_(c, d);
    const double add = minTravel_(a, c) + minTravel_(b, d);
    return add - rem;
}

void LkLiteLocalSearch::improve(aco::Solution& s, double& cost) const
{
    const int N = inst_->n();
    if (N < 4) return;
    if ((int)s.node.size() != N) throw std::runtime_error("LkLiteLocalSearch: node size mismatch");
    if (s.node[0] != 0) throw std::runtime_error("LkLiteLocalSearch: node[0] must be 0");
    if ((int)s.car.size() != N) s.car.assign((size_t)N, 0);
    if (!std::isfinite(cost)) return;

    const double startCost = cost;

    aco::Solution bestSol = s;
    double bestCost = cost;

    for (int att = 0; att < opt_.attempts; ++att)
    {
        aco::Solution cand = s;

        std::vector<int> pos((size_t)N, -1);
        for (int i = 0; i < N; ++i) pos[(size_t)cand.node[(size_t)i]] = i;

        // chain surrogate 2-opt
        for (int step = 0; step < opt_.chainLen; ++step)
        {
            int bestI = -1, bestK = -1;
            double bestDelta = std::numeric_limits<double>::infinity();

            for (int t = 0; t < opt_.triesPerStep; ++t)
            {
                const int i = rng_.uniformInt(1, N - 2);

                const int anchors[2] = {
                    cand.node[(size_t)(i - 1)],
                    cand.node[(size_t)i]
                };

                for (int aidx = 0; aidx < 2; ++aidx)
                {
                    const auto& neigh = opt_.cand->candidates(anchors[aidx]);
                    int taken = 0;

                    for (int v : neigh)
                    {
                        const int k = pos[(size_t)v];
                        if (!isValid2OptMove(i, k, N)) continue;

                        const double delta = surrogate2OptDelta_(cand.node, i, k);
                        if (delta < bestDelta) {
                            bestDelta = delta;
                            bestI = i;
                            bestK = k;
                        }

                        if (++taken >= opt_.topKPerI) break;
                    }
                }
            }

            if (bestI < 0) break;
            if (opt_.requireNegativeSurrogate && !(bestDelta < -opt_.eps)) break;

            applyReverseAndUpdatePos(cand.node, pos, bestI, bestK);
        }

        // DP verify: 1x po attemptu
        double candCost = dp_->reassignCars(cand.node, cand.car);
        if (!std::isfinite(candCost)) continue;

        // polish samo ako je već bolje
        if (candCost + opt_.eps < bestCost)
        {
            if (opt_.polish) {
                opt_.polish->improve(cand, candCost);
                // (polish je LS chain koji može mijenjati turu)
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

    if (bestCost + opt_.eps < startCost) {
        s = std::move(bestSol);
        cost = bestCost;
    }
}

} // namespace localSearch
