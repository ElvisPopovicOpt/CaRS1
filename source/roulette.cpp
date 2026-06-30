#include <cmath>
#include <limits>
#include <stdexcept>
#include <parser.hpp>
#include <rng.hpp>
#include <interfaces.hpp>
#include <params_cli.hpp>
#include <candidate_list.hpp>
#include <roulette.hpp>

namespace aco 
{

AntPolicyCandidateListRoulette::AntPolicyCandidateListRoulette(
    std::shared_ptr<const CandidateListCache> cl,
    double epsilonEta)
    : cl_(std::move(cl)), eps_(epsilonEta) {}

static inline double eta(double eps, double cost) 
{
    return 1.0 / (eps + cost);
}

Solution AntPolicyCandidateListRoulette::construct(const AntContext& ctx) const {
    if (!ctx.inst || !ctx.params || !ctx.pher || !ctx.rng) {
        throw std::runtime_error("AntPolicyCandidateListRoulette: ctx.inst/params/pher/rng must be set.");
    }
    auto& rng = *ctx.rng;

    const auto& inst = *ctx.inst;
    const int N = inst.n();
    const int C = inst.cars();

    Solution s;
    s.node.resize(static_cast<size_t>(N));
    s.car.resize(static_cast<size_t>(N));

    std::vector<uint8_t> visited(static_cast<size_t>(N), 0);
    s.node[0] = 0;
    visited[0] = 1;

    const double alphaNodes = ctx.params->alphaNodes;
    const double betaNodes  = ctx.params->betaNodes;
    const double alphaCars  = ctx.params->alphaCars;
    const double betaCars   = ctx.params->betaCars;

    auto choose_from_pairs = [&](int i, const std::vector<int>& nextCandidates, int forcedNextOrNeg1) -> std::pair<int,int> 
    {
        struct Item { int j; int c; double w; };
        std::vector<Item> items;
        items.reserve(static_cast<size_t>(nextCandidates.size() * std::max(1, C)));

        double sum = 0.0;

        for (int j : nextCandidates) {
            if (forcedNextOrNeg1 >= 0 && j != forcedNextOrNeg1) continue;

            // Weight for node j (tour movement): alphaNodes, betaNodes
            double tauNode = 0.0;
            double minCost = std::numeric_limits<double>::max();
            for (int c = 0; c < C; ++c) {
                tauNode += ctx.pher->tauMove(c, i, j);
                minCost = std::min(minCost, inst.travelCost(c, i, j));
            }
            const double etaNode = 1.0 / (eps_ + minCost);
            const double W_node = std::pow(tauNode, alphaNodes) * std::pow(etaNode, betaNodes);
            if (!(W_node > 0.0) || !std::isfinite(W_node)) continue;

            for (int c = 0; c < C; ++c) {
                const double t = ctx.pher->tauMove(c, i, j);
                const double h = 1.0 / (eps_ + inst.travelCost(c, i, j));
                double w = W_node * (std::pow(t, alphaCars) * std::pow(h, betaCars));

                if (!(w > 0.0) || !std::isfinite(w)) continue;

                items.push_back({j, c, w});
                sum += w;
            }
        }

        // No valid weights: fall back to a uniform, stable choice.
        if (items.empty() || !(sum > 0.0) || !std::isfinite(sum)) {
            int j = forcedNextOrNeg1 >= 0 ? forcedNextOrNeg1
                                          : nextCandidates[rng.uniformInt((int)nextCandidates.size())];
            int c = rng.uniformInt(C);
            return {j, c};
        }

        double r = rng.uniform01() * sum;
        double acc = 0.0;
        for (const auto& it : items) {
            acc += it.w;
            if (acc >= r) return {it.j, it.c};
        }
        return {items.back().j, items.back().c};
    };

    for (int pos = 0; pos < N - 1; ++pos) {
        const int i = s.node[(size_t)pos];

        std::vector<int> cand;
        cand.reserve((size_t)std::max(0, ctx.params->favorites));

        if (cl_ && cl_->N() == N) {
            for (int j : cl_->candidates(i)) {
                if (!visited[(size_t)j]) cand.push_back(j);
                if ((int)cand.size() >= ctx.params->favorites) break;
            }
        }

        if (cand.empty()) {
            for (int j = 0; j < N; ++j)
                if (!visited[(size_t)j]) cand.push_back(j);
        }

        auto [nextNode, chosenCar] = choose_from_pairs(i, cand, -1);
        s.car[(size_t)pos] = chosenCar;
        s.node[(size_t)pos + 1] = nextNode;
        visited[(size_t)nextNode] = 1;
    }

    const int lastNode = s.node[(size_t)N - 1];
    {
        std::vector<int> onlyStart = {0};
        auto [_, chosenCar] = choose_from_pairs(lastNode, onlyStart, 0);
        (void)_;
        s.car[(size_t)N - 1] = chosenCar;
    }

    return s;
}

} // namespace aco
