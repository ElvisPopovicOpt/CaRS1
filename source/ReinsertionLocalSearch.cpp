#include <algorithm>
#include <limits>
#include <stdexcept>
#include "parser.hpp"
#include "ReinsertionLocalSearch.hpp"
#include "FixedTourCarAssignerDP.hpp"
#include "candidate_list.hpp"

namespace localSearch
{

static void normalizeTour(std::vector<int>& nodes)
{
    if (nodes.empty()) return;
    const int N = static_cast<int>(nodes.size());
    int pos0 = -1;
    for (int idx = 0; idx < N; ++idx) {
        if (nodes[static_cast<size_t>(idx)] == 0) {
            pos0 = idx;
            break;
        }
    }
    if (pos0 > 0 && pos0 < N) {
        std::rotate(nodes.begin(), nodes.begin() + pos0, nodes.end());
    }
}

ReinsertionLocalSearch::ReinsertionLocalSearch(
    std::shared_ptr<const cars_tsplib::Instance> inst,
    std::shared_ptr<FixedTourCarAssignerDP> dp,
    std::shared_ptr<const aco::CandidateListCache> cand,
    ReinsertionOptions opt)
    : inst_(std::move(inst))
    , dp_(std::move(dp))
    , cand_(std::move(cand))
    , opt_(opt)
{
    if (!inst_) throw std::runtime_error("ReinsertionLocalSearch: inst is null");
    if (!dp_)  throw std::runtime_error("ReinsertionLocalSearch: dp is null");
}

double ReinsertionLocalSearch::evaluateReinsertionMove(
    const std::vector<int>& nodes, int nodeIdx, int newPos) const
{
    const int N = static_cast<int>(nodes.size());
    if (nodeIdx < 1 || nodeIdx >= N || newPos < 1 || newPos >= N)
        return std::numeric_limits<double>::infinity();
    if (nodeIdx == newPos || nodeIdx == newPos - 1)
        return std::numeric_limits<double>::infinity();

    const int nodeToMove = nodes[static_cast<size_t>(nodeIdx)];
    int insertPos = newPos;
    if (newPos > nodeIdx)
        insertPos = newPos - 1;

    std::vector<int> newNodes;
    newNodes.reserve(static_cast<size_t>(N));
    for (int i = 0; i < N; ++i) {
        if (i == nodeIdx) continue;
        if (static_cast<int>(newNodes.size()) == insertPos)
            newNodes.push_back(nodeToMove);
        newNodes.push_back(nodes[static_cast<size_t>(i)]);
    }
    if (static_cast<int>(newNodes.size()) == insertPos)
        newNodes.push_back(nodeToMove);

    if (static_cast<int>(newNodes.size()) != N || newNodes[0] != 0)
        return std::numeric_limits<double>::infinity();

    return dp_->evaluateCost(newNodes);
}

void ReinsertionLocalSearch::applyReinsertionMove(
    std::vector<int>& nodes, int nodeIdx, int newPos) const
{
    const int N = static_cast<int>(nodes.size());
    if (nodeIdx < 1 || nodeIdx >= N || newPos < 1 || newPos >= N) return;
    if (nodeIdx == newPos || nodeIdx == newPos - 1) return;

    const int nodeToMove = nodes[static_cast<size_t>(nodeIdx)];
    nodes.erase(nodes.begin() + nodeIdx);

    int insertPos = newPos;
    if (newPos > nodeIdx)
        insertPos = newPos - 1;
    nodes.insert(nodes.begin() + insertPos, nodeToMove);

    normalizeTour(nodes);
}

void ReinsertionLocalSearch::improve(aco::Solution& s, double& cost) const
{
    const int N = inst_->n();
    if (N < 4) return;
    if ((int)s.node.size() != N)
        throw std::runtime_error("ReinsertionLocalSearch: node size mismatch");
    if (s.node[0] != 0)
        throw std::runtime_error("ReinsertionLocalSearch: node[0] must be 0");
    if ((int)s.car.size() != N)
        s.car.assign(static_cast<size_t>(N), 0);

    const double eps = 1e-9;
    int passes = 0;
    int consecutiveNoImprovement = 0;
    const int maxPasses = opt_.maxPasses;
    const bool useCand = (cand_ && cand_->N() == N && cand_->K() > 0);

    while (passes < maxPasses) {
        ++passes;
        bool improved = false;

        for (int nodeIdx = 1; nodeIdx < N; ++nodeIdx) {
            const int nodeValue = s.node[static_cast<size_t>(nodeIdx)];

            std::vector<int> targetPositions;
            if (useCand) {
                for (int cand : cand_->candidates(nodeValue)) {
                    for (int pos = 1; pos < N; ++pos) {
                        if (s.node[static_cast<size_t>(pos)] == cand) {
                            if (pos > 0) targetPositions.push_back(pos);
                            if (pos + 1 < N) targetPositions.push_back(pos + 1);
                            break;
                        }
                    }
                }
            } else {
                for (int pos = 1; pos < N; ++pos)
                    targetPositions.push_back(pos);
            }

            std::sort(targetPositions.begin(), targetPositions.end());
            targetPositions.erase(
                std::unique(targetPositions.begin(), targetPositions.end()),
                targetPositions.end());

            for (int newPos : targetPositions) {
                if (newPos < 1 || newPos >= N) continue;
                if (newPos == nodeIdx || newPos == nodeIdx + 1) continue;

                const double newCost = evaluateReinsertionMove(s.node, nodeIdx, newPos);
                if (!std::isfinite(newCost) || newCost >= cost - eps)
                    continue;

                std::vector<int> savedNodes = s.node;
                std::vector<int> savedCars = s.car;

                applyReinsertionMove(s.node, nodeIdx, newPos);
                if (s.node.empty() || s.node[0] != 0) {
                    s.node = savedNodes;
                    s.car = savedCars;
                    continue;
                }

                if ((int)s.car.size() != N)
                    s.car.assign(static_cast<size_t>(N), 0);
                const double reassignedCost = dp_->reassignCars(s.node, s.car);

                if (!std::isfinite(reassignedCost) || reassignedCost >= cost - eps) {
                    s.node = savedNodes;
                    s.car = savedCars;
                    continue;
                }

                cost = reassignedCost;
                improved = true;
                consecutiveNoImprovement = 0;
                break;
            }
            if (improved) break;
        }

        if (!improved) {
            consecutiveNoImprovement++;
            if (consecutiveNoImprovement >= 2)
                break;
        }
    }
}

} // namespace localSearch
