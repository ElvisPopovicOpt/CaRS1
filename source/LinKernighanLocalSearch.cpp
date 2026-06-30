#include "LinKernighanLocalSearch.hpp"
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include "parser.hpp"

namespace localSearch
{

LinKernighanLocalSearch::LinKernighanLocalSearch(std::shared_ptr<const cars_tsplib::Instance> inst,
                                                 LinKernighanOptions opt,
                                                 std::shared_ptr<FixedTourCarAssignerDP> dp)
    : inst_(std::move(inst)), opt_(std::move(opt)), dp_(std::move(dp)), rng_(0)
{
    if (!inst_) throw std::runtime_error("LinKernighanLocalSearch: inst is null");
    if (!dp_)   throw std::runtime_error("LinKernighanLocalSearch: dp is null");
    if (!opt_.cand) throw std::runtime_error("LinKernighanLocalSearch: cand is null");

    const int N = inst_->n();
    const int C = inst_->cars();
    if (N <= 0) throw std::runtime_error("LinKernighanLocalSearch: inst->n() <= 0");
    if (C <= 0) throw std::runtime_error("LinKernighanLocalSearch: inst->cars() <= 0");

    // Precompute min travel cost
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

void LinKernighanLocalSearch::resetSeed(uint64_t seed) const
{
    rng_.reseed(seed);
}

bool LinKernighanLocalSearch::isValid2OptMove(int i, int k, int N)
{
    if (i < 1) return false;        // never move node[0]==0
    if (k <= i) return false;
    if (k >= N) return false;
    if (k == i + 1) return false;   // no-op
    if (i == 0 && k == N - 1) return false;
    return true;
}

void LinKernighanLocalSearch::applyReverseAndUpdatePos(std::vector<int>& tour,
                                                       std::vector<int>& pos,
                                                       int i, int k)
{
    std::reverse(tour.begin() + i, tour.begin() + (k + 1));
    for (int idx = i; idx <= k; ++idx) {
        pos[(std::size_t)tour[(std::size_t)idx]] = idx;
    }
}

double LinKernighanLocalSearch::minTravelFast_(int u, int v) const
{
    const int N = inst_->n();
    return minTravel_[(std::size_t)u * (std::size_t)N + (std::size_t)v];
}

double LinKernighanLocalSearch::surrogate2OptDelta_(const std::vector<int>& tour, int i, int k) const
{
    const int N = (int)tour.size();
    const int a = tour[(std::size_t)(i - 1)];
    const int b = tour[(std::size_t)i];
    const int c = tour[(std::size_t)k];
    const int d = tour[(std::size_t)((k + 1) % N)];

    const double rem = minTravelFast_(a, b) + minTravelFast_(c, d);
    const double add = minTravelFast_(a, c) + minTravelFast_(b, d);
    return add - rem;
}

bool LinKernighanLocalSearch::findBestNextMove(const LKState& state, int lastRemoved, Move& bestMove) const
{
    const int N = (int)state.tour.size();
    bestMove.gain = std::numeric_limits<double>::infinity();
    bestMove.i = -1;
    bestMove.k = -1;

    // After removing an edge, search for a new edge connecting the last node
    // to some other node, reversing the segment between them.

    if (lastRemoved < 0 || lastRemoved >= N)
    {
        // Start: search all possible 2-opt moves
        for (int i = 1; i < N - 1; ++i)
        {
            const int fromNode = state.tour[(std::size_t)(i - 1)];
            const auto& candidates = opt_.cand->candidates(fromNode);
            
            for (int candidate : candidates)
            {
                const int k = state.pos[(std::size_t)candidate];
                if (k <= i || k >= N) continue;
                
                if (!isValid2OptMove(i, k, N)) continue;
                
                // Skip if already used
                bool used = false;
                for (const auto& m : state.sequence)
                {
                    if ((m.i == i && m.k == k) || (m.i == k && m.k == i))
                    {
                        used = true;
                        break;
                    }
                }
                if (used) continue;

                const double gain = surrogate2OptDelta_(state.tour, i, k);
                if (gain < bestMove.gain)
                {
                    bestMove.gain = gain;
                    bestMove.i = i;
                    bestMove.k = k;
                }
            }
        }
    }
    else
    {
        // Continuation: search for a move that connects the last node
        const int fromNode = state.tour[(std::size_t)lastRemoved];
        const auto& candidates = opt_.cand->candidates(fromNode);
        
        for (int candidate : candidates)
        {
            const int k = state.pos[(std::size_t)candidate];
            if (k < 0 || k >= N) continue;
            
            // Find i such that reversing [i, k] connects fromNode to candidate,
            // i.e. tour[i-1] == fromNode or tour[i] == fromNode
            for (int i = 1; i < N; ++i)
            {
                if (!isValid2OptMove(i, k, N)) continue;

                // Skip if already used
                bool used = false;
                for (const auto& m : state.sequence)
                {
                    if ((m.i == i && m.k == k) || (m.i == k && m.k == i))
                    {
                        used = true;
                        break;
                    }
                }
                if (used) continue;

                // Check whether this move connects fromNode to candidate
                const int prevNode = state.tour[(std::size_t)(i - 1)];
                const int nextNode = state.tour[(std::size_t)((k + 1) % N)];
                
                if (prevNode == fromNode || state.tour[(std::size_t)i] == fromNode)
                {
                    const double gain = surrogate2OptDelta_(state.tour, i, k);
                    if (gain < bestMove.gain)
                    {
                        bestMove.gain = gain;
                        bestMove.i = i;
                        bestMove.k = k;
                    }
                }
            }
        }
    }
    
    return bestMove.i >= 0;
}

bool LinKernighanLocalSearch::searchLKSequence(LKState& state, int depth, int lastRemoved) const
{
    if (depth >= opt_.maxDepth) return false;
    if (state.sequence.size() >= (std::size_t)opt_.maxSequences) return false;
    
    Move bestMove;
    if (!findBestNextMove(state, lastRemoved, bestMove)) return false;
    
    // Only continue down this branch if it can still lead to a net negative total gain
    if (state.totalGain + bestMove.gain < -opt_.eps)
    {
        // Apply the move
        applyReverseAndUpdatePos(state.tour, state.pos, bestMove.i, bestMove.k);
        state.totalGain += bestMove.gain;
        state.sequence.push_back(bestMove);

        if (state.totalGain < -opt_.eps)
        {
            return true;  // Improving sequence found, accept it
        }

        // Keep extending the sequence
        const int newLastRemoved = bestMove.k;
        if (searchLKSequence(state, depth + 1, newLastRemoved))
        {
            return true;
        }

        // Backtrack: undo the move
        applyReverseAndUpdatePos(state.tour, state.pos, bestMove.i, bestMove.k);
        state.totalGain -= bestMove.gain;
        state.sequence.pop_back();
    }

    return false;
}

void LinKernighanLocalSearch::improve(aco::Solution& s, double& cost) const
{
    const int N = inst_->n();
    if (N < 4) return;
    if ((int)s.node.size() != N) throw std::runtime_error("LinKernighanLocalSearch: node size mismatch");
    if (s.node[0] != 0) throw std::runtime_error("LinKernighanLocalSearch: node[0] must be 0");
    if ((int)s.car.size() != N) s.car.assign((std::size_t)N, 0);
    if (!std::isfinite(cost)) return;

    const double startCost = cost;

    aco::Solution bestSol = s;
    double bestCost = cost;

    for (int att = 0; att < opt_.attempts; ++att)
    {
        aco::Solution cand = s;

        // Initialize LK state
        LKState state;
        state.tour = cand.node;
        state.pos.assign((std::size_t)N, -1);
        for (int i = 0; i < N; ++i)
            state.pos[(std::size_t)state.tour[(std::size_t)i]] = i;
        state.totalGain = 0.0;
        state.used.assign((std::size_t)N, false);

        // Try to find an LK sequence, starting from a random position
        const int startPos = rng_.uniformInt(1, N - 2);
        bool found = searchLKSequence(state, 0, startPos);

        if (found && !state.sequence.empty())
        {
            // Apply the best sequence found
            cand.node = state.tour;

            // DP verify
            double candCost = dp_->reassignCars(cand.node, cand.car);
            if (!std::isfinite(candCost)) continue;

            // Polish if configured
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
    }

    if (bestCost + opt_.eps < startCost) {
        s = std::move(bestSol);
        cost = bestCost;
    }
}

} // namespace localSearch
