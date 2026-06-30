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

void LinKernighanLocalSearch::resetSeed(uint64_t seed) const
{
    rng_.reseed(seed);
}

bool LinKernighanLocalSearch::isValid2OptMove(int i, int k, int N)
{
    if (i < 1) return false;        // ne diraj node[0]==0
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

    // U LK algoritmu, nakon što uklonimo brid, tražimo novi brid koji će povezati
    // posljednji čvor s nekim drugim čvorom i obrnuti segment između njih
    
    if (lastRemoved < 0 || lastRemoved >= N) 
    {
        // Početak: tražimo sve moguće 2-opt poteze
        for (int i = 1; i < N - 1; ++i)
        {
            const int fromNode = state.tour[(std::size_t)(i - 1)];
            const auto& candidates = opt_.cand->candidates(fromNode);
            
            for (int candidate : candidates)
            {
                const int k = state.pos[(std::size_t)candidate];
                if (k <= i || k >= N) continue;
                
                if (!isValid2OptMove(i, k, N)) continue;
                
                // Provjeri je li već korišten
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
        // Nastavak: tražimo potez koji povezuje posljednji čvor
        const int fromNode = state.tour[(std::size_t)lastRemoved];
        const auto& candidates = opt_.cand->candidates(fromNode);
        
        for (int candidate : candidates)
        {
            const int k = state.pos[(std::size_t)candidate];
            if (k < 0 || k >= N) continue;
            
            // Pronađi i tako da reverse [i, k] povezuje fromNode s candidate
            // Tražimo i gdje je tour[i-1] == fromNode ili tour[i] == fromNode
            for (int i = 1; i < N; ++i)
            {
                if (!isValid2OptMove(i, k, N)) continue;
                
                // Provjeri je li već korišten
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
                
                // Provjeri je li ovaj potez povezuje fromNode s candidate
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
    
    // Ako je gain pozitivan, možemo prihvatiti sekvencu
    if (state.totalGain + bestMove.gain < -opt_.eps)
    {
        // Primijeni potez
        applyReverseAndUpdatePos(state.tour, state.pos, bestMove.i, bestMove.k);
        state.totalGain += bestMove.gain;
        state.sequence.push_back(bestMove);
        
        // Provjeri je li ovo bolje rješenje
        if (state.totalGain < -opt_.eps)
        {
            return true;  // Prihvati sekvencu
        }
        
        // Nastavi tražiti dalje
        const int newLastRemoved = bestMove.k;
        if (searchLKSequence(state, depth + 1, newLastRemoved))
        {
            return true;
        }
        
        // Backtrack: vrati potez
        applyReverseAndUpdatePos(state.tour, state.pos, bestMove.i, bestMove.k);
        state.totalGain -= bestMove.gain;
        state.sequence.pop_back();
    }
    
    // Pokušaj s drugim potezom (ako postoji)
    // U punom LK algoritmu, testiramo sve moguće poteze
    // Ovdje koristimo greedy pristup s backtracking-om
    
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

        // Inicijaliziraj LK state
        LKState state;
        state.tour = cand.node;
        state.pos.assign((std::size_t)N, -1);
        for (int i = 0; i < N; ++i)
            state.pos[(std::size_t)state.tour[(std::size_t)i]] = i;
        state.totalGain = 0.0;
        state.used.assign((std::size_t)N, false);

        // Pokušaj pronaći LK sekvencu
        // Počni s random početnom pozicijom
        const int startPos = rng_.uniformInt(1, N - 2);
        bool found = searchLKSequence(state, 0, startPos);

        if (found && !state.sequence.empty())
        {
            // Primijeni najbolju sekvencu
            cand.node = state.tour;
            
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
    }

    if (bestCost + opt_.eps < startCost) {
        s = std::move(bestSol);
        cost = bestCost;
    }
}

} // namespace localSearch
