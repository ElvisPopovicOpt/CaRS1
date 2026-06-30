#pragma once
#include <memory>
#include <vector>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <algorithm>

#include "interfaces.hpp"          // aco::ILocalSearch, aco::Solution
#include "rng.hpp"                 // aco::Rng
#include "candidate_list.hpp"  // CandidateListCache
#include "FixedTourCarAssignerDP.hpp"

namespace localSearch
{

struct LinKernighanOptions
{
    // Number of attempts (each attempt = one full LK search)
    int attempts = 3;

    // Maximum backtracking depth (max sequence length)
    int maxDepth = 20;

    // Maximum number of sequences tested per attempt
    int maxSequences = 1000;

    // Candidate list (required)
    std::shared_ptr<const aco::CandidateListCache> cand = nullptr;

    // Optional: after a successful LK verify, run a polish LS (e.g. 2opt+reloc+2opt)
    std::shared_ptr<const aco::ILocalSearch> polish = nullptr;

    double eps = 1e-12;
};

// Full Lin-Kernighan algorithm with backtracking.
// Builds sequences of 2-opt moves and accepts them if total gain > 0.
class LinKernighanLocalSearch final : public aco::ILocalSearch
{
public:
    LinKernighanLocalSearch(std::shared_ptr<const cars_tsplib::Instance> inst,
                           LinKernighanOptions opt,
                           std::shared_ptr<FixedTourCarAssignerDP> dp);

    // Reset RNG state (call at the start of a stagnation event with a derived seed)
    void resetSeed(uint64_t seed) const override;

    void improve(aco::Solution& s, double& cost) const override;

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    LinKernighanOptions opt_;
    std::shared_ptr<FixedTourCarAssignerDP> dp_;

    mutable aco::Rng rng_;

    // Precomputed min travel cost
    std::vector<double> minTravel_;

    // 2-opt move structure
    struct Move
    {
        int i, k;  // reverse segment [i, k]
        double gain;  // gain from this move
    };

    // LK state for backtracking
    struct LKState
    {
        std::vector<int> tour;
        std::vector<int> pos;  // position map
        double totalGain;
        std::vector<Move> sequence;
        std::vector<bool> used;  // used[i] = true if edge i was removed
    };

    static bool isValid2OptMove(int i, int k, int N);
    static void applyReverseAndUpdatePos(std::vector<int>& tour, std::vector<int>& pos, int i, int k);

    double minTravelFast_(int u, int v) const;
    double surrogate2OptDelta_(const std::vector<int>& tour, int i, int k) const;

    // Full LK search with backtracking
    bool searchLKSequence(LKState& state, int depth, int lastRemoved) const;

    // Find the best next move
    bool findBestNextMove(const LKState& state, int lastRemoved, Move& bestMove) const;
};

} // namespace localSearch
