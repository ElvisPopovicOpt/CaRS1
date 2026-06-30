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
    // Koliko pokušaja (svaki pokušaj = 1x puni LK search)
    int attempts = 3;

    // Maksimalna dubina backtracking-a (maksimalna duljina sekvence)
    int maxDepth = 20;

    // Maksimalan broj sekvenci za testiranje po pokušaju
    int maxSequences = 1000;

    // Candidate list (obavezno)
    std::shared_ptr<const aco::CandidateListCache> cand = nullptr;

    // Opcionalno: nakon uspješnog LK verify, pusti polish LS (npr. 2opt+reloc+2opt)
    std::shared_ptr<const aco::ILocalSearch> polish = nullptr;

    double eps = 1e-12;
};

// Puni Lin-Kernighan algoritam s backtracking-om
// Gradi sekvence 2-opt poteza i prihvaća ih ako je ukupni gain > 0
class LinKernighanLocalSearch final : public aco::ILocalSearch
{
public:
    LinKernighanLocalSearch(std::shared_ptr<const cars_tsplib::Instance> inst,
                           LinKernighanOptions opt,
                           std::shared_ptr<FixedTourCarAssignerDP> dp);

    // Reset RNG stanja (pozovi na početku stagnation-eventa s deriviranim seedom)
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
    
    // Puni LK search s backtracking-om
    bool searchLKSequence(LKState& state, int depth, int lastRemoved) const;
    
    // Pronađi najbolji sljedeći potez
    bool findBestNextMove(const LKState& state, int lastRemoved, Move& bestMove) const;
};

} // namespace localSearch
