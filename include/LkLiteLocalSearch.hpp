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

struct LkLiteOptions
{
    // Number of attempts (each attempt = one dp->reassignCars call).
    int attempts = 6;

    // Length of the surrogate-guided 2-opt move chain.
    int chainLen = 12;

    // Random i-positions tried per step (each looks at topKPerI candidates).
    int triesPerStep = 40;

    // Candidates taken from the candidate list per anchor for a given i.
    int topKPerI = 8;

    // If true, only apply a move when surrogate delta < -eps.
    bool requireNegativeSurrogate = true;

    double eps = 1e-12;

    // Candidate list (required).
    std::shared_ptr<const aco::CandidateListCache> cand = nullptr;

    // Optional polish LS (e.g. 2opt+reloc+2opt) run after a successful LK-lite verify.
    std::shared_ptr<const aco::ILocalSearch> polish = nullptr;
};

// LK-lite: runs a chain of 2-opt moves guided by a cheap surrogate (minTravel),
// then verifies once via the DP (reassignCars). Triggered on stagnation.
class LkLiteLocalSearch final : public aco::ILocalSearch
{
public:
    LkLiteLocalSearch(std::shared_ptr<const cars_tsplib::Instance> inst,
                      LkLiteOptions opt,
                      std::shared_ptr<FixedTourCarAssignerDP> dp);

    // Resets RNG state; call at the start of a stagnation event with a derived seed.
    void resetSeed(uint64_t seed) const override;

    void improve(aco::Solution& s, double& cost) const override;

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    LkLiteOptions opt_;
    std::shared_ptr<FixedTourCarAssignerDP> dp_;

    mutable aco::Rng rng_;

    static bool isValid2OptMove(int i, int k, int N);
    static void applyReverseAndUpdatePos(std::vector<int>& tour, std::vector<int>& pos, int i, int k);

    double minTravel_(int i, int j) const;
    double surrogate2OptDelta_(const std::vector<int>& tour, int i, int k) const;
};

} // namespace localSearch
