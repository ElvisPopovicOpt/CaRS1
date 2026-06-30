#pragma once

#include <memory>
#include <vector>
#include <limits>
#include <cstdint>

#include "interfaces.hpp"

namespace cars_tsplib { struct Instance; }
namespace aco { class CandidateListCache; class Rng; }

namespace localSearch { class FixedTourCarAssignerDP; }

namespace localSearch
{

struct ThreeOptLiteOptions
{
    std::shared_ptr<const aco::CandidateListCache> cand;

    int attempts = 6;
    int chainLen = 8;
    int triesPerStep = 30;

    // Max candidates taken from the candidate list at each level (for j and k).
    int topKPerI = 8;

    double eps = 1e-12;
    bool requireNegativeSurrogate = true;

    // If set, invoked only after the DP-verified solution is already better.
    std::shared_ptr<const aco::ILocalSearch> polish;
};

class ThreeOptLiteLocalSearch final : public aco::ILocalSearch
{
public:
    ThreeOptLiteLocalSearch(std::shared_ptr<const cars_tsplib::Instance> inst,
                            ThreeOptLiteOptions opt,
                            std::shared_ptr<FixedTourCarAssignerDP> dp);

    void resetSeed(uint64_t seed) const override;

    void improve(aco::Solution& s, double& cost) const override;

private:
    // Precomputed min_c travelCost(c,u,v)
    double minTravelFast_(int u, int v) const;

    // Surrogate delta for a 3-opt segment swap (S1 <-> S2), with optional segment reversal.
    double surrogate3OptDeltaSwap_(const std::vector<int>& tour,
                                   int i, int j, int k,
                                   bool revS1, bool revS2) const;

    // Applies the swap move and rebuilds pos (simple and fast enough for N<=300).
    void applySwap_(std::vector<int>& tour,
                    std::vector<int>& pos,
                    int i, int j, int k,
                    bool revS1, bool revS2) const;

    static bool isValidCut_(int i, int j, int k, int N);

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    ThreeOptLiteOptions opt_;
    std::shared_ptr<localSearch::FixedTourCarAssignerDP> dp_;

    // minTravel[u*N + v]: min over cars c of travelCost(c,u,v)
    std::vector<double> minTravel_;

    mutable aco::Rng rng_;
};

} // namespace localSearch
