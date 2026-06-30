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

struct ThreeOptOptions
{
    std::shared_ptr<const aco::CandidateListCache> cand;

    int attempts = 3;
    int maxPasses = 5;  // max number of passes over all candidate 3-opt moves

    double eps = 1e-12;

    // Optional: run this LS as a polish step after a successful 3-opt verify
    std::shared_ptr<const aco::ILocalSearch> polish = nullptr;
};

// Full 3-opt: tests all 7 reconnection variants for each cut
class ThreeOptLocalSearch final : public aco::ILocalSearch
{
public:
    ThreeOptLocalSearch(std::shared_ptr<const cars_tsplib::Instance> inst,
                        ThreeOptOptions opt,
                        std::shared_ptr<FixedTourCarAssignerDP> dp);

    void resetSeed(uint64_t seed) const override;

    void improve(aco::Solution& s, double& cost) const override;

private:
    // Precomputed min over cars of travelCost(c,u,v)
    double minTravelFast_(int u, int v) const;

    // Surrogate cost delta for a 3-opt move (variants 0-6)
    double surrogate3OptDelta_(const std::vector<int>& tour,
                                int i, int j, int k,
                                int variant) const;

    // Apply a 3-opt move (variant 0-6)
    void apply3OptMove_(std::vector<int>& tour,
                        int i, int j, int k,
                        int variant) const;

    static bool isValidCut_(int i, int j, int k, int N);

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    ThreeOptOptions opt_;
    std::shared_ptr<localSearch::FixedTourCarAssignerDP> dp_;

    // Flattened minTravel[u*N + v]
    std::vector<double> minTravel_;

    mutable aco::Rng rng_;
};

} // namespace localSearch
