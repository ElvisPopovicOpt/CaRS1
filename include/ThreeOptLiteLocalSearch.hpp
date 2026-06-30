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

    // Limit koliko kandidata uzimamo iz candidate liste na svakoj razini (za j i za k)
    int topKPerI = 8;

    double eps = 1e-12;
    bool requireNegativeSurrogate = true;

    // Ako je postavljeno, zove se samo kad je DP-verified rješenje već bolje.
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

    // Surrogate delta za 3-opt segment-swap (swap S1 i S2), s opcionalnim reverzijama segmenata.
    double surrogate3OptDeltaSwap_(const std::vector<int>& tour,
                                   int i, int j, int k,
                                   bool revS1, bool revS2) const;

    // Primijeni swap potez i rebuild pos (jednostavno i sigurno za N<=300)
    void applySwap_(std::vector<int>& tour,
                    std::vector<int>& pos,
                    int i, int j, int k,
                    bool revS1, bool revS2) const;

    static bool isValidCut_(int i, int j, int k, int N);

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    ThreeOptLiteOptions opt_;
    std::shared_ptr<localSearch::FixedTourCarAssignerDP> dp_;

    // minTravel[u*N + v]
    std::vector<double> minTravel_;

    mutable aco::Rng rng_;
};

} // namespace localSearch
