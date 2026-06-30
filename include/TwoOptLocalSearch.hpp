#pragma once

#include <memory>
#include <optional>
#include <vector>
#include <cstdint>

#include "interfaces.hpp"              // defines aco::ILocalSearch and aco::Solution
#include "FixedTourCarAssignerDP.hpp"  // needed for the pointer and evaluateCostView

namespace cars_tsplib { struct Instance; }


namespace aco {
class CandidateListCache;
class ILocalSearch;
struct Solution;
} // namespace aco

namespace localSearch
{

struct CarAssignmentDPOptions;
class FixedTourCarAssignerDP;

struct TwoOptOptions
{
    int maxPasses = 20;
    bool firstImprovement = true;
    std::int64_t maxMoveEvaluations = -1; // -1 = unlimited
    double timeLimitMs = 0.0;
    double eps = 1e-12;
    bool reassignCarsAtStart = true;
    CarAssignmentDPOptions dpOpt;

    std::shared_ptr<const aco::CandidateListCache> cand = nullptr;

    // Don't-look bits (DLB)
    bool useDontLookBits = false;

    // How many positions around each cut to wake after an accepted move (0 = cuts only). Typically 1 or 2.
    int dontLookWakeRadius = 1;
};

class TwoOptLocalSearch final : public aco::ILocalSearch
{
public:
    TwoOptLocalSearch(std::shared_ptr<const cars_tsplib::Instance> inst,
                      TwoOptOptions opt = {},
                      std::shared_ptr<FixedTourCarAssignerDP> dp = nullptr);

    void improve(aco::Solution& s, double& cost) const override;

private:
    // --- core helpers ---
    bool timeExceeded_(std::uint64_t startNs) const;
    double singleCarCycleCost_(const std::vector<int>& nodes) const;

    bool likelySymmetricTsp_() const; // distinguishes symmetric vs asymmetric TSP

    // --- per-branch implementations ---
    void improveFirst_(aco::Solution& s,
                       double& cost,
                       bool useCarDP,
                       bool useCand,
                       bool useSymTspDelta,   
                       std::vector<int>& pos,
                       std::vector<std::uint8_t>& dontLookNode,
                       std::int64_t& evals,
                       std::uint64_t t0,
                       bool& anyAccepted) const;

    void improveBest_(aco::Solution& s,
                      double& cost,
                      bool useCarDP,
                      bool useCand,
                      bool useSymTspDelta,
                      std::vector<int>& pos,
                      std::vector<std::uint8_t>& dontLookNode,
                      std::int64_t& evals,
                      std::uint64_t t0,
                      bool& anyAccepted) const;

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    TwoOptOptions opt_;

    // If nullptr, behavior falls back to plain TSP.
    std::shared_ptr<FixedTourCarAssignerDP> dp_;

    // Scratch buffer for DP view evaluation (avoids realloc), same pattern as RelocationLocalSearch.
    mutable Scratch scratch_;
};

} // namespace localSearch


