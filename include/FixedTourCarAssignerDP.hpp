#pragma once

#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cars_tsplib { struct Instance; }

namespace localSearch 
{

// Result of the DP: minimal cost and the optimal car per edge (size N).
struct CarAssignmentDPResult
{
    double cost = std::numeric_limits<double>::infinity();
    std::vector<int> carPerEdge; // size = N
};

struct CarAssignmentDPOptions
{
    // If true: verify that tour is a permutation of [0..N-1] and node[0]==0.
    // Useful while developing 2-opt/LK (catches silent bugs).
    bool validateTour = false;

    // Heuristic cap on segment length (for speed). -1 => no limit (exact optimum).
    int maxSegmentLen = -1;

    // Robust fallback: if maxSegmentLen > 0 and the DP returns INF, automatically retry with a larger L, then finally with -1.
    bool fallbackToUnlimitedOnInf = true;

    // How many times to grow L before the final -1 attempt (if enabled).
    // Example: maxSegmentLen=20, growth=2, retries=3 => 20,40,80, then -1.
    int fallbackMaxRetries = 3;

    // Growth factor for the segment length on retry (>=2 recommended).
    int fallbackGrowthFactor = 2;
};

// Avoids allocating vectors during execution.
// N - number of nodes, C - number of cars
struct Scratch
{
    int N = 0;
    int C = 0;
    int M = 0; // 1<<C

    // pref[c*(N+1) + t]
    std::vector<double> pref;

    // dp[t*M + mask]
    std::vector<double> dp;

    // Materialized tour (N+1, includes closure tour[N]=tour[0])
    std::vector<int> tour;

    // For validateTour without allocations/clearing.
    std::vector<int> seenStamp;
    int stamp = 1;

    void ensure(int N_, int C_);
};




// DP optimizer that, for a fixed Hamiltonian cycle of nodes (node[0]==0),
// finds the optimal segmentation across cars under CaRS "switch-return" semantics.
//
// Semantics:
// - segment [s, t) uses car c on edges s..t-1
// - segment cost is:
//      sum travelCost(c, v[k], v[k+1]) for k=s..t-1
//    + returnCost(c, v[t], v[s])   // closing the segment (switch/final)
//
// This way the DP naturally models the last car's final return too, since v[N]=v[0].

class FixedTourCarAssignerDP
{
public:
    explicit FixedTourCarAssignerDP(std::shared_ptr<const cars_tsplib::Instance> inst,
                                    CarAssignmentDPOptions opt = {});

    CarAssignmentDPResult optimize(const std::vector<int>& nodes) const;

    double reassignCars(const std::vector<int>& nodes,
                        std::vector<int>& carPerEdgeOut) const;

    double evaluateCost(const std::vector<int>& nodes) const;

    // Cost-only evaluation of a tour via a "view".
    // getNode(idx) must return the node at position idx (0..N-1); the cycle is closed internally.
    double evaluateCostView(const std::function<int(int)>& getNode) const;

    // Same as evaluateCostView but uses a caller-supplied scratch (no reallocation).
    double evaluateCostViewScratch(const std::function<int(int)>& getNode, Scratch& scratch) const;

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    CarAssignmentDPOptions opt_;

    void validateTourOrThrow_(const std::vector<int>& nodes) const;

    // Internal: single run without fallback, with an explicit maxSegmentLen parameter.
    CarAssignmentDPResult optimizeWithMaxSegLen_(const std::vector<int>& nodes,
                                                int maxSegmentLen) const;
    
    double evaluateCostViewWithMaxSegLen_(const std::function<int(int)>& getNode,
                                         int maxSegmentLen) const;

    // Scratch-based variant of evaluateCostViewWithMaxSegLen_.
    double evaluateCostViewWithMaxSegLenScratch_(const std::function<int(int)>& getNode,
                                                int maxSegmentLen,
                                                Scratch& scratch) const;

    // Internal: fallback wrapper for optimize() / evaluateCostView()
    CarAssignmentDPResult optimizeWithFallback_(const std::vector<int>& nodes) const;

    double evaluateCostViewWithFallback_(const std::function<int(int)>& getNode) const;

    int normalizeMaxSegLen_(int L) const;

    // Per-instance scratch buffer; mutable because evaluateCostView is const.
    mutable Scratch scratch_;
};

} // namespace localSearch
