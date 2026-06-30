#pragma once
#include <memory>
#include <vector>
#include <interfaces.hpp>
#include <candidate_list.hpp>

namespace aco 
{

class AntPolicyCandidateListRoulette final : public IAntPolicy 
{
public:
    AntPolicyCandidateListRoulette(std::shared_ptr<const CandidateListCache> cl,
                                   double epsilonEta,
                                   double returnFactor = 1.0,
                                   bool sampleCars = true, 
                                   bool enforceNoRerent = false,
                                   double q0 = 0.0);

    Solution construct(const AntContext& ctx) const override;

private:
    std::shared_ptr<const CandidateListCache> cl_;
    double eps_;

    // returnFactor_ only affects construction (ant policy):
    // 0.0 skips returnCost()/tauReturn() entirely; >0.0 lets return heuristic/pheromones
    // influence selection weight. Pheromone deposit still uses tauReturn from DP-consistent
    // solutions independently of this.
    double returnFactor_;

    bool sampleCars_;          // false => build node permutation only
    bool enforceNoRerent_;     // true => forbid re-renting a car after it has been returned
    double q0_;
};

} // namespace aco
