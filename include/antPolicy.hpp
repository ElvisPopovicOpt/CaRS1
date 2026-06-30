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

    // returnFactor_ utječe SAMO na konstrukciju (ant-policy).
    // - returnFactor_ == 0.0: policy uopće ne koristi returnCost() ni tauReturn() (taj kod se preskače).
    // - returnFactor_ > 0.0 : policy može uključiti return heuristiku/feromone u težinu izbora.
    // Pheromone model i dalje deponira tauReturn iz DP-konzistentnih rješenja (to je odvojeno).
    double returnFactor_;

    bool sampleCars_;          // false => gradi samo node permutaciju
    bool enforceNoRerent_;     // true => zabrani ponovno iznajmljivanje auta nakon povrata
    double q0_; 
};

} // namespace aco
