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
                                   double epsilonEta = 1e-12);

    Solution construct(const AntContext& ctx) const override;

private:
    std::shared_ptr<const CandidateListCache> cl_;
    double eps_;
};

} // namespace aco
