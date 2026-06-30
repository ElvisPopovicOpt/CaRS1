#pragma once

#include <memory>
#include "interfaces.hpp"

namespace cars_tsplib { struct Instance; }
namespace localSearch { class FixedTourCarAssignerDP; }
namespace aco { class CandidateListCache; }

namespace localSearch
{

// Reinsertion Local Search: single-node relocation with optional candidate list.
// Same algorithm as in ACO_CarSP (2-opt -> Or-opt -> Reinsertion -> 2-opt), CaRS variant
// has no passengers; only improve(Solution&, double&) is used.
struct ReinsertionOptions
{
    int maxPasses = 5;
    bool firstImprovement = true;
};

class ReinsertionLocalSearch final : public aco::ILocalSearch
{
public:
    ReinsertionLocalSearch(std::shared_ptr<const cars_tsplib::Instance> inst,
                           std::shared_ptr<FixedTourCarAssignerDP> dp,
                           std::shared_ptr<const aco::CandidateListCache> cand,
                           ReinsertionOptions opt = {});

    void improve(aco::Solution& s, double& cost) const override;

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    std::shared_ptr<FixedTourCarAssignerDP> dp_;
    std::shared_ptr<const aco::CandidateListCache> cand_;
    ReinsertionOptions opt_;

    double evaluateReinsertionMove(const std::vector<int>& nodes, int nodeIdx, int newPos) const;
    void applyReinsertionMove(std::vector<int>& nodes, int nodeIdx, int newPos) const;
};

} // namespace localSearch
