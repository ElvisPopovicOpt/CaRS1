#pragma once

#include <memory>
#include <cstdint>
#include <limits>
#include <vector>
#include "interfaces.hpp"
#include "FixedTourCarAssignerDP.hpp" 

namespace cars_tsplib { struct Instance; }
namespace localSearch { class FixedTourCarAssignerDP; }
namespace aco { class CandidateListCache; }

namespace localSearch
{

struct RelocationOptions
{
    int maxPasses = 10;
    bool firstImprovement = true;
    std::int64_t maxMoveEvaluations = 2000;
    double timeLimitMs = 0.0;
    double eps = 1e-12;

    bool reassignCarsAtStart = true;

    // Candidate list, indexed per node
    std::shared_ptr<const aco::CandidateListCache> cand = nullptr;

    // Surrogate filter threshold: 0 = only call DP when surrogate shows strict improvement
    double minSurrogateGain = 0.0;
};

class RelocationLocalSearch final : public aco::ILocalSearch
{
public:
    RelocationLocalSearch(std::shared_ptr<const cars_tsplib::Instance> inst,
                          RelocationOptions opt,
                          std::shared_ptr<FixedTourCarAssignerDP> dp = nullptr);

    void improve(aco::Solution& s, double& cost) const override;

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    RelocationOptions opt_;
    std::shared_ptr<FixedTourCarAssignerDP> dp_;

    // Scratch buffer for DP view-evaluation, avoids reallocation
    mutable localSearch::Scratch scratch_;
};

class ChainedLocalSearch final : public aco::ILocalSearch
{
public:
    explicit ChainedLocalSearch(std::vector<std::shared_ptr<const aco::ILocalSearch>> seq)
        : seq_(std::move(seq)) {}

    void improve(aco::Solution& s, double& cost) const override
    {
        for (const auto& ls : seq_) {
            ls->improve(s, cost);
        }
    }

private:
    std::vector<std::shared_ptr<const aco::ILocalSearch>> seq_;
};


} // namespace localSearch
