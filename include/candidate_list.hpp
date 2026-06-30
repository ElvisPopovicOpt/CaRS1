#pragma once
#include <vector>
#include <memory>
#include <algorithm>
#include <limits>

// forward declaration of Instance
namespace cars_tsplib { struct Instance; }

namespace aco 
{

class CandidateListCache 
{
public:
    CandidateListCache(std::shared_ptr<const cars_tsplib::Instance> inst, int candidateListSize);

    int N() const { return N_; }
    int K() const { return K_; }

    // candidate list for node i (excluding i itself)
    const std::vector<int>& candidates(int i) const { return cand_[i]; }

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    int N_ = 0;
    int C_ = 0;
    int K_ = 0;

    std::vector<std::vector<int>> cand_;
};

} // namespace aco
