#pragma once
#include <vector>
#include <memory>
#include <algorithm>
#include <limits>

// forward declare (ti već imaš Instance u svom projektu)
namespace cars_tsplib { struct Instance; }

namespace aco 
{

class CandidateListCache 
{
public:
    CandidateListCache(std::shared_ptr<const cars_tsplib::Instance> inst, int candidateListSize);

    int N() const { return N_; }
    int K() const { return K_; }

    // lista kandidata za i (bez i samog)
    const std::vector<int>& candidates(int i) const { return cand_[i]; }

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    int N_ = 0;
    int C_ = 0;
    int K_ = 0;

    std::vector<std::vector<int>> cand_;
};

} // namespace aco
