#include <numeric>
#include <stdexcept>
#include <parser.hpp>
#include <candidate_list.hpp>

namespace aco 
{

CandidateListCache::CandidateListCache(
    std::shared_ptr<const cars_tsplib::Instance> inst, int candidateListSize)
    : inst_(std::move(inst))
{
    if (!inst_) throw std::runtime_error("CandidateListCache: inst is null");
    N_ = inst_->n();
    C_ = inst_->cars();
    K_ = std::max(0, candidateListSize);

    cand_.assign(static_cast<size_t>(N_), {});

    // For each node i, sort candidates j by min-over-cars travelCost(car,i,j) and keep top-K
    std::vector<int> nodes;
    nodes.resize(static_cast<size_t>(N_));

    for (int i = 0; i < N_; ++i) 
    {
        std::iota(nodes.begin(), nodes.end(), 0);

        auto score = [&](int j) -> double 
        {
            if (j == i) return std::numeric_limits<double>::infinity();
            double best = std::numeric_limits<double>::infinity();
            for (int c = 0; c < C_; ++c) 
            {
                best = std::min(best, inst_->travelCost(c, i, j));
            }
            return best;
        };

        std::partial_sort(nodes.begin(), nodes.begin() + std::min(N_, K_ + 1), nodes.end(),
            [&](int a, int b) { return score(a) < score(b); });

        auto& v = cand_[static_cast<size_t>(i)];
        v.clear();
        v.reserve(static_cast<size_t>(std::min(K_, N_ - 1)));
        for (int j : nodes) 
        {
            if (j == i) continue;
            v.push_back(j);
            if ((int)v.size() >= K_) break;
        }
    }
}

} // namespace aco
