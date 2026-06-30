#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <chrono>
#include "parser.hpp"
#include "ReallocationLocalSearch.hpp"
#include "candidate_list.hpp"

namespace localSearch
{

static inline int prevIdx(int i, int N) { return (i == 0) ? (N - 1) : (i - 1); }
static inline int nextIdx(int i, int N) { return (i + 1 == N) ? 0 : (i + 1); }

struct RelocateBeforeView
{
    const std::vector<int>& base;
    int r; // remove index
    int k; // insert-before index

    int operator()(int idx) const
    {
        if (k < r) {
            if (idx < k) return base[(std::size_t)idx];
            if (idx == k) return base[(std::size_t)r];
            if (idx >= k + 1 && idx <= r) return base[(std::size_t)(idx - 1)];
            return base[(std::size_t)idx];
        } else { // k > r
            if (idx < r) return base[(std::size_t)idx];
            if (idx >= r && idx <= k - 2) return base[(std::size_t)(idx + 1)];
            if (idx == k - 1) return base[(std::size_t)r];
            return base[(std::size_t)idx];
        }
    }
};

static inline void applyRelocateBeforeAndUpdatePos(std::vector<int>& tour,
                                                   std::vector<int>& pos,
                                                   int r, int k)
{
    if (k < r) {
        // move tour[r] to position k by rotating [k..r]
        std::rotate(tour.begin() + k,
                    tour.begin() + r,
                    tour.begin() + (r + 1));
        for (int i = k; i <= r; ++i) pos[tour[(std::size_t)i]] = i;
    } else {
        // move tour[r] to position k-1 by rotating [r..k-1]
        std::rotate(tour.begin() + r,
                    tour.begin() + (r + 1),
                    tour.begin() + k);
        for (int i = r; i <= k - 1; ++i) pos[tour[(std::size_t)i]] = i;
    }
}

// Surrogate: cheapest travel cost over cars between a and b
static inline double minTravel(const cars_tsplib::Instance& inst, int a, int b)
{
    double best = std::numeric_limits<double>::infinity();
    for (int c = 0; c < inst.cars(); ++c)
        best = std::min(best, inst.travelCost(c, a, b));
    return best;
}

// Surrogate delta for insert-before: remove u from (a-u-b), insert between (p-v) giving p-u-v
static inline double surrogateRelocateBeforeDelta(const cars_tsplib::Instance& inst,
                                                  const std::vector<int>& tour,
                                                  int r, int k)
{
    const int N = (int)tour.size();
    const int u = tour[(std::size_t)r];
    const int a = tour[(std::size_t)prevIdx(r, N)];
    const int b = tour[(std::size_t)nextIdx(r, N)];
    const int v = tour[(std::size_t)k];
    const int p = tour[(std::size_t)prevIdx(k, N)];

    const double add = minTravel(inst, a, b) + minTravel(inst, p, u) + minTravel(inst, u, v);
    const double rem = minTravel(inst, a, u) + minTravel(inst, u, b) + minTravel(inst, p, v);
    return add - rem;
}

RelocationLocalSearch::RelocationLocalSearch(std::shared_ptr<const cars_tsplib::Instance> inst,
                                             RelocationOptions opt,
                                             std::shared_ptr<FixedTourCarAssignerDP> dp)
    : inst_(std::move(inst)), opt_(opt), dp_(std::move(dp))
{
    if (!inst_) throw std::runtime_error("RelocationLocalSearch: inst is null");
}

void RelocationLocalSearch::improve(aco::Solution& s, double& cost) const
{
    const int N = inst_->n();
    const bool useCarDP = (dp_ != nullptr);
    if (N < 4) return;
    if ((int)s.node.size() != N) throw std::runtime_error("RelocationLocalSearch: node size mismatch");
    if (s.node[0] != 0) throw std::runtime_error("RelocationLocalSearch: node[0] must be 0");
    if ((int)s.car.size() != N) s.car.assign((std::size_t)N, 0);

    // Optional baseline car reassignment
    if (opt_.reassignCarsAtStart && useCarDP) {
        cost = dp_->reassignCars(s.node, s.car);
    }

    // pos[node] = its index in the tour
    std::vector<int> pos((std::size_t)N, -1);
    for (int i = 0; i < N; ++i) pos[s.node[(std::size_t)i]] = i;

    const bool useCand = (opt_.cand && opt_.cand->N() == N && opt_.cand->K() > 0);

    std::int64_t evals = 0;
    bool anyAccepted = false;

    auto isValidMove = [&](int r, int k) -> bool {
        if (r == 0) return false;     // don't move the depot
        if (k == 0) return false;     // don't insert before the depot
        if (k == r) return false;
        if (k == r + 1) return false; // no-op: u is already right before v
        return true;
    };

    // Stop early if no improvement for several consecutive passes (scaled to instance size)
    int noImproveCount = 0;
    const int maxNoImprove = std::max(10, N / 4);

    for (int pass = 0; pass < opt_.maxPasses; ++pass)
    {
        bool improvedThisPass = false;

        for (int r = 1; r <= N - 1; ++r)
        {
            const int u = s.node[(std::size_t)r];

            auto tryMove = [&](int k) -> bool
            {
                if (!isValidMove(r, k)) return false;

                if (opt_.maxMoveEvaluations > 0 && evals >= opt_.maxMoveEvaluations) return true;

                // Surrogate filter
                const double sur = surrogateRelocateBeforeDelta(*inst_, s.node, r, k);
                if (!(sur < -opt_.minSurrogateGain)) return false;

                double candCost = 0.0;
                if (useCarDP) {
                    RelocateBeforeView view{s.node, r, k};
                    candCost = dp_->evaluateCostViewScratch([&](int idx){ return view(idx); }, scratch_);
                } else {
                    return false; // TODO: TSP fallback without car DP
                }

                ++evals;

                if (std::isfinite(candCost) && candCost + opt_.eps < cost)
                {
                    applyRelocateBeforeAndUpdatePos(s.node, pos, r, k);
                    cost = candCost;
                    anyAccepted = useCarDP;
                    improvedThisPass = true;
                    return opt_.firstImprovement; // true => stop scanning and go to next pass
                }

                return false;
            };

            if (useCand) {
                for (int v : opt_.cand->candidates(u)) {
                    const int k = pos[(std::size_t)v];
                    const bool stop = tryMove(k);
                    if (stop) goto next_pass;
                    if (opt_.firstImprovement && improvedThisPass) goto next_pass;
                }
            } else {
                for (int k = 1; k < N; ++k) { // k=0 not allowed
                    const bool stop = tryMove(k);
                    if (stop) goto next_pass;
                    if (opt_.firstImprovement && improvedThisPass) goto next_pass;
                }
            }
        }

    next_pass:
        if (improvedThisPass) {
            noImproveCount = 0;
        } else {
            noImproveCount++;
            if (noImproveCount >= maxNoImprove) {
                break; // No point searching further
            }
        }
    }

    if (useCarDP && anyAccepted) {
        cost = dp_->reassignCars(s.node, s.car);
    }
}

}