#include <intensifier.hpp>
#include <parser.hpp>
#include <limits>
#include <algorithm>

namespace aco 
{

// Constructor with a DP component
Intensifier::Intensifier(std::shared_ptr<const ILocalSearch> localSearch,
                        std::shared_ptr<const IFixedTourCarAssigner> dp,
                        std::shared_ptr<const cars_tsplib::Instance> inst,
                        int maxStagnationIterations,
                        int resetPheromonesOnDeactivationCount)
    : localSearch_(std::move(localSearch))
    , dp_(std::move(dp))
    , inst_(std::move(inst))
    , maxStagnationIterations_(maxStagnationIterations)
    , resetPheromonesOnDeactivationCount_(resetPheromonesOnDeactivationCount)
{
}

// Constructor without a DP component (LS only)
Intensifier::Intensifier(std::shared_ptr<const ILocalSearch> localSearch,
                        std::shared_ptr<const cars_tsplib::Instance> inst,
                        int maxStagnationIterations,
                        int resetPheromonesOnDeactivationCount)
    : localSearch_(std::move(localSearch))
    , dp_(nullptr)
    , inst_(std::move(inst))
    , maxStagnationIterations_(maxStagnationIterations)
    , resetPheromonesOnDeactivationCount_(resetPheromonesOnDeactivationCount)
{
}

// Requires localSearch_ and a finite best.cost; activates on improvedByIter or while
// within maxStagnationIterations_ after stagnation. On activation, runs LS on a copy
// of best, optionally reassigns cars via DP, and commits the copy back to best only if
// it improves cost by more than improveEps. Exceptions are swallowed since
// intensification is an optional optimization that must not block the algorithm.
double Intensifier::intensify(EvaluatedSolution& best,
                              bool improvedByIter,
                              int itersSinceGlobalImprovement,
                              double improveEps) const
{
    if (!localSearch_ || !std::isfinite(best.cost))
    {
        return std::numeric_limits<double>::infinity();
    }

    // Activate if best was just improved, or we're still within the post-stagnation window
    const bool shouldIntensify = improvedByIter ||
                                (itersSinceGlobalImprovement > 0 &&
                                 itersSinceGlobalImprovement <= maxStagnationIterations_);

    // Count a transition from inactive to active as one activation
    if (shouldIntensify && !wasActiveLastTime_)
    {
        ++activationCount_;
    }

    // wasActiveThisTime_ is set now; wasActiveLastTime_ is updated later via updateState(),
    // after deactivation/pheromone-reset checks have had a chance to run
    wasActiveThisTime_ = shouldIntensify;

    if (!shouldIntensify)
    {
        return std::numeric_limits<double>::infinity();
    }

    // Work on a copy so the original best is untouched until we confirm an improvement
    Solution bestCopy = best.sol;
    double bestCostCopy = best.cost;

    try
    {
        // LS may change the tour (node[]) and cost, but not car[] directly
        localSearch_->improve(bestCopy, bestCostCopy);

        // LS may have altered the tour, so the car assignment may no longer be optimal;
        // DP recomputes it optimally for the (possibly new) tour
        if (dp_ && std::isfinite(bestCostCopy))
        {
            const int N = inst_->n();
            if ((int)bestCopy.car.size() != N)
                bestCopy.car.assign((size_t)N, 0);

            const double newCost = dp_->reassignCars(bestCopy.node, bestCopy.car);
            if (std::isfinite(newCost))
                bestCostCopy = newCost;
        }

        if (std::isfinite(bestCostCopy) && bestCostCopy < best.cost - improveEps)
        {
            best.sol = std::move(bestCopy);
            best.cost = bestCostCopy;
            return bestCostCopy;
        }
    }
    catch (const std::exception& ex)
    {
        // Intensification is an optional optimization; on failure just skip it
    }

    return std::numeric_limits<double>::infinity();
}

bool Intensifier::wasJustDeactivated() const
{
    return wasActiveLastTime_ && !wasActiveThisTime_;
}

bool Intensifier::shouldResetPheromones() const
{
    // Reset only if resets are enabled, the intensifier just deactivated, and we're still
    // within the allowed number of deactivations (e.g. param 1 = reset after the 1st only)
    return resetPheromonesOnDeactivationCount_ > 0 &&
           wasJustDeactivated() &&
           activationCount_ <= resetPheromonesOnDeactivationCount_;
}

void Intensifier::updateState() const
{
    wasActiveLastTime_ = wasActiveThisTime_;
}

} // namespace aco
