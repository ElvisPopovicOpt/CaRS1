#pragma once
#include <interfaces.hpp>
#include <memory>

namespace cars_tsplib { struct Instance; }

namespace aco 
{

/**
 * Intensifier: applies an extra Local Search pass to the global-best solution for
 * deeper exploitation without much added cost (it only ever works on the best solution).
 *
 * Activates when:
 * 1. The best was just improved this iteration (improvedByIter == true), or
 * 2. We're within the first maxStagnationIterations_ iterations after stagnation
 *    (default 3, overridable via constructor) — gives a few extra chances to improve
 *    even without a fresh global improvement.
 *
 * Wiring: created in main.cpp's factory lambda before the Colony, then passed to the
 * Colony constructor (optional, default nullptr). Colony will auto-create one if none
 * is supplied but the required components (localSearch_, dp_) exist. Used inside
 * Colony::run() after rr.best is updated with iterBestCost.
 *
 * Flow: check activation conditions -> copy best and run an LS pass -> if DP is
 * available, reassign cars optimally for the (possibly LS-modified) tour -> if the
 * result improves on best.cost by more than improveEps, update best and return the
 * new cost; otherwise return infinity.
 */
class Intensifier final : public IIntensifier
{
public:
    /**
     * Constructor with a DP component.
     *
     * @param localSearch Local Search component used to improve the solution
     * @param dp Dynamic Programming component for optimal car assignment
     * @param inst Problem instance
     * @param maxStagnationIterations Iterations after stagnation during which intensification still triggers (default: 3)
     * @param resetPheromonesOnDeactivationCount Number of deactivations after which pheromones are reset (<=0 = never, 1 = after 1st deactivation, 2 = after the 1st and 2nd, etc.; default: 1)
     */
    Intensifier(std::shared_ptr<const ILocalSearch> localSearch,
                std::shared_ptr<const IFixedTourCarAssigner> dp,
                std::shared_ptr<const cars_tsplib::Instance> inst,
                int maxStagnationIterations = 3,
                int resetPheromonesOnDeactivationCount = 1);
    
    /**
     * Constructor without a DP component (LS only).
     *
     * @param localSearch Local Search component used to improve the solution
     * @param inst Problem instance
     * @param maxStagnationIterations Iterations after stagnation during which intensification still triggers (default: 3)
     * @param resetPheromonesOnDeactivationCount Number of deactivations after which pheromones are reset (<=0 = never, 1 = after 1st deactivation, etc.; default: 1)
     */
    Intensifier(std::shared_ptr<const ILocalSearch> localSearch,
                std::shared_ptr<const cars_tsplib::Instance> inst,
                int maxStagnationIterations = 3,
                int resetPheromonesOnDeactivationCount = 1);
    
    /**
     * Applies intensification to the best solution.
     *
     * @param best Global best solution (updated in place if intensification improves it)
     * @param improvedByIter Whether best was just updated this iteration
     * @param itersSinceGlobalImprovement Iterations since the last global improvement
     * @param improveEps Epsilon used to judge whether a result counts as an improvement
     * @return The improved cost, or infinity if no improvement was found
     */
    double intensify(EvaluatedSolution& best,
                    bool improvedByIter,
                    int itersSinceGlobalImprovement,
                    double improveEps) const override;

    /**
     * Returns how many times the intensifier has switched on (incremented each
     * time it transitions from inactive to active).
     *
     * @return Number of activations
     */
    int getActivationCount() const override { return activationCount_; }

    /**
     * @return True if the intensifier was active last call but is not active now
     */
    bool wasJustDeactivated() const override;

    /**
     * True if the intensifier was just deactivated, resetPheromonesOnDeactivationCount_ > 0,
     * and activationCount_ <= resetPheromonesOnDeactivationCount_ (e.g. parameter 1 = reset
     * only after the 1st deactivation, 2 = after the 1st and 2nd).
     *
     * @return True if pheromones should be reset
     */
    bool shouldResetPheromones() const override;

    /**
     * Updates internal state for the next call. Must be called after deactivation
     * and pheromone-reset checks have been made.
     */
    void updateState() const override;

private:
    std::shared_ptr<const ILocalSearch> localSearch_;   // LS component used for improvement
    std::shared_ptr<const IFixedTourCarAssigner> dp_;   // optional DP component for car assignment
    std::shared_ptr<const cars_tsplib::Instance> inst_; // problem instance
    int maxStagnationIterations_;            // iterations after stagnation during which intensification still triggers
    int resetPheromonesOnDeactivationCount_; // number of deactivations after which pheromones are reset (<=0 = never)

    mutable int activationCount_ = 0;         // number of times the intensifier has switched on
    mutable bool wasActiveLastTime_ = false;  // active on the previous call (for detecting inactive->active transitions)
    mutable bool wasActiveThisTime_ = false;  // active on this call (for detecting deactivation)
};

} // namespace aco
