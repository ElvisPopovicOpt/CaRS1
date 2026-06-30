#pragma once
#include <solution.hpp>
#include <rng.hpp>  
#include <mutex>
#include <memory>
#include <limits>

namespace cars_tsplib { struct Instance; }
namespace aco_cli { struct ParamsData; }

namespace aco 
{

struct IPheromoneAccessor
{
    virtual ~IPheromoneAccessor() = default;

    // MOVE pheromones (always present)
    virtual double tauMove(int car, int i, int j) const = 0;

    // RETURN pheromones (optional). Default: none -> neutral 1.0 (no effect on weights).
    virtual double tauReturn(int /*car*/, int /*from*/, int /*to*/) const { return 1.0; }

    // Lets ant policies know whether return pheromones are meaningful
    virtual bool hasReturnPheromones() const { return false; }
};

class UniformPheromoneAccessor final : public IPheromoneAccessor
{
public:
    explicit UniformPheromoneAccessor(double tau = 1.0) : tau_(tau) {}
    double tauMove(int, int, int) const override { return tau_; }

    // uniform accessor has no return pheromones
    double tauReturn(int, int, int) const override { return 1.0; }
    bool hasReturnPheromones() const override { return false; }

private:
    double tau_;
};


struct IPheromoneModel : public IPheromoneAccessor
{
    virtual ~IPheromoneModel() = default;

    virtual void reset() = 0;

    virtual void onIterationEnd(int iterationIndex,
                                const EvaluatedSolution& iterBest,
                                const EvaluatedSolution& globalBest,
                                bool globalImproved) = 0;

    // Ranked (top-W) update. Default implementation just uses ranked[0] as iterBest.
    // Colony currently ties W to eliteKAnts (a dedicated parameter could be added later).
    virtual void onIterationEndRanked(int iterationIndex,
                                      const std::vector<EvaluatedSolution>& ranked,
                                      const EvaluatedSolution& globalBest,
                                      bool globalImproved)
    {
        if (ranked.empty()) return;
        onIterationEnd(iterationIndex, ranked.front(), globalBest, globalImproved);
    }

    // Check if smoothing was applied in the last iteration (default: false)
    virtual bool wasSmoothingApplied() const { return false; }
    
    // Get stagnation counter (noImprove_) for smoothing calculation (default: 0)
    virtual int getStagnationCounter() const { return 0; }
    
    /**
     * Selectively reduces pheromones on a tour's edges (smoothing) instead of a full
     * reset: tau = (1 - gamma) * tau + gamma * tauBase.
     *
     * @param solution Solution whose tour edges should be smoothed
     * @param gamma Reduction factor (0.0 = fully reset to tauBase, 1.0 = no change)
     */
    virtual void reducePheromonesOnTour(const Solution& solution, double gamma) { (void)solution; (void)gamma; }

    // Adds a solution to the archive (no-op if archiving is disabled)
    virtual void addToArchive(const EvaluatedSolution& solution) { (void)solution; }

    // @return Best archived cost, or infinity if archiving is disabled/empty
    virtual double getBestArchivedCost() const { return std::numeric_limits<double>::infinity(); }

    // @return true and fills out with the best archived solution, if any exists
    virtual bool getBestArchivedSolution(EvaluatedSolution& out) const { (void)out; return false; }
    
};


struct ICostModel 
{
    virtual ~ICostModel() = default;
    virtual double evaluate(const Solution& s) const = 0;
};

struct ILocalSearch 
{
    virtual ~ILocalSearch() = default;
    virtual void improve(Solution& s, double& cost) const = 0;
    // optional seeding (default: no-op)
    virtual void resetSeed(uint64_t /*seed*/) const {}
};



struct AntContext 
{
    std::shared_ptr<const cars_tsplib::Instance> inst;
    const aco_cli::ParamsData* params = nullptr;
    const IPheromoneAccessor* pher = nullptr;

    // per-run RNG
    Rng* rng = nullptr;

    // Dynamic q0 override: -1.0 = use the ant policy's fixed q0, >=0.0 = override q0
    double adaptiveQ0 = -1.0;
};

struct IAntPolicy 
{
    virtual ~IAntPolicy() = default;
    virtual Solution construct(const AntContext& ctx) const = 0;
};

struct IFixedTourCarAssigner
{
    virtual ~IFixedTourCarAssigner() = default;

    // For a given Hamiltonian cycle (nodes[0]==0, size N), fills carPerEdgeOut (size N)
    // and returns the minimal cost
    virtual double reassignCars(const std::vector<int>& nodes,
                                std::vector<int>& carPerEdgeOut) const = 0;

    // Cost-only variant
    virtual double evaluateCost(const std::vector<int>& nodes) const = 0;
};

/**
 * Interface for the intensifier component, which applies an extra Local Search pass
 * to the global-best solution for deeper exploitation without much added cost.
 *
 * Activates when best was just updated this iteration (improvedByIter == true), or
 * during the first few iterations after stagnation.
 */
struct IIntensifier
{
    virtual ~IIntensifier() = default;

    /**
     * Applies intensification to the best solution.
     *
     * @param best Global best solution (updated in place if intensification improves it)
     * @param improvedByIter Whether best was just updated this iteration
     * @param itersSinceGlobalImprovement Iterations since the last global improvement
     * @param improveEps Epsilon used to judge whether a result counts as an improvement
     * @return The improved cost, or infinity if no improvement was found
     */
    virtual double intensify(EvaluatedSolution& best,
                            bool improvedByIter,
                            int itersSinceGlobalImprovement,
                            double improveEps) const = 0;

    /**
     * True if the intensifier was just deactivated and activationCount_ is within
     * the allowed threshold (e.g. param 1 = reset after the 1st deactivation, 2 = after
     * the 1st and 2nd). Threshold <= 0 means never reset.
     */
    virtual bool shouldResetPheromones() const = 0;

    // @return True if the intensifier was active last call but is not active now
    virtual bool wasJustDeactivated() const = 0;

    // @return Number of times the intensifier has switched on
    virtual int getActivationCount() const = 0;

    // Updates internal state for the next call; call after deactivation/reset checks.
    virtual void updateState() const = 0;
};

} // namespace aco
