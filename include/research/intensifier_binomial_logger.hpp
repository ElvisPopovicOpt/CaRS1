#pragma once

/**
 * Research module: logs intensifier activations for binomial analysis.
 *
 * Records per activation: cost at start, cost at deactivation, whether the
 * best solution improved, and depth (number of global-improvement stagnation
 * resets that occurred while the intensifier was active). From this, per-run
 * (n, X) pairs can be derived for a Bin(n, p) model, where n = number of
 * activations and X = number of successes (improvements).
 *
 * Used for ablations and statistical validation of the intensifier's
 * binomial behavior. CSV output is written for later analysis.
 */

#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace aco
{

struct IntensifierActivationRecord
{
    int runIndex = -1;
    int activationIndex = 0;   // 1-based index of the activation within the run
    double costAtStart = 0.0;  // best cost at the start of this activation
    double costAtEnd = 0.0;    // best cost at deactivation
    bool improved = false;     // whether the best solution improved during this activation
    int depth = 0;             // number of global-improvement stagnation resets during this activation
};

class IntensifierBinomialLogger
{
public:
    IntensifierBinomialLogger() = default;

    /** Called when the intensifier starts a new activation. */
    void onActivationStart(int runIdx, int activationIndex, double costAtStart);

    /** Resets the depth counter for a run at the start of a new activation. */
    void resetActivationDepth(int runIdx);

    /** Records a global-improvement stagnation reset while the intensifier is active; increments depth for the run's current activation. */
    void noteGlobalStagnationResetDuringIntensifier(int runIdx);

    /** Called on each intensifier deactivation, recording the cost at start/end and whether the best solution improved. */
    void onDeactivation(int runIdx, int activationIndex,
                        double costAtStart, double costAtEnd, bool improved);

    /** Writes CSV: run;activation_index;cost_start;cost_end;improved;depth (one row per deactivation). */
    void writeCsv(const std::string& filename, char sep = ';') const;

    /** Returns the number of runs with at least one recorded activation. */
    size_t numRunsWithData() const;

    /** Returns all records (thread-safe copy). */
    std::vector<IntensifierActivationRecord> getRecords() const;

private:
    mutable std::mutex mtx_;
    std::vector<IntensifierActivationRecord> records_;
    std::unordered_map<int, int> activationDepthByRun_;
};

} // namespace aco
