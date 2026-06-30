#pragma once

/**
 * Research module: logs surrogate vs. full DP cost for correlation analysis.
 *
 * Records (surrogate, dp_cost) pairs per constructed solution (per ant per
 * iteration) so the Pearson correlation between the surrogate (min travel
 * cost per edge) and the true cost (CarAssignerDP) can be computed and
 * validated. Intended for small instances and few iterations, since DP
 * evaluateCost is called for every ant.
 */

#include <mutex>
#include <string>
#include <vector>

namespace aco
{

struct SurrogateDpRecord
{
    int runIndex = -1;
    int iteration = -1;
    int antIndex = -1;       // index of the ant within the iteration (0..antsN-1)
    double surrogate = 0.0;  // sum of min travel cost per edge (surrogate for LS)
    double dpCost = 0.0;     // full cost from CarAssignerDP (evaluateCost)
};

class SurrogateCorrelationLogger
{
public:
    SurrogateCorrelationLogger() = default;

    /** Records one (surrogate, dp_cost) pair for the given ant and iteration. */
    void log(int runIdx, int iter, int antIdx, double surrogate, double dpCost);

    /** Writes CSV: run;iteration;ant_index;surrogate;dp_cost. */
    void writeCsv(const std::string& filename, char sep = ';') const;

    /** Number of recorded pairs. */
    size_t numRecords() const;

    /** Thread-safe copy of all records, for later correlation analysis. */
    std::vector<SurrogateDpRecord> getRecords() const;

private:
    mutable std::mutex mtx_;
    std::vector<SurrogateDpRecord> records_;
};

} // namespace aco
