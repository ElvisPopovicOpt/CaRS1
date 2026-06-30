#pragma once

#include <iosfwd>
#include <memory>

namespace cars_tsplib { struct Instance; }

namespace aco { struct Solution; }

namespace solution_io 
{

// Brief: prints node[] and car[] (per edge) on one or two lines.
void printSolutionBrief(std::ostream& os,
                        const aco::Solution& s,
                        double cost,
                        int maxElems = 40);

// Detailed: prints each edge and its cost contribution, including switch return costs
// and the final return. Useful for debugging DP/local search.
void printSolutionDetailed(std::ostream& os,
                           std::shared_ptr<const cars_tsplib::Instance> inst,
                           const aco::Solution& s,
                           double cost);

} // namespace solution_io
