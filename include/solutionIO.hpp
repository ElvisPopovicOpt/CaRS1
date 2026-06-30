#pragma once

#include <iosfwd>
#include <memory>

namespace cars_tsplib { struct Instance; }

namespace aco { struct Solution; }

namespace solution_io 
{

// Kratko: ispiše node[] i car[] (po edge-u) u jednoj liniji (ili više).
void printSolutionBrief(std::ostream& os,
                        const aco::Solution& s,
                        double cost,
                        int maxElems = 40);

// Detaljno: ispiše svaki edge i doprinos trošku uključujući return na switch i finalni return.
// Ovo je super za debug DP/LS.
void printSolutionDetailed(std::ostream& os,
                           std::shared_ptr<const cars_tsplib::Instance> inst,
                           const aco::Solution& s,
                           double cost);

} // namespace solution_io
