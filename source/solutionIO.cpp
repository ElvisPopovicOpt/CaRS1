#include <ostream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>

#include <parser.hpp>         // cars_tsplib::Instance (prilagodi include)
#include <solution.hpp>    // aco::Solution (prilagodi include)
#include <solutionIO.hpp>

namespace solution_io 
{

void printSolutionBrief(std::ostream& os,
                        const aco::Solution& s,
                        double cost,
                        int maxElems)
{
    const int N = static_cast<int>(s.node.size());
    os << "cost=" << cost << ", N=" << N << "\n";
    os << "nodes: ";
    for (int i = 0; i < N && i < maxElems; ++i) os << s.node[static_cast<std::size_t>(i)] << ' ';
    if (N > maxElems) os << "...";
    os << "\n";

    os << "cars : ";
    for (int i = 0; i < static_cast<int>(s.car.size()) && i < maxElems; ++i) os << s.car[static_cast<std::size_t>(i)] << ' ';
    if (static_cast<int>(s.car.size()) > maxElems) os << "...";
    os << "\n";
}

void printSolutionDetailed(std::ostream& os,
                           std::shared_ptr<const cars_tsplib::Instance> inst,
                           const aco::Solution& s,
                           double cost)
{
    if (!inst) throw std::runtime_error("printSolutionDetailed: inst is null");

    const int N = inst->n();
    if ((int)s.node.size() != N) {
        throw std::runtime_error("printSolutionDetailed: solution node size mismatch with instance N");
    }

    // Dopusti TSP/single-car slučaj gdje car vektor može biti prazan ili krive veličine:
    const bool hasCarsVector = ((int)s.car.size() == N);
    const bool hasRet = inst->hasReturnCosts();

    if (N <= 0) {
        os << "Empty solution.\n";
        return;
    }

    os << "=== Solution detailed ===\n";
    os << "Total cost (reported): " << std::setprecision(15) << cost << "\n";
    os << "Tour start fixed: node[0]" << "\n";
    os << "Return costs present: " << (hasRet ? "YES" : "NO (TSP-like)") << "\n\n";

    double accum = 0.0;

    int currentCar = hasCarsVector ? s.car[0] : 0;
    int lastCarNode = s.node[0];

    for (int i = 0; i < N; ++i)
    {
        const int u = s.node[i];
        const int v = s.node[(i + 1) % N];

        const int lastCar = (i > 0)
            ? (hasCarsVector ? s.car[i - 1] : 0)
            : currentCar;

        currentCar = hasCarsVector ? s.car[i] : 0;

        // Switch trošak ima smisla samo ako postoje return costs.
        if (hasRet && currentCar != lastCar) {
            const double ret = inst->returnCost(lastCar, u, lastCarNode);
            accum += ret;
            os << "[switch] at node " << u
               << " return car " << lastCar
               << " from " << u << " -> " << lastCarNode
               << " : +" << ret << " (acc=" << accum << ")\n";
            lastCarNode = u;
        }

        const double tr = inst->travelCost(currentCar, u, v);
        accum += tr;

        os << "edge " << i << ": " << u << " -> " << v
           << "  car=" << currentCar
           << "  travel=+" << tr
           << "  (acc=" << accum << ")\n";
    }

    // Final return ima smisla samo ako postoje return costs.
    if (hasRet) 
    {
        const double finalRet = inst->returnCost(currentCar, s.node[0], lastCarNode);
        accum += finalRet;
        os << "[final] return car " << currentCar
           << " from " << s.node[0] << " -> " << lastCarNode
           << " : +" << finalRet << " (acc=" << accum << ")\n";
    }
    if (!hasRet) 
    {
        for (int c : s.car) if (c != 0)
            throw std::runtime_error("printSolutionDetailed: TSP instance expects all cars == 0");
    }


    os << "Total cost (reconstructed): " << std::setprecision(15) << accum << "\n";
    os << "=========================\n";
}


} // namespace solution_io
