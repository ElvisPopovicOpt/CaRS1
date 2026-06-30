#include <stdexcept>
#include <limits>
#include <iostream>
#include <parser.hpp> 
#include <cost_model.hpp>

namespace aco 
{

CostModelCars::CostModelCars(std::shared_ptr<const cars_tsplib::Instance> inst)
    : inst_(std::move(inst))
{
    if (!inst_) throw std::runtime_error("CostModelCars: inst is null");
}

double CostModelCars::evaluate(const Solution& s) const {
    const int N = inst_->n();
    if ((int)s.node.size() != N || (int)s.car.size() != N) {
        throw std::runtime_error("CostModelCars: Solution node/car size mismatch with instance N.");
    }
    if (N <= 0) return 0.0;

    // start fiksan
    if (s.node[0] != 0) 
    {
        throw std::runtime_error("CostModelCars: node[0] must be 0 (fixed start).");
    }

    double result = 0.0;

    int currentCar = s.car[0];
    int lastCarNode = s.node[0];
    for (int i = 0; i < N; ++i) {
        const int currentNode = s.node[i];
        const int nextNode    = s.node[(i + 1) % N];

        const int lastCar = (i > 0) ? s.car[i - 1] : currentCar;
        currentCar = s.car[i];
        // CaRS => passNumber=0 => / (passNumber+1) == /1, pa ne dijelimo
        if (currentCar != lastCar) 
        {
            result += inst_->returnCost(lastCar, currentNode, lastCarNode);
            lastCarNode = currentNode;
        }

        result += inst_->travelCost(currentCar, currentNode, nextNode);
    }
    // finalni return trenutnog auta u start node-u prema lastCarNode
    // ako je samo jedan auto to je tsp, a i mora cvor vracanja biti isti kao i cvor iznajmljivanja pa se ne naplacuje
    // vec je u returnCost funkciji ispitano hasReturnCost i vraca nulu ako nema matrice returnCosts
    result += inst_->returnCost(currentCar, s.node[0], lastCarNode);

    return result;
}

} // namespace aco
