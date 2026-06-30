#pragma once
#include <memory>
#include <interfaces.hpp> // ICostModel + Solution

namespace cars_tsplib { struct Instance; }

namespace aco 
{

// Cost model za CaRS (bez putnika).
// Semantika je 1:1 s tvojim referentnim CalculateCostLikeThis (passNumber=0 => /1).
class CostModelCars final : public ICostModel 
{
public:
    explicit CostModelCars(std::shared_ptr<const cars_tsplib::Instance> inst);

    double evaluate(const Solution& s) const override;

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
};

} // namespace aco
