#pragma once
#include <memory>
#include <interfaces.hpp> // ICostModel + Solution

namespace cars_tsplib { struct Instance; }

namespace aco
{

// Cost model for CaRS (no passengers); matches the reference cost calculation with passNumber=0 (i.e. no division).
class CostModelCars final : public ICostModel
{
public:
    explicit CostModelCars(std::shared_ptr<const cars_tsplib::Instance> inst);

    double evaluate(const Solution& s) const override;

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
};

} // namespace aco
