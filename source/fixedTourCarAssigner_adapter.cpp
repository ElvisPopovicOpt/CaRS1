#include "fixedTourCarAssigner_adapter.hpp"
#include <FixedTourCarAssignerDP.hpp> // tvoj pravi header (ime prilagodi)

namespace aco 
{

double FixedTourCarAssignerDPAdapter::reassignCars(const std::vector<int>& nodes,
                                                   std::vector<int>& carPerEdgeOut) const
{
    return impl_->reassignCars(nodes, carPerEdgeOut);
}

double FixedTourCarAssignerDPAdapter::evaluateCost(const std::vector<int>& nodes) const
{
    return impl_->evaluateCost(nodes);
}

} // namespace aco
