#pragma once
#include <memory>
#include <stdexcept>
#include <vector>
#include <interfaces.hpp>

// forward declaration da header ostane lagan
namespace localSearch { class FixedTourCarAssignerDP; }

namespace aco 
{

class FixedTourCarAssignerDPAdapter final : public IFixedTourCarAssigner
{
public:
    explicit FixedTourCarAssignerDPAdapter(
        std::shared_ptr<const localSearch::FixedTourCarAssignerDP> impl)
        : impl_(std::move(impl))
    {
        if (!impl_) throw std::runtime_error("FixedTourCarAssignerDPAdapter: impl is null");
    }

    double reassignCars(const std::vector<int>& nodes,
                        std::vector<int>& carPerEdgeOut) const override;

    double evaluateCost(const std::vector<int>& nodes) const override;

private:
    std::shared_ptr<const localSearch::FixedTourCarAssignerDP> impl_;
};

} // namespace aco
