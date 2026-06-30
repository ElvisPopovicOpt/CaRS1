#pragma once
#include <interfaces.hpp>

namespace aco { struct Solution; }

namespace localSearch 
{

struct NoLocalSearch final : aco::ILocalSearch
{
    void improve(aco::Solution&, double&) const override {}
};

} // namespace localSearch
