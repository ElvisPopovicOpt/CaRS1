#pragma once

#include <memory>
#include <cstdint>

namespace cars_tsplib { struct Instance; }

namespace aco 
{
struct Solution;
class CostModelCars;
}

namespace selftest 
{

// Throws std::runtime_error if the test fails.
void runDpConsistencyTest(std::shared_ptr<const cars_tsplib::Instance> inst,
                          std::uint64_t seed,
                          int toursToTest = 200,
                          double eps = 1e-9);

} // namespace selftest
