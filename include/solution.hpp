#pragma once
#include <vector>

namespace aco 
{

// Hamiltonian cycle: node.size()==N, car.size()==N, node[0]==0
struct Solution 
{
    std::vector<int> node;
    std::vector<int> car;
};

struct EvaluatedSolution 
{
    Solution sol;
    double cost = 0.0;
};

} // namespace aco
