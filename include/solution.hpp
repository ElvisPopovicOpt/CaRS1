#pragma once
#include <vector>

namespace aco 
{

// Hamiltonov ciklus: node.size()==N, car.size()==N, node<a href="" class="citation-link" target="_blank" style="vertical-align: super; font-size: 0.8em; margin-left: 3px;">[0]</a>==0
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
