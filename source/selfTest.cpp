#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <parser.hpp>
#include <FixedTourCarAssignerDP.hpp>   // localSearch::FixedTourCarAssignerDP
#include <cost_model.hpp>            // aco::CostModelCars
#include <solution.hpp>            // aco::Solution (prilagodi imena include-a)
#include <selfTest.hpp>

namespace selftest 
{

static bool approxEqual(double a, double b, double eps)
{
    // relativno + apsolutno (robustno za male i velike vrijednosti)
    const double diff = std::fabs(a - b);
    const double scale = std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
    return diff <= eps * scale;
}

void runDpConsistencyTest(std::shared_ptr<const cars_tsplib::Instance> inst,
                          std::uint64_t seed,
                          int toursToTest,
                          double eps)
{
    if (!inst) throw std::runtime_error("runDpConsistencyTest: inst is null");

    const int N = inst->n();
    const int C = inst->cars();
    if (N <= 1) throw std::runtime_error("runDpConsistencyTest: N must be > 1");
    if (C <= 0) throw std::runtime_error("runDpConsistencyTest: C must be > 0");
    if (!inst->hasReturnCosts()) throw std::runtime_error("runDpConsistencyTest: instance has no return costs");

    aco::CostModelCars evaluator(inst);

    localSearch::CarAssignmentDPOptions opt;
    opt.validateTour = true;    // hvata bugove (duplikati, out-of-range, node<a href="" class="citation-link" target="_blank" style="vertical-align: super; font-size: 0.8em; margin-left: 3px;">[0]</a>!=0)
    opt.maxSegmentLen = -1;     // exact optimal

    localSearch::FixedTourCarAssignerDP dp(inst, opt);

    std::mt19937_64 rng(seed);

    std::vector<int> nodes(static_cast<std::size_t>(N));
    nodes[0] = 0;

    // bazna permutacija 1..N-1
    std::vector<int> perm(static_cast<std::size_t>(N - 1));
    std::iota(perm.begin(), perm.end(), 1);

    for (int t = 0; t < toursToTest; ++t) {
        std::shuffle(perm.begin(), perm.end(), rng);
        for (int i = 1; i < N; ++i) nodes[static_cast<std::size_t>(i)] = perm[static_cast<std::size_t>(i - 1)];

        // DP optimalni car assignment
        auto dpRes = dp.optimize(nodes);

        aco::Solution sol;
        sol.node = nodes;
        sol.car  = std::move(dpRes.carPerEdge);

        const double costEval = evaluator.evaluate(sol);

        if (!approxEqual(dpRes.cost, costEval, eps)) {
            std::cerr << "DP consistency test FAILED on tour " << t << "\n";
            std::cerr << "DP cost      = " << dpRes.cost << "\n";
            std::cerr << "Evaluator    = " << costEval << "\n";
            std::cerr << "Abs diff     = " << std::fabs(dpRes.cost - costEval) << "\n";
            throw std::runtime_error("runDpConsistencyTest: DP cost != evaluator cost");
        }
    }

    std::cout << "[SelfTest] DP consistency OK (" << toursToTest << " random tours), eps=" << eps << "\n";
}

} // namespace selftest
