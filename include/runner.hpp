#pragma once
#include <memory>
#include <vector>
#include <thread>
#include <functional>
#include <mutex>
#include <atomic>
#include <colony.hpp>
#include <rng.hpp>

// forward declaration
namespace cars_tsplib { struct Instance; }
namespace aco_cli { struct ParamsData; }

namespace aco 
{

struct RunnerConfig 
{
    int threads = 0; // 0 => hardware_concurrency
};

class Runner 
{
public:
    explicit Runner(RunnerConfig cfg = {});

    // ColonyFactory: kreira fully-wired Colony za runIndex/seed.
    using ColonyFactory = std::function<std::unique_ptr<Colony>(int runIndex, uint64_t seed)>;

    std::vector<RunResult> runAll(const aco_cli::ParamsData& params, uint64_t baseSeed, ColonyFactory factory);

private:
    RunnerConfig cfg_;
};

} // namespace aco
