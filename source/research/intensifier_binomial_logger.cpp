#include <research/intensifier_binomial_logger.hpp>
#include <fstream>
#include <set>
#include <string>

namespace aco
{

void IntensifierBinomialLogger::onActivationStart(int runIdx, int activationIndex, double costAtStart)
{
    std::lock_guard<std::mutex> lock(mtx_);
    records_.push_back({ runIdx, activationIndex, costAtStart, costAtStart, false, 0 });
    // improved = false since this is the start; costAtEnd = costAtStart for consistency
}

void IntensifierBinomialLogger::resetActivationDepth(int runIdx)
{
    std::lock_guard<std::mutex> lock(mtx_);
    activationDepthByRun_[runIdx] = 0;
}

void IntensifierBinomialLogger::noteGlobalStagnationResetDuringIntensifier(int runIdx)
{
    std::lock_guard<std::mutex> lock(mtx_);
    ++activationDepthByRun_[runIdx];
}

void IntensifierBinomialLogger::onDeactivation(int runIdx, int activationIndex,
                                               double costAtStart, double costAtEnd, bool improved)
{
    std::lock_guard<std::mutex> lock(mtx_);
    int depth = 0;
    if (const auto it = activationDepthByRun_.find(runIdx); it != activationDepthByRun_.end()) {
        depth = it->second;
        activationDepthByRun_.erase(it);
    }
    records_.push_back({ runIdx, activationIndex, costAtStart, costAtEnd, improved, depth });
}

void IntensifierBinomialLogger::writeCsv(const std::string& filename, char sep) const
{
    std::vector<IntensifierActivationRecord> copy;
    {
        std::lock_guard<std::mutex> lock(mtx_);
        copy = records_;
    }

    std::ofstream f(filename);
    if (!f)
        return;

    f << "run" << sep << "activation_index" << sep << "cost_start" << sep << "cost_end" << sep
      << "improved" << sep << "depth\n";
    for (const auto& r : copy)
    {
        f << r.runIndex << sep << r.activationIndex << sep
          << r.costAtStart << sep << r.costAtEnd << sep
          << (r.improved ? 1 : 0) << sep << r.depth << "\n";
    }
}

size_t IntensifierBinomialLogger::numRunsWithData() const
{
    std::lock_guard<std::mutex> lock(mtx_);
    std::set<int> runs;
    for (const auto& r : records_)
        runs.insert(r.runIndex);
    return runs.size();
}

std::vector<IntensifierActivationRecord> IntensifierBinomialLogger::getRecords() const
{
    std::lock_guard<std::mutex> lock(mtx_);
    return records_;
}

} // namespace aco
