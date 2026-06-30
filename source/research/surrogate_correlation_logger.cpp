#include <research/surrogate_correlation_logger.hpp>
#include <fstream>
#include <string>

namespace aco
{

void SurrogateCorrelationLogger::log(int runIdx, int iter, int antIdx, double surrogate, double dpCost)
{
    std::lock_guard<std::mutex> lock(mtx_);
    records_.push_back({ runIdx, iter, antIdx, surrogate, dpCost });
}

void SurrogateCorrelationLogger::writeCsv(const std::string& filename, char sep) const
{
    std::vector<SurrogateDpRecord> copy;
    {
        std::lock_guard<std::mutex> lock(mtx_);
        copy = records_;
    }

    std::ofstream f(filename);
    if (!f)
        return;

    f << "run" << sep << "iteration" << sep << "ant_index" << sep << "surrogate" << sep << "dp_cost\n";
    for (const auto& r : copy)
        f << r.runIndex << sep << r.iteration << sep << r.antIndex << sep << r.surrogate << sep << r.dpCost << "\n";
}

size_t SurrogateCorrelationLogger::numRecords() const
{
    std::lock_guard<std::mutex> lock(mtx_);
    return records_.size();
}

std::vector<SurrogateDpRecord> SurrogateCorrelationLogger::getRecords() const
{
    std::lock_guard<std::mutex> lock(mtx_);
    return records_;
}

} // namespace aco
