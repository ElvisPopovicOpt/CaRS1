#pragma once
// Cache wrapper around a DP-based IFixedTourCarAssigner's evaluateCost/reassignCars.
// Canonical key: nodes[0] == 0, plus the reversed form if the instance is symmetric.

#include <interfaces.hpp>
#include <parser.hpp>

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>
#include <mutex>
#include <optional>

namespace aco
{

struct DPCacheOptions
{
    // Max number of entries (FIFO eviction). 0 disables the cache (wrapper just forwards calls).
    std::size_t capacity = 5000;

    // If true, also cache carPerEdge so reassignCars can hit.
    bool storeCars = true;

    // If true, auto-detect symmetry by sampling when symmetry tags are missing.
    bool autoDetectSymmetryWhenUnknown = true;

    // Number of random samples to check for symmetry when unknown. 0 skips auto-detect.
    int symmetryDetectSamples = 2000;

    // Epsilon for floating-point symmetry comparison.
    double symmetryEps = 1e-12;

    // If true, the cache assumes nodes[0] == 0 and skips rotation (holds everywhere; keep true for speed).
    bool assumeStartIsZero = true;

    // Basic thread-safety (mutex around the map). Safe to disable if usage is single-threaded per Colony.
    bool threadSafe = false;
};

class CachedFixedTourCarAssigner final : public IFixedTourCarAssigner
{
public:
    CachedFixedTourCarAssigner(std::shared_ptr<const cars_tsplib::Instance> inst,
                               std::shared_ptr<const IFixedTourCarAssigner> inner,
                               DPCacheOptions opt = {});

    void clear();

    // IFixedTourCarAssigner
    double reassignCars(const std::vector<int>& nodes,
                        std::vector<int>& carPerEdgeOut) const override;

    double evaluateCost(const std::vector<int>& nodes) const override;

    // Optional stats
    struct Stats
    {
        std::uint64_t evalHits = 0;
        std::uint64_t evalMiss = 0;
        std::uint64_t reasHits = 0;
        std::uint64_t reasMiss = 0;
        std::uint64_t inserts  = 0;
        bool allowReverse = false;
    };
    Stats stats() const;

private:
    struct Entry
    {
        bool valid = false;
        std::uint64_t h = 0;
        std::vector<int> key;      // canonical nodes
        double cost = 0.0;

        bool hasCars = false;
        std::vector<int> cars;     // carPerEdge
    };

private:
    bool computeAllowReverse_(const cars_tsplib::Instance& inst) const;
    bool detectSymmetryBySampling_(const cars_tsplib::Instance& inst) const;

    void buildCanonicalKey_(const std::vector<int>& nodes,
                            std::vector<int>& out) const;

    std::uint64_t hashKey_(const std::vector<int>& key) const;
    static std::uint64_t splitmix64_(std::uint64_t x);

    // Lookup/insert
    const Entry* findEntry_(const std::vector<int>& key, std::uint64_t h) const;
    Entry*       findEntryMut_(const std::vector<int>& key, std::uint64_t h) const;

    void insertOrUpdate_(const std::vector<int>& key,
                         std::uint64_t h,
                         double cost,
                         const std::vector<int>* carsOrNull) const;

    void eraseIndexFromBucket_(std::uint64_t h, std::size_t index) const;

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    std::shared_ptr<const IFixedTourCarAssigner> inner_;
    DPCacheOptions opt_;

    bool allowReverse_ = false;

    // FIFO storage
    mutable std::vector<Entry> entries_;
    mutable std::size_t cursor_ = 0;

    // hash -> list of indices (for hash collisions)
    mutable std::unordered_map<std::uint64_t, std::vector<std::size_t>> buckets_;

    // Scratch key buffer (avoids reallocations in the hot path)
    mutable std::vector<int> scratchKey_;

    // stats
    mutable Stats stats_;

    // Optional thread-safety
    mutable std::mutex mx_;
};

} // namespace aco
