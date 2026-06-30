#pragma once
// NOVO: Cache wrapper za DP evaluateCost/reassignCars.
// - koristi IFixedTourCarAssigner (ne dira DP implementaciju)
// - canonical key: nodes (node<a href="" class="citation-link" target="_blank" style="vertical-align: super; font-size: 0.8em; margin-left: 3px;">[0]</a>==0) + optional reverse ako je instanca simetrična.

#include <interfaces.hpp> // ako nemaš ovaj file, zamijeni sa stvarnim headerom gdje su IFixedTourCarAssigner/Solution/...
#include <parser.hpp>    // ili stvarni header za cars_tsplib::Instance

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
    // Maks broj zapisa u cacheu (FIFO eviction).
    // 0 => cache isključen (wrapper samo forwarda pozive).
    std::size_t capacity = 5000;

    // Ako true: cache pamti i carPerEdge za reassignCars hit.
    bool storeCars = true;

    // Ako true: ako tagovi simetrije nisu prisutni, probaj auto-detect samplingom.
    bool autoDetectSymmetryWhenUnknown = true;

    // Koliko random uzoraka provjeriti za simetriju (ako je unknown).
    // 0 => preskoči auto-detect.
    int symmetryDetectSamples = 2000;

    // Epsilon za usporedbu simetrije (floating).
    double symmetryEps = 1e-12;

    // Ako true: cache pretpostavlja da je nodes<a href="" class="citation-link" target="_blank" style="vertical-align: super; font-size: 0.8em; margin-left: 3px;">[0]</a>==0 i ne radi rotaciju.
    // (kod tebe to vrijedi svugdje; drži true radi brzine)
    bool assumeStartIsZero = true;

    // Basic thread-safety (mutex oko map-a). Ako znaš da je single-thread per Colony, može false.
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

    // (opcionalno) statistika
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

    // hash -> list of indices (za kolizije hash-a)
    mutable std::unordered_map<std::uint64_t, std::vector<std::size_t>> buckets_;

    // scratch za key (izbjegava re-alokacije u hot pathu)
    mutable std::vector<int> scratchKey_;

    // stats
    mutable Stats stats_;

    // thread-safety (opcionalno)
    mutable std::mutex mx_;
};

} // namespace aco
