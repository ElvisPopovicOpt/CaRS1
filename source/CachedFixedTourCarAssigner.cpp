#include "CachedFixedTourCarAssigner.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

namespace aco
{

static bool isFinite(double x) { return std::isfinite(x); }

CachedFixedTourCarAssigner::CachedFixedTourCarAssigner(
    std::shared_ptr<const cars_tsplib::Instance> inst,
    std::shared_ptr<const IFixedTourCarAssigner> inner,
    DPCacheOptions opt)
    : inst_(std::move(inst))
    , inner_(std::move(inner))
    , opt_(opt)
{
    if (!inst_)  throw std::runtime_error("CachedFixedTourCarAssigner: inst is null");
    if (!inner_) throw std::runtime_error("CachedFixedTourCarAssigner: inner is null");

    allowReverse_ = computeAllowReverse_(*inst_);
    stats_.allowReverse = allowReverse_;

    if (opt_.capacity > 0)
    {
        entries_.resize(opt_.capacity);
        buckets_.reserve(opt_.capacity * 2);
        scratchKey_.reserve(inst_->n());
    }
}

void CachedFixedTourCarAssigner::clear()
{
    if (opt_.threadSafe)
    {
        std::lock_guard<std::mutex> lk(mx_);
        buckets_.clear();
        for (auto& e : entries_) e = Entry{};
        cursor_ = 0;
        stats_ = Stats{};
        stats_.allowReverse = allowReverse_;
        return;
    }

    buckets_.clear();
    for (auto& e : entries_) e = Entry{};
    cursor_ = 0;
    stats_ = Stats{};
    stats_.allowReverse = allowReverse_;
}

CachedFixedTourCarAssigner::Stats CachedFixedTourCarAssigner::stats() const
{
    if (opt_.threadSafe)
    {
        std::lock_guard<std::mutex> lk(mx_);
        return stats_;
    }
    return stats_;
}

bool CachedFixedTourCarAssigner::computeAllowReverse_(const cars_tsplib::Instance& inst) const
{
    // Sigurno reverse canonicaliziranje samo ako je instanca stvarno simetrična.

    // 1) Travel symmetry mora biti eksplicitno poznata i true, ili auto-detect mora potvrditi.
    bool travelSymKnownTrue = inst.edgeWeightIsSymmetric.has_value() && *inst.edgeWeightIsSymmetric;
    bool travelSymKnownFalse = inst.edgeWeightIsSymmetric.has_value() && !*inst.edgeWeightIsSymmetric;

    if (travelSymKnownFalse) return false;

    // 2) Return: ako postoji i tag kaže ASYMMETRIC => nema reverse.
    if (inst.hasReturnCosts())
    {
        if (inst.returnRateIsAsymmetric.has_value() && *inst.returnRateIsAsymmetric)
            return false;
    }

    // Ako su tagovi kompletni i podržavaju simetriju:
    if (travelSymKnownTrue)
    {
        if (!inst.hasReturnCosts()) return true;

        // Return costs prisutni: reverse je siguran samo ako imamo eksplicitno "nije asimetrično"
        // (ili ako parser ima nekakav "symmetric" tag; ti imaš samo returnRateIsAsymmetric optional).
        if (inst.returnRateIsAsymmetric.has_value() && !*inst.returnRateIsAsymmetric)
            return true;

        // Unknown return symmetry -> konzervativno false (osim auto-detect).
    }

    // Unknown: opcionalno auto-detect
    if (opt_.autoDetectSymmetryWhenUnknown && opt_.symmetryDetectSamples > 0)
        return detectSymmetryBySampling_(inst);

    return false;
}

bool CachedFixedTourCarAssigner::detectSymmetryBySampling_(const cars_tsplib::Instance& inst) const
{
    const int N = inst.n();
    const int C = inst.cars();
    if (N <= 1 || C <= 0) return true;

    // Deterministički RNG (ne želimo ovisiti o global RNG-u)
    std::uint64_t seed = 0x9e3779b97f4a7c15ULL;
    seed ^= (std::uint64_t)N + 0xBF58476D1CE4E5B9ULL;
    seed ^= (std::uint64_t)C + 0x94D049BB133111EBULL;

    auto nextU64 = [&]() {
        seed = splitmix64_(seed);
        return seed;
    };

    auto rndInt = [&](int hi) -> int {
        return (int)(nextU64() % (std::uint64_t)std::max(1, hi));
    };

    const int samples = std::max(1, opt_.symmetryDetectSamples);
    const double eps = opt_.symmetryEps;

    for (int t = 0; t < samples; ++t)
    {
        const int i = rndInt(N);
        const int j = rndInt(N);
        const int c = rndInt(C);

        const double a = inst.travelCost(c, i, j);
        const double b = inst.travelCost(c, j, i);
        if (std::fabs(a - b) > eps) return false;

        if (inst.hasReturnCosts())
        {
            const double r1 = inst.returnCost(c, i, j);
            const double r2 = inst.returnCost(c, j, i);
            if (std::fabs(r1 - r2) > eps) return false;
        }
    }

    return true;
}

void CachedFixedTourCarAssigner::buildCanonicalKey_(const std::vector<int>& nodes,
                                                   std::vector<int>& out) const
{
    out = nodes;

    if (!allowReverse_) return;
    const int N = (int)nodes.size();
    if (N <= 1) return;

    // nodes<a href="" class="citation-link" target="_blank" style="vertical-align: super; font-size: 0.8em; margin-left: 3px;">[0]</a>==0 pretpostavka: reverse forma je [0, nodes[N-1], nodes[N-2], ..., nodes<a href="" class="citation-link" target="_blank" style="vertical-align: super; font-size: 0.8em; margin-left: 3px;">[1]</a>]
    // Odluči leksikografski bez alokacije, pa tek onda eventualno napravi reverse.
    bool useRev = false;
    for (int i = 0; i < N; ++i)
    {
        const int f = nodes[(std::size_t)i];
        const int r = (i == 0) ? 0 : nodes[(std::size_t)(N - i)];
        if (r < f) { useRev = true; break; }
        if (r > f) { useRev = false; break; }
    }

    if (!useRev) return;

    out.resize((std::size_t)N);
    out[0] = 0;
    for (int i = 1; i < N; ++i)
        out[(std::size_t)i] = nodes[(std::size_t)(N - i)];
}

std::uint64_t CachedFixedTourCarAssigner::splitmix64_(std::uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

std::uint64_t CachedFixedTourCarAssigner::hashKey_(const std::vector<int>& key) const
{
    std::uint64_t h = 0xD6E8FEB86659FD93ULL;
    for (int v : key)
    {
        std::uint64_t x = (std::uint64_t)(std::uint32_t)v;
        h ^= splitmix64_(x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
    }
    // dodatno ubaci duljinu
    h ^= splitmix64_((std::uint64_t)key.size());
    return h;
}

const CachedFixedTourCarAssigner::Entry*
CachedFixedTourCarAssigner::findEntry_(const std::vector<int>& key, std::uint64_t h) const
{
    auto it = buckets_.find(h);
    if (it == buckets_.end()) return nullptr;

    const auto& idxs = it->second;
    for (std::size_t idx : idxs)
    {
        const Entry& e = entries_[idx];
        if (!e.valid) continue;
        if (e.h != h) continue;
        if (e.key == key) return &e;
    }
    return nullptr;
}

CachedFixedTourCarAssigner::Entry*
CachedFixedTourCarAssigner::findEntryMut_(const std::vector<int>& key, std::uint64_t h) const
{
    auto it = buckets_.find(h);
    if (it == buckets_.end()) return nullptr;

    auto& idxs = it->second;
    for (std::size_t idx : idxs)
    {
        Entry& e = entries_[idx];
        if (!e.valid) continue;
        if (e.h != h) continue;
        if (e.key == key) return &e;
    }
    return nullptr;
}

void CachedFixedTourCarAssigner::eraseIndexFromBucket_(std::uint64_t h, std::size_t index) const
{
    auto it = buckets_.find(h);
    if (it == buckets_.end()) return;

    auto& v = it->second;
    for (std::size_t k = 0; k < v.size(); ++k)
    {
        if (v[k] == index)
        {
            v[k] = v.back();
            v.pop_back();
            break;
        }
    }
    if (v.empty()) buckets_.erase(it);
}

void CachedFixedTourCarAssigner::insertOrUpdate_(const std::vector<int>& key,
                                                std::uint64_t h,
                                                double cost,
                                                const std::vector<int>* carsOrNull) const
{
    if (opt_.capacity == 0) return;

    // update ako postoji
    if (Entry* ex = findEntryMut_(key, h))
    {
        ex->cost = cost;
        if (carsOrNull && opt_.storeCars)
        {
            ex->cars = *carsOrNull;
            ex->hasCars = true;
        }
        return;
    }

    // insert FIFO (overwrite na cursor_)
    const std::size_t index = cursor_;
    cursor_ = (cursor_ + 1) % entries_.size();

    Entry& slot = entries_[index];

    // izbaci stari zapis iz bucket-a
    if (slot.valid)
        eraseIndexFromBucket_(slot.h, index);

    slot.valid = true;
    slot.h = h;
    slot.key = key;
    slot.cost = cost;

    if (carsOrNull && opt_.storeCars)
    {
        slot.cars = *carsOrNull;
        slot.hasCars = true;
    }
    else
    {
        slot.cars.clear();
        slot.hasCars = false;
    }

    buckets_[h].push_back(index);
    ++stats_.inserts;
}

double CachedFixedTourCarAssigner::evaluateCost(const std::vector<int>& nodes) const
{
    if (opt_.threadSafe)
    {
        std::lock_guard<std::mutex> lk(mx_);
        // padamo na “unsafe” implementaciju unutar lock-a
    }

    if (opt_.capacity == 0)
        return inner_->evaluateCost(nodes);

    if (opt_.assumeStartIsZero && !nodes.empty() && nodes[0] != 0)
        throw std::runtime_error("CachedFixedTourCarAssigner: nodes[0] must be 0");

    buildCanonicalKey_(nodes, scratchKey_);
    const std::uint64_t h = hashKey_(scratchKey_);

    if (const Entry* e = findEntry_(scratchKey_, h))
    {
        ++stats_.evalHits;
        return e->cost;
    }

    ++stats_.evalMiss;
    const double c = inner_->evaluateCost(nodes);

    if (isFinite(c))
        insertOrUpdate_(scratchKey_, h, c, /*carsOrNull*/ nullptr);

    return c;
}

double CachedFixedTourCarAssigner::reassignCars(const std::vector<int>& nodes,
                                               std::vector<int>& carPerEdgeOut) const
{
    if (opt_.threadSafe)
    {
        std::lock_guard<std::mutex> lk(mx_);
        // padamo na “unsafe” implementaciju unutar lock-a
    }

    if (opt_.capacity == 0)
        return inner_->reassignCars(nodes, carPerEdgeOut);

    if (opt_.assumeStartIsZero && !nodes.empty() && nodes[0] != 0)
        throw std::runtime_error("CachedFixedTourCarAssigner: nodes[0] must be 0");

    buildCanonicalKey_(nodes, scratchKey_);
    const std::uint64_t h = hashKey_(scratchKey_);

    if (const Entry* e = findEntry_(scratchKey_, h))
    {
        if (e->hasCars && opt_.storeCars)
        {
            ++stats_.reasHits;
            carPerEdgeOut = e->cars;
            return e->cost;
        }
        // ima cost, ali nema cars: i dalje je hit za cost, ali moramo izračunati cars
        // (ne brojimo kao reasHits jer nije full hit)
    }

    ++stats_.reasMiss;
    const double c = inner_->reassignCars(nodes, carPerEdgeOut);

    if (isFinite(c))
        insertOrUpdate_(scratchKey_, h, c, opt_.storeCars ? &carPerEdgeOut : nullptr);

    return c;
}

} // namespace aco
