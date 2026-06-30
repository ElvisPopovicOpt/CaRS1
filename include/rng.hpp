#pragma once
#include <cstdint>
#include <random>
#include <algorithm>
#include <type_traits>

namespace aco 
{

// splitmix64: odličan za deterministički "seed mix"
inline uint64_t splitmix64(uint64_t& x) 
{
    uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

inline uint64_t mixSeed(uint64_t baseSeed, uint64_t runIndex) 
{
    uint64_t x = baseSeed ^ (runIndex + 0x9e3779b97f4a7c15ULL);
    // nekoliko rundi za stabilan mix
    uint64_t a = splitmix64(x);
    uint64_t b = splitmix64(x);
    return a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2));
}

// ----- SplitMix64 (deterministički seed-mixer) -----
inline uint64_t splitmix64_step(uint64_t& x)
{
    uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

// mix baseSeed + (a,b,c) u novi seed (stabilno i dobro raspršeno)
inline uint64_t mixSeed4(uint64_t base, uint64_t a, uint64_t b, uint64_t c)
{
    uint64_t x = base ^ 0x9e3779b97f4a7c15ULL;

    x ^= a + 0x9e3779b97f4a7c15ULL + (x << 6) + (x >> 2);
    uint64_t r1 = splitmix64_step(x);

    x ^= b + 0x9e3779b97f4a7c15ULL + (x << 6) + (x >> 2);
    uint64_t r2 = splitmix64_step(x);

    x ^= c + 0x9e3779b97f4a7c15ULL + (x << 6) + (x >> 2);
    uint64_t r3 = splitmix64_step(x);

    return r1 ^ (r2 + 0x9e3779b97f4a7c15ULL + (r1 << 6) + (r1 >> 2)) ^ (r3 << 1);
}

// najčešći slučaj: baseSeed + runIndex
inline uint64_t mixSeedRun(uint64_t baseSeed, uint64_t runIndex)
{
    // b,c = 0 => deterministički
    return mixSeed4(baseSeed, 0xA11CE000ULL, runIndex, 0);
}

// ----- Komponentni tagovi (odvojeni RNG tokovi) -----
static constexpr uint64_t TAG_ANTS   = 0xA11CE001ULL;
static constexpr uint64_t TAG_KICK   = 0xA11CE002ULL;
static constexpr uint64_t TAG_LKLITE = 0xA11CE003ULL;
static constexpr uint64_t TAG_THREEOPT = 0xA11CE004ULL;

// ----- RNG wrapper (mt19937_64) -----
class Rng
{
public:
    Rng() : eng_(0) {}
    explicit Rng(uint64_t seed) : eng_(static_cast<std::mt19937_64::result_type>(seed)) {}

    void reseed(uint64_t seed)
    {
        eng_.seed(static_cast<std::mt19937_64::result_type>(seed));
    }

    // [0,1)
    double uniform01()
    {
        return std::generate_canonical<double, 64>(eng_);
    }

    // int in [0, hiExclusive)
    int uniformInt(int hiExclusive)
    {
        std::uniform_int_distribution<int> dist(0, hiExclusive - 1);
        return dist(eng_);
    }

    // int in [loInclusive, hiInclusive]
    int uniformInt(int loInclusive, int hiInclusive)
    {
        return loInclusive + uniformInt(hiInclusive - loInclusive + 1);
    }

    template <class It>
    void shuffle(It b, It e)
    {
        std::shuffle(b, e, eng_);
    }

    std::mt19937_64& engine() { return eng_; }
    const std::mt19937_64& engine() const { return eng_; }

private:
    std::mt19937_64 eng_;
};

} // namespace aco
