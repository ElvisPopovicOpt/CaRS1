#pragma once
#include <memory>
#include <vector>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <algorithm>

#include "interfaces.hpp"          // aco::ILocalSearch, aco::Solution
#include "rng.hpp"                 // aco::Rng
#include "candidate_list.hpp"  // CandidateListCache
#include "FixedTourCarAssignerDP.hpp"

namespace localSearch
{

struct LkLiteOptions
{
    // Koliko pokušaja (svaki pokušaj = 1x dp->reassignCars).
    int attempts = 6;

    // Duljina lanca 2-opt poteza (surrogate vođeno).
    int chainLen = 12;

    // Koliko random i-pozicija po koraku (svaki puta gleda topKPerI kandidata).
    int triesPerStep = 40;

    // Koliko kandidata (iz candidate liste) uzeti po anchoru za jedan i.
    int topKPerI = 8;

    // Ako true: radi potez samo ako surrogate delta < -eps
    bool requireNegativeSurrogate = true;

    double eps = 1e-12;

    // Candidate list (obavezno)
    std::shared_ptr<const aco::CandidateListCache> cand = nullptr;

    // Opcionalno: nakon uspješnog LK-lite verify, pusti polish LS (npr. 2opt+reloc+2opt)
    std::shared_ptr<const aco::ILocalSearch> polish = nullptr;
};

// LK-lite: radi chain 2-opt poteza vođenih jeftinim surrogate-om (minTravel),
// pa 1x DP verify (reassignCars). Aktivira se na stagnaciju.
class LkLiteLocalSearch final : public aco::ILocalSearch
{
public:
    LkLiteLocalSearch(std::shared_ptr<const cars_tsplib::Instance> inst,
                      LkLiteOptions opt,
                      std::shared_ptr<FixedTourCarAssignerDP> dp);

    // Reset RNG stanja (pozovi na početku stagnation-eventa s deriviranim seedom)
    void resetSeed(uint64_t seed) const override;

    void improve(aco::Solution& s, double& cost) const override;

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    LkLiteOptions opt_;
    std::shared_ptr<FixedTourCarAssignerDP> dp_;

    mutable aco::Rng rng_;

    static bool isValid2OptMove(int i, int k, int N);
    static void applyReverseAndUpdatePos(std::vector<int>& tour, std::vector<int>& pos, int i, int k);

    double minTravel_(int i, int j) const;
    double surrogate2OptDelta_(const std::vector<int>& tour, int i, int k) const;
};

} // namespace localSearch
