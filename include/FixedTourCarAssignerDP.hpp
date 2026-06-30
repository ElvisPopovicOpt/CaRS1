#pragma once

#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cars_tsplib { struct Instance; }

namespace localSearch 
{

// Rezultat DP-a: minimalni cost i optimalni car po edge-u (size N).
struct CarAssignmentDPResult
{
    double cost = std::numeric_limits<double>::infinity();
    std::vector<int> carPerEdge; // size = N
};

struct CarAssignmentDPOptions
{
    // Ako true: provjeri da je tour permutacija [0..N-1] i da je node<a href="" class="citation-link" target="_blank" style="vertical-align: super; font-size: 0.8em; margin-left: 3px;">[0]</a>==0.
    // Korisno dok razvijaš 2-opt/LK (hvata tihe bugove).
    bool validateTour = false;

    // Heurističko ograničenje duljine segmenta (radi brzine).
    // -1 => bez ograničenja (točno optimalno).
    int maxSegmentLen = -1;

    // ROBUST fallback:
    // Ako je maxSegmentLen > 0 i DP vrati INF, automatski retry s većim L, pa na kraju s -1.
    bool fallbackToUnlimitedOnInf = true;

    // Koliko puta pokušati povećati L prije finalnog -1 (ako je enabled).
    // Primjer: maxSegmentLen=20, growth=2, retries=3 => 20,40,80 pa -1.
    int fallbackMaxRetries = 3;

    // Faktor rasta segmenta pri retryu (>=2 preporučeno).
    int fallbackGrowthFactor = 2;
};

// Da se ne alocira vektor tijekom izvršavanja
// N - broj cvorova, C broj automobila
struct Scratch 
{
    int N = 0;
    int C = 0;
    int M = 0; // 1<<C

    // pref[c*(N+1) + t]
    std::vector<double> pref;

    // dp[t*M + mask]
    std::vector<double> dp;

    // materijalizirana tura (N+1, uključuje closure tour[N]=tour<a href="" class="citation-link" target="_blank" style="vertical-align: super; font-size: 0.8em; margin-left: 3px;">[0]</a>)
    std::vector<int> tour;

    // NOVO: za validateTour bez alokacija/clearanja
    std::vector<int> seenStamp;
    int stamp = 1;

    void ensure(int N_, int C_);
};




// DP optimizator koji za fiksni Hamiltonov ciklus čvorova (node<a href="" class="citation-link" target="_blank" style="vertical-align: super; font-size: 0.8em; margin-left: 3px;">[0]</a>==0)
// pronalazi optimalnu segmentaciju po autima pod CaRS "switch-return" semantikom.
//
// Semantika:
// - segment [s, t) koristi auto c na edge-ovima s..t-1
// - cijena segmenta je:
//      sum travelCost(c, v[k], v[k+1]) za k=s..t-1
//    + returnCost(c, v[t], v[s])   // zatvaranje segmenta (switch/final)
//
// Time DP prirodno modelira i finalni return zadnjeg auta jer je v[N]=v[0].
// update - robust fallback
// Semantika:
// - segment [s, t) koristi auto c na edge-ovima s..t-1
// - cijena segmenta je:
//      sum travelCost(c, v[k], v[k+1]) za k=s..t-1
//    + returnCost(c, v[t], v[s])   // zatvaranje segmenta (switch/final)
//
// Time DP prirodno modelira i finalni return zadnjeg auta jer je v[N]=v[0].

class FixedTourCarAssignerDP
{
public:
    explicit FixedTourCarAssignerDP(std::shared_ptr<const cars_tsplib::Instance> inst,
                                    CarAssignmentDPOptions opt = {});

    CarAssignmentDPResult optimize(const std::vector<int>& nodes) const;

    double reassignCars(const std::vector<int>& nodes,
                        std::vector<int>& carPerEdgeOut) const;

    double evaluateCost(const std::vector<int>& nodes) const;

    // cost-only evaluacija ture preko "view"-a.
    // getNode(idx) mora vratiti čvor na poziciji idx (0..N-1), a ciklus se zatvara interno.
    double evaluateCostView(const std::function<int(int)>& getNode) const;

     // NOVO: koristi scratch (bez realokacije)
    double evaluateCostViewScratch(const std::function<int(int)>& getNode, Scratch& scratch) const;

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    CarAssignmentDPOptions opt_;

    void validateTourOrThrow_(const std::vector<int>& nodes) const;

    // Internal: single run bez fallbacka, s eksplicitnim maxSegmentLen parametrom.
    CarAssignmentDPResult optimizeWithMaxSegLen_(const std::vector<int>& nodes,
                                                int maxSegmentLen) const;
    
    double evaluateCostViewWithMaxSegLen_(const std::function<int(int)>& getNode,
                                         int maxSegmentLen) const;

    // postojeće (ali sada implementiraj varijantu koja prima scratch)
    double evaluateCostViewWithMaxSegLenScratch_(const std::function<int(int)>& getNode,
                                                int maxSegmentLen,
                                                Scratch& scratch) const;

                                         // Internal: fallback wrapper za optimize() / evaluateCostView()
    CarAssignmentDPResult optimizeWithFallback_(const std::vector<int>& nodes) const;

    double evaluateCostViewWithFallback_(const std::function<int(int)>& getNode) const;

    int normalizeMaxSegLen_(int L) const;

    // Opcija 1 (najjednostavnije): scratch po DP instanci
    // mutable jer evaluateCostView je const
    mutable Scratch scratch_;
};

} // namespace localSearch
