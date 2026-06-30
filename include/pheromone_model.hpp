#pragma once
#include <memory>
#include <vector>
#include <cstdint>
#include <interfaces.hpp>

namespace cars_tsplib { struct Instance; }
namespace aco_cli { struct ParamsData; }

namespace aco 
{

struct ArchiveOptions
{
    bool enabled = true;
    int size = 30;               // typical: 20-40 for N~300
    int depositPerIter = 1;      // 1 or 2
    double depositWeight = 0.25; // small memory signal
    bool hashIncludesCars = false; // recommended: false (node-only)
    bool canonicalizeReverseIfSymmetric = true; // canonicalize for symmetric instances
};

struct ExplorationOptions
{
    bool enabled = false;
    int iters = 20;                 // typical: 10-30
    double rhoMultiplier = 2.0;     // temporarily boosts evaporation (e.g. 1.5-3.0)
};

class PheromoneModelMMASMove final : public IPheromoneModel
{
public:
    PheromoneModelMMASMove(std::shared_ptr<const cars_tsplib::Instance> inst,
                           const aco_cli::ParamsData& params,
                           ExplorationOptions exploration = {},
                           ArchiveOptions archive = {});

    void reset() override;
    double tauMove(int car, int i, int j) const override;
    double tauReturn(int car, int from, int to) const override;
    bool hasReturnPheromones() const override { return hasReturn_; }

    void onIterationEnd(int iterationIndex,
                        const EvaluatedSolution& iterBest,
                        const EvaluatedSolution& globalBest,
                        bool globalImproved) override;

    void onIterationEndRanked(int iterationIndex,
                          const std::vector<EvaluatedSolution>& ranked,
                          const EvaluatedSolution& globalBest,
                          bool globalImproved) override;
    
    // Check if smoothing was applied in the last iteration
    bool wasSmoothingApplied() const override { return lastSmoothingApplied_; }
    
    // Get stagnation counter (noImprove_) for smoothing calculation
    int getStagnationCounter() const override { return noImprove_; }
    
    // Selectively reduces pheromones on a given tour (smoothing)
    void reducePheromonesOnTour(const Solution& solution, double gamma) override;

    // Adds a solution to the archive
    void addToArchive(const EvaluatedSolution& solution) override;

    // Returns the best archived cost
    double getBestArchivedCost() const override;

    // Returns the best archived solution
    bool getBestArchivedSolution(EvaluatedSolution& out) const override;

private:
    // Archive memory (global best pool)
    ArchiveOptions archiveOpt_;

    struct ArchiveEntry
    {
        EvaluatedSolution s;
        uint64_t h = 0;
    };
    std::vector<ArchiveEntry> archive_;
    size_t archiveCursor_ = 0; // round-robin deposit cursor

    bool treatAsAsymmetric_() const;   // mirrors Colony::preferThreeOptPolish
    bool treatAsSymmetric_() const { return !treatAsAsymmetric_(); }

    void buildCanonicalNodeKey_(const std::vector<int>& nodes,
                                std::vector<int>& out) const;

    uint64_t hashSolution_(const Solution& s) const;
    void tryAddToArchive_(const EvaluatedSolution& s);
    void trimArchive_();
    void depositArchiveMemory_();

    // hash helper
    static uint64_t splitmix64_(uint64_t x);


    size_t idxMove(int car, int i, int j) const;
    size_t idxReturn(int car, int from, int to) const;

    int computeGbPeriodAuto() const;
    bool isGbIteration(int iterationIndex) const;

    void ensureTauBoundsKnownFromGlobalBest(double Lbest);
    void clampAll();
    inline double clampOne(double x) const;

    void evaporate(double rhoEff);
    void depositSolution(const Solution& s, double solCost, double weight = 1.0);

    void trailSmoothing(double gamma);

    static inline double clamp01(double x);
    void updateTauBaseFromRestartTarget();

    double computeTauMinFromPBest(double pBest) const;

    bool hasReturn_ = false;
    std::vector<double> tauReturn_; // allocated only when hasReturn_ == true

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    const aco_cli::ParamsData& params_;

    int N_ = 0;
    int C_ = 0;

    std::vector<double> tau_;

    bool boundsKnown_ = false;
    double tauMax_ = 1.0;
    double tauMin_ = 0.0;
    double tauBase_ = 1.0;

    int noImprove_ = 0;

    // smoothing
    double smoothingGamma_ = 0.3;

    // exploration phase
    ExplorationOptions exploration_;
    int explorationLeft_ = 0;
    
    // Track if smoothing was applied in the last iteration
    bool lastSmoothingApplied_ = false;
};

class PheromoneModelTBASMove final : public IPheromoneModel
{
public:
    PheromoneModelTBASMove(std::shared_ptr<const cars_tsplib::Instance> inst,
                           const aco_cli::ParamsData& params,
                           ExplorationOptions exploration = {},
                           ArchiveOptions archive = {});

    void reset() override;
    double tauMove(int car, int i, int j) const override;
    double tauReturn(int car, int from, int to) const override;
    bool hasReturnPheromones() const override { return hasReturn_; }

    void onIterationEnd(int iterationIndex,
                        const EvaluatedSolution& iterBest,
                        const EvaluatedSolution& globalBest,
                        bool globalImproved) override;

    void onIterationEndRanked(int iterationIndex,
                          const std::vector<EvaluatedSolution>& ranked,
                          const EvaluatedSolution& globalBest,
                          bool globalImproved) override;
    
    bool wasSmoothingApplied() const override { return lastSmoothingApplied_; }
    int getStagnationCounter() const override { return noImprove_; }
    void reducePheromonesOnTour(const Solution& solution, double gamma) override;
    void addToArchive(const EvaluatedSolution& solution) override;
    double getBestArchivedCost() const override;
    bool getBestArchivedSolution(EvaluatedSolution& out) const override;

private:
    // Archive memory (shared with MMAS)
    ArchiveOptions archiveOpt_;
    struct ArchiveEntry
    {
        EvaluatedSolution s;
        uint64_t h = 0;
    };
    std::vector<ArchiveEntry> archive_;
    size_t archiveCursor_ = 0;

    bool treatAsAsymmetric_() const;
    bool treatAsSymmetric_() const { return !treatAsAsymmetric_(); }
    void buildCanonicalNodeKey_(const std::vector<int>& nodes,
                                std::vector<int>& out) const;
    uint64_t hashSolution_(const Solution& s) const;
    void tryAddToArchive_(const EvaluatedSolution& s);
    void trimArchive_();
    void depositArchiveMemory_();
    static uint64_t splitmix64_(uint64_t x);

    size_t idxMove(int car, int i, int j) const;
    size_t idxReturn(int car, int from, int to) const;

    // TBAS-specific methods
    void ensureTauBoundsKnownFromGlobalBest(double Lbest);
    void clampAll();
    inline double clampOne(double x) const;
    void depositSolution(const Solution& s, double solCost, double weight = 1.0);
    void performClipping();  // TBAS clipping procedure

    // Clipping helper functions (notation from paper)
    inline double clampLower(double x, double l) const { return std::max(x, l); }
    inline double clampUpper(double x, double u) const { return std::min(x, u); }
    inline double clampBounds(double x, double l, double u) const { return clampUpper(clampLower(x, l), u); }

    static inline double clamp01(double x);

    bool hasReturn_ = false;
    std::vector<double> tauReturn_;

private:
    std::shared_ptr<const cars_tsplib::Instance> inst_;
    const aco_cli::ParamsData& params_;

    int N_ = 0;
    int C_ = 0;
    std::vector<double> tau_;

    // TBAS: Three bounds
    bool boundsKnown_ = false;
    double tauLB_ = 0.0;  // Lower bound
    double tauUB_ = 1.0;  // Upper bound
    double tauCB_ = 1.0;  // Clipping bound (τ_CB = ω · τ_UB, where ω ∈ [τ_LB/τ_UB, 1])

    // TBAS: Q variable (Q_i = Q_{i-1} / (1 - ρ))
    double Q_ = 1.0;

    int noImprove_ = 0;
    double smoothingGamma_ = 0.3;
    ExplorationOptions exploration_;
    int explorationLeft_ = 0;
    bool lastSmoothingApplied_ = false;
};

} // namespace aco
