
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <parser.hpp>
#include <params_cli.hpp>
#include <pheromone_model.hpp>
#include <parser.hpp>
#include <params_cli.hpp>
#include <cassert>

namespace aco 
{

PheromoneModelMMASMove::PheromoneModelMMASMove(std::shared_ptr<const cars_tsplib::Instance> inst,
                                               const aco_cli::ParamsData& params,
                                               ExplorationOptions exploration,
                                               ArchiveOptions archive)
    : inst_(std::move(inst))
    , params_(params)
    , exploration_(exploration)
    , archiveOpt_(archive)
{
    if (!inst_) throw std::runtime_error("PheromoneModelMMASMove: inst is null");
    N_ = inst_->n();
    C_ = inst_->cars();

    // TSP: returnCostPerCar is empty => false
    hasReturn_ = inst_->hasReturnCosts();

    smoothingGamma_ = clamp01(params_.smoothingGammaCars);

    reset();
}

uint64_t PheromoneModelMMASMove::splitmix64_(uint64_t x)
{
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31);
}

void PheromoneModelMMASMove::buildCanonicalNodeKey_(const std::vector<int>& nodes,
                                                    std::vector<int>& out) const
{
    if (nodes.empty()) {
        out.clear();
        return;
    }

    out = nodes;

    if (!archiveOpt_.canonicalizeReverseIfSymmetric) return;
    if (!treatAsSymmetric_()) return;

    // node[0]==0 is already the canonical start.
    std::vector<int> rev(nodes.size());
    rev[0] = nodes[0];

    // Reverse the cycle with fixed start 0:
    // forward: 0, a1, a2, ..., a_{N-1}
    // reverse: 0, a_{N-1}, ..., a1
    for (size_t i = 1; i < nodes.size(); ++i)
        rev[i] = nodes[nodes.size() - i];

    if (rev < out) out.swap(rev);
}

bool PheromoneModelMMASMove::treatAsAsymmetric_() const
{
    // Thread-safety: inst_ may be null in a multi-threaded context.
    if (!inst_) {
        // Assume symmetric (default) if inst_ is null.
        return false;
    }

    // Same logic as Colony::preferThreeOptPolish(), without the polish context.
    const bool isATSP = (inst_->type == cars_tsplib::ProblemType::ATSP);

    const bool ewAsym =
        inst_->edgeWeightIsSymmetric.has_value() && !*inst_->edgeWeightIsSymmetric;

    const bool rrAsym =
        inst_->returnRateIsAsymmetric.has_value() && *inst_->returnRateIsAsymmetric;

    return isATSP || ewAsym || rrAsym;
}


uint64_t PheromoneModelMMASMove::hashSolution_(const Solution& s) const
{
    if (s.node.empty()) {
        // Default hash for an empty solution.
        return 0x123456789abcdef0ull;
    }

    uint64_t h = 0x123456789abcdef0ull;

    std::vector<int> key;
    buildCanonicalNodeKey_(s.node, key);

    for (size_t i = 0; i < key.size(); ++i)
    {
        uint64_t v = (uint64_t)key[i] + 0x9e3779b97f4a7c15ull * (i + 1);
        h ^= splitmix64_(v + (h << 1));
    }

    if (archiveOpt_.hashIncludesCars)
    {
        for (size_t i = 0; i < s.car.size(); ++i)
        {
            uint64_t v = (uint64_t)s.car[i] + 0xbf58476d1ce4e5b9ull * (i + 1);
            h ^= splitmix64_(v + (h << 1));
        }
    }

    return h;
}

void PheromoneModelMMASMove::tryAddToArchive_(const EvaluatedSolution& es)
{
    if (!archiveOpt_.enabled) return;
    if (!(es.cost > 0.0) || !std::isfinite(es.cost)) return;

    // Skip solutions inconsistent with N_ (don't add to archive)
    if ((int)es.sol.node.size() != N_ || (int)es.sol.car.size() != N_)
    {
        return;
    }

    uint64_t h = hashSolution_(es.sol);

    // Linear scan is fine for archive size <= 40 (simpler/faster than unordered_set here)
    for (const auto& e : archive_)
        if (e.h == h) return;

    archive_.push_back(ArchiveEntry{es, h});
}

void PheromoneModelMMASMove::trimArchive_()
{
    if (!archiveOpt_.enabled) return;
    if ((int)archive_.size() <= archiveOpt_.size) return;

    std::sort(archive_.begin(), archive_.end(),
              [](const ArchiveEntry& a, const ArchiveEntry& b) { return a.s.cost < b.s.cost; });

    archive_.resize((size_t)archiveOpt_.size);
    if (archiveCursor_ >= archive_.size()) archiveCursor_ = 0;
}

void PheromoneModelMMASMove::depositArchiveMemory_()
{
    if (!archiveOpt_.enabled) return;
    if (archive_.empty()) return;

    // Keep archive sorted (best first)
    std::sort(archive_.begin(), archive_.end(),
              [](const ArchiveEntry& a, const ArchiveEntry& b) { return a.s.cost < b.s.cost; });

    int k = std::max(1, archiveOpt_.depositPerIter);
    for (int t = 0; t < k; ++t)
    {
        const auto& es = archive_[archiveCursor_];
        depositSolution(es.s.sol, es.s.cost, archiveOpt_.depositWeight);

        archiveCursor_++;
        if (archiveCursor_ >= archive_.size()) archiveCursor_ = 0;
    }
}

inline double PheromoneModelMMASMove::clamp01(double x)
{
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

void PheromoneModelMMASMove::updateTauBaseFromRestartTarget()
{
    // Only called once boundsKnown_/tauMin_/tauMax_ are valid.
    if (params_.restartTargetCars == aco_cli::RestartTarget::TauMax) {
        tauBase_ = tauMax_;
    } else { // Mid (default)
        tauBase_ = (tauMin_ + tauMax_) * 0.5;
    }
}

void PheromoneModelMMASMove::reset()
{
    // Main pheromone matrices
    tau_.assign((size_t)C_ * (size_t)N_ * (size_t)N_, 1.0);

    if (hasReturn_)
    {
        tauReturn_.assign((size_t)C_ * (size_t)N_ * (size_t)N_, 1.0);
    }
    else
    {
        tauReturn_.clear();
    }

    // MMAS bounds / state
    boundsKnown_ = false;
    tauMax_ = 1.0;
    tauMin_ = 0.0;
    tauBase_ = 1.0;

    noImprove_ = 0;
    explorationLeft_ = 0;
    lastSmoothingApplied_ = false;

    smoothingGamma_ = clamp01(params_.smoothingGammaCars);

    // Archive memory state
    archive_.clear();
    archiveCursor_ = 0;
}


size_t PheromoneModelMMASMove::idxMove(int car, int i, int j) const 
{
    return ((size_t)car * (size_t)N_ + (size_t)i) * (size_t)N_ + (size_t)j;
}

size_t PheromoneModelMMASMove::idxReturn(int car, int from, int to) const 
{
    return ((size_t)car * (size_t)N_ + (size_t)from) * (size_t)N_ + (size_t)to;
}

double PheromoneModelMMASMove::tauMove(int car, int i, int j) const 
{
    return tau_[idxMove(car, i, j)];
}

double PheromoneModelMMASMove::tauReturn(int car, int from, int to) const
{
    if (!hasReturn_) return 1.0; // TSP: neutral value
    return tauReturn_[idxReturn(car, from, to)];
}

int PheromoneModelMMASMove::computeGbPeriodAuto() const
{
    // Good default for stagnation ~30-50 (and LK)
    int x = params_.stagnation / 3;
    x = std::max(5, std::min(20, x));
    return x;
}

bool PheromoneModelMMASMove::isGbIteration(int iterationIndex) const
{
    const int gb = params_.gbPeriod;
    if (gb < 0) return false;            // off
    int period = (gb == 0) ? computeGbPeriodAuto() : gb;
    if (period <= 0) return false;
    return (iterationIndex % period) == 0;
}

void PheromoneModelMMASMove::ensureTauBoundsKnownFromGlobalBest(double Lbest)
{
    if (!(Lbest > 0.0) || !std::isfinite(Lbest)) return;

    const double rho = params_.rhoCars;
    if (!(rho > 0.0) || !(rho < 1.0)) return;

    const double newTauMax = 1.0 / (rho * Lbest);
    if (!(newTauMax > 0.0) || !std::isfinite(newTauMax)) return;

    boundsKnown_ = true;
    tauMax_ = newTauMax;

    // tauMin: pBest overrides maxMin (if enabled)
    if (params_.pBestCars > 0.0) {
        tauMin_ = computeTauMinFromPBest(params_.pBestCars);
    }
    else if (params_.maxMinCars > 0.0) {
        tauMin_ = tauMax_ / params_.maxMinCars;
    }
    else {
        tauMin_ = 0.0; // max-min disabled and pBest disabled => floor 0
    }

    // Restart/smoothing target depends on restartTargetCars
    updateTauBaseFromRestartTarget();

    clampAll();
}


double PheromoneModelMMASMove::computeTauMinFromPBest(double pBest) const
{
    // pBest expected in (0,1)
    if (!(pBest > 0.0 && pBest < 1.0) || !std::isfinite(pBest)) return 0.0;
    if (!boundsKnown_) return 0.0;

    const double n = (N_ > 0) ? (double)N_ : 1.0;

    // avg number of choices per decision:
    // - if candidate list used, favorites is a good proxy
    // - otherwise fallback to N/2 (classic MMAS heuristic)
    double avg = 0.0;
    if (params_.favorites > 1) avg = (double)params_.favorites;
    else                       avg = std::max(2.0, n * 0.5);

    // p^(1/n)
    const double pn = std::pow(pBest, 1.0 / n);

    const double denom = (avg - 1.0) * pn;
    if (!(denom > 0.0) || !std::isfinite(denom)) return 0.0;

    double tmin = tauMax_ * (1.0 - pn) / denom;

    // clamp to [0, tauMax] just in case
    if (!std::isfinite(tmin) || tmin < 0.0) tmin = 0.0;
    if (tmin > tauMax_) tmin = tauMax_;
    return tmin;
}


inline double PheromoneModelMMASMove::clampOne(double x) const 
{
    if (!boundsKnown_) return x;
    if (x < tauMin_) return tauMin_;
    if (x > tauMax_) return tauMax_;
    return x;
}

void PheromoneModelMMASMove::clampAll() 
{
    if (!boundsKnown_) return;
    for (double& t : tau_) t = clampOne(t);
    if (hasReturn_) 
    {
        for (double& t : tauReturn_)
            t = clampOne(t);
    }
}

void PheromoneModelMMASMove::evaporate(double rhoEff) 
{
    const double oneMinus = 1.0 - rhoEff;
    for (double& t : tau_) 
    {
        t *= oneMinus;
        if (boundsKnown_) t = clampOne(t);
    }
    if (hasReturn_) 
    {
        for (double& t : tauReturn_) {
            t *= oneMinus;
            if (boundsKnown_) t = clampOne(t);
        }
    }
}

void PheromoneModelMMASMove::depositSolution(const Solution& s, double solCost, double weight) 
{
    if (!(solCost > 0.0) || !std::isfinite(solCost)) return;
    if (!(weight > 0.0) || !std::isfinite(weight)) return;

    if ((int)s.node.size() != N_ || (int)s.car.size() != N_) 
    {
        throw std::runtime_error("depositSolution: node/car size mismatch (node.size()=" + 
                                 std::to_string(s.node.size()) + ", car.size()=" + 
                                 std::to_string(s.car.size()) + ", N_=" + std::to_string(N_) + ")");
    }

    const double delta = weight / solCost;

    for (int k = 0; k < N_; ++k) 
    {
        const int i = s.node[k];
        const int j = s.node[(k + 1) % N_];
        const int c = s.car[k];

        if (i < 0 || i >= N_ || j < 0 || j >= N_ || c < 0 || c >= C_) {
            throw std::runtime_error("depositSolution: invalid i/j/c index");
        }

        const size_t im = idxMove(c, i, j);
        if (im >= tau_.size()) {
            throw std::runtime_error("depositSolution: idxMove OOB");
        }

        double& t = tau_[im];
        t += delta;
        if (boundsKnown_) t = clampOne(t);
    }

    if (!hasReturn_) return;

    if (tauReturn_.size() != (size_t)C_ * (size_t)N_ * (size_t)N_) {
        throw std::runtime_error("depositSolution: tauReturn_ bad size");
    }

    int currentCar  = s.car[0];
    int lastCarNode = s.node[0];

    if (currentCar < 0 || currentCar >= C_) {
        throw std::runtime_error("depositSolution: car[0] OOB");
    }
    if (lastCarNode < 0 || lastCarNode >= N_) {
        throw std::runtime_error("depositSolution: node[0] OOB");
    }

    for (int i = 0; i < N_; ++i) {
        const int currentNode = s.node[i];
        const int lastCar = (i > 0) ? s.car[i - 1] : currentCar;
        currentCar = s.car[i];

        if (currentNode < 0 || currentNode >= N_) {
            throw std::runtime_error("depositSolution: currentNode OOB");
        }
        if (lastCar < 0 || lastCar >= C_) {
            throw std::runtime_error("depositSolution: lastCar OOB");
        }
        if (currentCar < 0 || currentCar >= C_) {
            throw std::runtime_error("depositSolution: currentCar OOB");
        }

        if (currentCar != lastCar) 
        {
            const size_t ir = idxReturn(lastCar, currentNode, lastCarNode);
            if (ir >= tauReturn_.size()) {
                throw std::runtime_error("depositSolution: idxReturn OOB");
            }

            double& tr = tauReturn_[idxReturn(lastCar, currentNode, lastCarNode)];
            tr += delta;
            if (boundsKnown_) tr = clampOne(tr);
            lastCarNode = currentNode;
        }
    }

    {
        const int startNode = s.node[0];
        const size_t ir = idxReturn(currentCar, startNode, lastCarNode);
        if (ir >= tauReturn_.size()) {
            throw std::runtime_error("depositSolution: final idxReturn OOB");
        }

        double& tr = tauReturn_[idxReturn(currentCar, startNode, lastCarNode)];
        tr += delta;
        if (boundsKnown_) tr = clampOne(tr);
    }
}

void PheromoneModelMMASMove::trailSmoothing(double gamma) 
{
    if (!boundsKnown_) return; // nema smisla prije nego znamo tauMax
    gamma = std::max(0.0, std::min(1.0, gamma));

    for (double& t : tau_) 
    {
        t = (1.0 - gamma) * t + gamma * tauBase_;
        t = clampOne(t);
    }
    if (hasReturn_) 
    {
        for (double& t : tauReturn_) {
            t = (1.0 - gamma) * t + gamma * tauBase_;
            t = clampOne(t);
        }
    }
}

void PheromoneModelMMASMove::onIterationEnd(int iterationIndex,
                                            const EvaluatedSolution& iterBest,
                                            const EvaluatedSolution& globalBest,
                                            bool globalImproved)
{
    std::vector<EvaluatedSolution> ranked;
    ranked.reserve(1);
    ranked.push_back(iterBest);
    onIterationEndRanked(iterationIndex, ranked, globalBest, globalImproved);
}


void PheromoneModelMMASMove::onIterationEndRanked(int iterationIndex,
                                                  const std::vector<EvaluatedSolution>& ranked,
                                                  const EvaluatedSolution& globalBest,
                                                  bool globalImproved)
{
    // As soon as we have the first globalBest, the bounds are known.
    if (!boundsKnown_ && std::isfinite(globalBest.cost) && globalBest.cost > 0.0) {
        ensureTauBoundsKnownFromGlobalBest(globalBest.cost);
    }

    // Archive (global memory): update pool (node-only hash recommended).
    // ranked is assumed DP-consistent (guaranteed by Colony before this call).
    if (archiveOpt_.enabled)
    {
        tryAddToArchive_(globalBest);
        for (const auto& es : ranked)
            tryAddToArchive_(es);
        trimArchive_();
    }

    // Update stagnation counter
    if (globalImproved) noImprove_ = 0;
    else               noImprove_++;

    // Stagnation handling: smoothing toward tauBase (+ optional exploration phase).
    // Trigger smoothing earlier (1.5x instead of 2x) for a faster reset in late stagnation.
    const int stagPher = (params_.stagnation > 0) ? (int)(1.5 * params_.stagnation) : 0;

    lastSmoothingApplied_ = false; // reset for this iteration
    if (stagPher > 0 && noImprove_ >= stagPher)
    {
        trailSmoothing(smoothingGamma_);
        lastSmoothingApplied_ = true;
        noImprove_ = 0;

        if (exploration_.enabled) {
            explorationLeft_ = exploration_.iters;
        }
    }

    // Effective rho (exploration phase temporarily boosts evaporation)
    double rhoEff = params_.rhoCars;
    if (explorationLeft_ > 0) {
        rhoEff = std::min(0.99, rhoEff * exploration_.rhoMultiplier);
        explorationLeft_--;
    }

    // MMAS update:
    // 1) evaporate
    evaporate(rhoEff);

    // 2) deposit (GB on schedule, otherwise ranked top-W).
    // If the global best improved this iteration, always deposit it immediately.
    const bool useGb = isGbIteration(iterationIndex);

    const int W = (int)ranked.size();

    // Reduce the normal bonus when W reaches 7-8 so the update doesn't get too wide.
    const double bonusNormal = (W <= 6) ? 0.5 : 0.35;

    auto isValidSolution = [this](const Solution& s) -> bool {
        return (int)s.node.size() == N_ && (int)s.car.size() == N_;
    };

    if (globalImproved)
    {
        // New global best is deposited immediately (only if consistent).
        if (isValidSolution(globalBest.sol)) {
            depositSolution(globalBest.sol, globalBest.cost, 1.0);
        }

        const double bonusOnGlobalImproved = 0.2;

        // Skip ranked[0] if it's effectively identical to globalBest (avoid double deposit).
        const double eps = 1e-12;
        int r0 = 0;

        if (W > 0 && std::isfinite(ranked[0].cost) && std::isfinite(globalBest.cost))
        {
            const double scale = std::max(1.0, std::abs(globalBest.cost));
            if (std::abs(ranked[0].cost - globalBest.cost) <= eps * scale)
                r0 = 1;
        }

        for (int r = r0; r < W; ++r)
        {
            if (r < (int)ranked.size() && isValidSolution(ranked[r].sol)) {
                const double w = bonusOnGlobalImproved * (double)(W - r) / (double)W;
                depositSolution(ranked[r].sol, ranked[r].cost, w);
            }
        }
    }
    else if (useGb || ranked.empty())
    {
        // Classic MMAS "anchor" deposit (only if consistent).
        if (isValidSolution(globalBest.sol)) {
            depositSolution(globalBest.sol, globalBest.cost, 1.0);
        }
    }
    else
    {
        // Ranked update: iteration-best + bonus to the rest (only if consistent).
        if (!ranked.empty() && isValidSolution(ranked[0].sol)) {
            depositSolution(ranked[0].sol, ranked[0].cost, 1.0);
        }

        for (int r = 1; r < W; ++r)
        {
            if (r < (int)ranked.size() && isValidSolution(ranked[r].sol)) {
                const double w = bonusNormal * (double)(W - r) / (double)W;
                depositSolution(ranked[r].sol, ranked[r].cost, w);
            }
        }
    }

    // If the global best improved, refresh tauMax/tauMin (standard MMAS).
    if (globalImproved) {
        ensureTauBoundsKnownFromGlobalBest(globalBest.cost);
    }

    if (archiveOpt_.enabled && boundsKnown_)
    {
        depositArchiveMemory_();
    }
}

void PheromoneModelMMASMove::reducePheromonesOnTour(const Solution& solution, double gamma)
{
    // Selectively reduces pheromones on a tour (smoothing):
    // tau = (1 - gamma) * tau + gamma * tauBase
    // gamma = 0.0 => full reset to tauBase; gamma = 1.0 => no change.

    if ((int)solution.node.size() != N_ || (int)solution.car.size() != N_) return;

    const double oneMinusGamma = 1.0 - gamma;
    const double gammaTauBase = gamma * tauBase_;

    // Reduce pheromones on move edges in the tour.
    for (int k = 0; k < N_; ++k)
    {
        const int i = solution.node[k];
        const int j = solution.node[(k + 1) % N_];
        const int c = solution.car[k];
        
        if (i < 0 || i >= N_ || j < 0 || j >= N_ || c < 0 || c >= C_) continue;
        
        const size_t im = idxMove(c, i, j);
        if (im < tau_.size())
        {
            double& t = tau_[im];
            t = oneMinusGamma * t + gammaTauBase;
            if (boundsKnown_) t = clampOne(t);
        }
    }
    
    // Reduce pheromones on return edges in the tour (if present).
    if (hasReturn_ && tauReturn_.size() == (size_t)C_ * (size_t)N_ * (size_t)N_)
    {
        int currentCar = solution.car[0];
        int lastCarNode = solution.node[0];
        
        for (int i = 0; i < N_; ++i)
        {
            const int currentNode = solution.node[i];
            const int lastCar = (i > 0) ? solution.car[i - 1] : currentCar;
            currentCar = solution.car[i];
            
            if (currentCar != lastCar)
            {
                const size_t ir = idxReturn(lastCar, currentNode, lastCarNode);
                if (ir < tauReturn_.size())
                {
                    double& tr = tauReturn_[ir];
                    tr = oneMinusGamma * tr + gammaTauBase;
                    if (boundsKnown_) tr = clampOne(tr);
                }
                lastCarNode = currentNode;
            }
        }
        
        // Final return edge
        const int startNode = solution.node[0];
        const size_t ir = idxReturn(currentCar, startNode, lastCarNode);
        if (ir < tauReturn_.size())
        {
            double& tr = tauReturn_[ir];
            tr = oneMinusGamma * tr + gammaTauBase;
            if (boundsKnown_) tr = clampOne(tr);
        }
    }
}

void PheromoneModelMMASMove::addToArchive(const EvaluatedSolution& solution)
{
    tryAddToArchive_(solution);
    trimArchive_();
}

double PheromoneModelMMASMove::getBestArchivedCost() const
{
    if (!archiveOpt_.enabled || archive_.empty())
        return std::numeric_limits<double>::infinity();

    // Re-scan rather than trust sort order, since a new solution may have
    // been added between calls.
    double bestCost = std::numeric_limits<double>::infinity();
    for (const auto& entry : archive_)
    {
        if ((int)entry.s.sol.node.size() != N_ || (int)entry.s.sol.car.size() != N_)
            continue;
        if (std::isfinite(entry.s.cost) && entry.s.cost < bestCost)
            bestCost = entry.s.cost;
    }
    
    return bestCost;
}

bool PheromoneModelMMASMove::getBestArchivedSolution(EvaluatedSolution& out) const
{
    if (!archiveOpt_.enabled || archive_.empty())
        return false;

    const ArchiveEntry* bestEntry = nullptr;
    double bestCost = std::numeric_limits<double>::infinity();

    for (const auto& entry : archive_)
    {
        // Skip entries with inconsistent size (guards against corruption/races).
        if ((int)entry.s.sol.node.size() != N_ || (int)entry.s.sol.car.size() != N_)
            continue;
        if (std::isfinite(entry.s.cost) && entry.s.cost < bestCost)
        {
            bestCost = entry.s.cost;
            bestEntry = &entry;
        }
    }
    
    if (bestEntry)
    {
        out = bestEntry->s;
        return true;
    }
    
    return false;
}




// ============================================================================
// TBAS (Three Bounds Ant System) Implementation
// ============================================================================

PheromoneModelTBASMove::PheromoneModelTBASMove(std::shared_ptr<const cars_tsplib::Instance> inst,
                                               const aco_cli::ParamsData& params,
                                               ExplorationOptions exploration,
                                               ArchiveOptions archive)
    : inst_(std::move(inst))
    , params_(params)
    , exploration_(exploration)
    , archiveOpt_(archive)
{
    if (!inst_) throw std::runtime_error("PheromoneModelTBASMove: inst is null");
    N_ = inst_->n();
    C_ = inst_->cars();

    hasReturn_ = inst_->hasReturnCosts();
    smoothingGamma_ = clamp01(params_.smoothingGammaCars);

    reset();
}

uint64_t PheromoneModelTBASMove::splitmix64_(uint64_t x)
{
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31);
}

void PheromoneModelTBASMove::buildCanonicalNodeKey_(const std::vector<int>& nodes,
                                                    std::vector<int>& out) const
{
    if (nodes.empty()) {
        out.clear();
        return;
    }

    out = nodes;

    if (!archiveOpt_.canonicalizeReverseIfSymmetric) return;
    if (!treatAsSymmetric_()) return;

    std::vector<int> rev(nodes.size());
    rev[0] = nodes[0];

    for (size_t i = 1; i < nodes.size(); ++i)
        rev[i] = nodes[nodes.size() - i];

    if (rev < out) out.swap(rev);
}

bool PheromoneModelTBASMove::treatAsAsymmetric_() const
{
    if (!inst_) return false;

    const bool isATSP = (inst_->type == cars_tsplib::ProblemType::ATSP);
    const bool ewAsym =
        inst_->edgeWeightIsSymmetric.has_value() && !*inst_->edgeWeightIsSymmetric;
    const bool rrAsym =
        inst_->returnRateIsAsymmetric.has_value() && *inst_->returnRateIsAsymmetric;

    return isATSP || ewAsym || rrAsym;
}

uint64_t PheromoneModelTBASMove::hashSolution_(const Solution& s) const
{
    if (s.node.empty()) {
        return 0x123456789abcdef0ull;
    }

    uint64_t h = 0x123456789abcdef0ull;

    std::vector<int> key;
    buildCanonicalNodeKey_(s.node, key);

    for (size_t i = 0; i < key.size(); ++i)
    {
        uint64_t v = (uint64_t)key[i] + 0x9e3779b97f4a7c15ull * (i + 1);
        h ^= splitmix64_(v + (h << 1));
    }

    if (archiveOpt_.hashIncludesCars)
    {
        for (size_t i = 0; i < s.car.size(); ++i)
        {
            uint64_t v = (uint64_t)s.car[i] + 0xbf58476d1ce4e5b9ull * (i + 1);
            h ^= splitmix64_(v + (h << 1));
        }
    }

    return h;
}

void PheromoneModelTBASMove::tryAddToArchive_(const EvaluatedSolution& es)
{
    if (!archiveOpt_.enabled) return;
    if (!(es.cost > 0.0) || !std::isfinite(es.cost)) return;

    if ((int)es.sol.node.size() != N_ || (int)es.sol.car.size() != N_) 
    {
        return;
    }

    uint64_t h = hashSolution_(es.sol);

    for (const auto& e : archive_)
        if (e.h == h) return;

    archive_.push_back(ArchiveEntry{es, h});
}

void PheromoneModelTBASMove::trimArchive_()
{
    if (!archiveOpt_.enabled) return;
    if ((int)archive_.size() <= archiveOpt_.size) return;

    std::sort(archive_.begin(), archive_.end(),
              [](const ArchiveEntry& a, const ArchiveEntry& b) { return a.s.cost < b.s.cost; });

    archive_.resize((size_t)archiveOpt_.size);
    if (archiveCursor_ >= archive_.size()) archiveCursor_ = 0;
}

void PheromoneModelTBASMove::depositArchiveMemory_()
{
    if (!archiveOpt_.enabled) return;
    if (archive_.empty()) return;

    std::sort(archive_.begin(), archive_.end(),
              [](const ArchiveEntry& a, const ArchiveEntry& b) { return a.s.cost < b.s.cost; });

    int k = std::max(1, archiveOpt_.depositPerIter);
    for (int t = 0; t < k; ++t)
    {
        const auto& es = archive_[archiveCursor_];
        depositSolution(es.s.sol, es.s.cost, archiveOpt_.depositWeight);

        archiveCursor_++;
        if (archiveCursor_ >= archive_.size()) archiveCursor_ = 0;
    }
}

inline double PheromoneModelTBASMove::clamp01(double x)
{
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

void PheromoneModelTBASMove::reset()
{
    // TBAS: Initial pheromone τ_0 = τ_LB
    // We'll set tauLB_ after bounds are known, for now use a default
    tau_.assign((size_t)C_ * (size_t)N_ * (size_t)N_, 1.0);

    if (hasReturn_)
    {
        tauReturn_.assign((size_t)C_ * (size_t)N_ * (size_t)N_, 1.0);
    }
    else
    {
        tauReturn_.clear();
    }

    boundsKnown_ = false;
    tauLB_ = 0.0;
    tauUB_ = 1.0;
    tauCB_ = 1.0;

    // TBAS: Q variable initialization
    Q_ = 1.0;

    noImprove_ = 0;
    explorationLeft_ = 0;
    lastSmoothingApplied_ = false;

    smoothingGamma_ = clamp01(params_.smoothingGammaCars);

    archive_.clear();
    archiveCursor_ = 0;
}

size_t PheromoneModelTBASMove::idxMove(int car, int i, int j) const 
{
    return ((size_t)car * (size_t)N_ + (size_t)i) * (size_t)N_ + (size_t)j;
}

size_t PheromoneModelTBASMove::idxReturn(int car, int from, int to) const 
{
    return ((size_t)car * (size_t)N_ + (size_t)from) * (size_t)N_ + (size_t)to;
}

double PheromoneModelTBASMove::tauMove(int car, int i, int j) const 
{
    return tau_[idxMove(car, i, j)];
}

double PheromoneModelTBASMove::tauReturn(int car, int from, int to) const 
{
    if (!hasReturn_) return 1.0;
    return tauReturn_[idxReturn(car, from, to)];
}

void PheromoneModelTBASMove::ensureTauBoundsKnownFromGlobalBest(double Lbest)
{
    if (!(Lbest > 0.0) || !std::isfinite(Lbest)) return;

    const double rho = params_.rhoCars;
    if (!(rho > 0.0) || !(rho < 1.0)) return;

    // TBAS: Similar to MMAS, but we use these as τ_LB and τ_UB
    // τ_UB = 1.0 / (ρ * Lbest) (similar to MMAS tauMax)
    const double newTauUB = 1.0 / (rho * Lbest);
    if (!(newTauUB > 0.0) || !std::isfinite(newTauUB)) return;

    boundsKnown_ = true;
    tauUB_ = newTauUB;

    // τ_LB: use pBest or maxMin if available, otherwise use a fraction of τ_UB
    if (params_.pBestCars > 0.0) {
        // Similar to MMAS pBest calculation
        const double n = (N_ > 0) ? (double)N_ : 1.0;
        double avg = 0.0;
        if (params_.favorites > 1) avg = (double)params_.favorites;
        else                       avg = std::max(2.0, n * 0.5);

        const double pn = std::pow(params_.pBestCars, 1.0 / n);
        const double denom = (avg - 1.0) * pn;
        if ((denom > 0.0) && std::isfinite(denom)) {
            tauLB_ = tauUB_ * (1.0 - pn) / denom;
        } else {
            tauLB_ = tauUB_ / 100.0;  // fallback
        }
    }
    else if (params_.maxMinCars > 0.0) {
        tauLB_ = tauUB_ / params_.maxMinCars;
    }
    else {
        tauLB_ = tauUB_ / 100.0;  // default fallback
    }

    // Ensure τ_LB < τ_UB
    if (tauLB_ >= tauUB_) {
        tauLB_ = tauUB_ * 0.01;
    }

    // TBAS: τ_CB = ω · τ_UB, where ω ∈ [τ_LB/τ_UB, 1]
    // Use ω = 1 - ρ (similar to MMAS relationship)
    const double omega = std::max(tauLB_ / tauUB_, 1.0 - rho);
    tauCB_ = omega * tauUB_;

    // TBAS: Initial pheromone τ_0 = τ_LB
    for (double& t : tau_) {
        t = tauLB_;
    }
    if (hasReturn_) {
        for (double& t : tauReturn_) {
            t = tauLB_;
        }
    }

    clampAll();
}

inline double PheromoneModelTBASMove::clampOne(double x) const 
{
    if (!boundsKnown_) return x;
    if (x < tauLB_) return tauLB_;
    if (x > tauUB_) return tauUB_;
    return x;
}

void PheromoneModelTBASMove::clampAll() 
{
    if (!boundsKnown_) return;
    for (double& t : tau_) t = clampOne(t);
    if (hasReturn_) 
    {
        for (double& t : tauReturn_)
            t = clampOne(t);
    }
}

void PheromoneModelTBASMove::depositSolution(const Solution& s, double solCost, double weight) 
{
    if (!(solCost > 0.0) || !std::isfinite(solCost)) return;
    if (!(weight > 0.0) || !std::isfinite(weight)) return;

    if ((int)s.node.size() != N_ || (int)s.car.size() != N_) 
    {
        throw std::runtime_error("depositSolution: node/car size mismatch");
    }

    // TBAS: Δτ_c = Q / f(s*) for components of solution s*
    // Q is updated each iteration: Q_i = Q_{i-1} / (1 - ρ)
    const double delta = (Q_ * weight) / solCost;

    for (int k = 0; k < N_; ++k) 
    {
        const int i = s.node[k];
        const int j = s.node[(k + 1) % N_];
        const int c = s.car[k];

        if (i < 0 || i >= N_ || j < 0 || j >= N_ || c < 0 || c >= C_) {
            throw std::runtime_error("depositSolution: invalid i/j/c index");
        }

        const size_t im = idxMove(c, i, j);
        if (im >= tau_.size()) {
            throw std::runtime_error("depositSolution: idxMove OOB");
        }

        double& t = tau_[im];
        t += delta;
        if (boundsKnown_) t = clampOne(t);
    }

    if (!hasReturn_) return;

    if (tauReturn_.size() != (size_t)C_ * (size_t)N_ * (size_t)N_) {
        throw std::runtime_error("depositSolution: tauReturn_ bad size");
    }

    int currentCar  = s.car[0];
    int lastCarNode = s.node[0];

    if (currentCar < 0 || currentCar >= C_) {
        throw std::runtime_error("depositSolution: car[0] OOB");
    }
    if (lastCarNode < 0 || lastCarNode >= N_) {
        throw std::runtime_error("depositSolution: node[0] OOB");
    }

    for (int i = 0; i < N_; ++i) {
        const int currentNode = s.node[i];
        const int lastCar = (i > 0) ? s.car[i - 1] : currentCar;
        currentCar = s.car[i];

        if (currentNode < 0 || currentNode >= N_) {
            throw std::runtime_error("depositSolution: currentNode OOB");
        }
        if (lastCar < 0 || lastCar >= C_) {
            throw std::runtime_error("depositSolution: lastCar OOB");
        }
        if (currentCar < 0 || currentCar >= C_) {
            throw std::runtime_error("depositSolution: currentCar OOB");
        }

        if (currentCar != lastCar) 
        {
            const size_t ir = idxReturn(lastCar, currentNode, lastCarNode);
            if (ir >= tauReturn_.size()) {
                throw std::runtime_error("depositSolution: idxReturn OOB");
            }

            double& tr = tauReturn_[idxReturn(lastCar, currentNode, lastCarNode)];
            tr += delta;
            if (boundsKnown_) tr = clampOne(tr);
            lastCarNode = currentNode;
        }
    }

    {
        const int startNode = s.node[0];
        const size_t ir = idxReturn(currentCar, startNode, lastCarNode);
        if (ir >= tauReturn_.size()) {
            throw std::runtime_error("depositSolution: final idxReturn OOB");
        }

        double& tr = tauReturn_[idxReturn(currentCar, startNode, lastCarNode)];
        tr += delta;
        if (boundsKnown_) tr = clampOne(tr);
    }
}

void PheromoneModelTBASMove::performClipping()
{
    if (!boundsKnown_) return;

    // Find maximum pheromone value
    double tauMax = tauLB_;
    for (const double& t : tau_) {
        if (t > tauMax) tauMax = t;
    }
    if (hasReturn_) {
        for (const double& t : tauReturn_) {
            if (t > tauMax) tauMax = t;
        }
    }

    // TBAS: Clipping only if τ_max > τ_UB
    if (tauMax <= tauUB_) return;

    // TBAS clipping procedure: τ'' = [ω · τ']_{τ_LB}^{τ_UB}
    // where ω is determined by τ_CB = ω · τ_UB
    const double omega = tauCB_ / tauUB_;
    
    // Apply clipping: multiply by ω, then clamp to [τ_LB, τ_UB]
    for (double& t : tau_) {
        t = omega * t;
        t = clampBounds(t, tauLB_, tauUB_);
    }
    
    if (hasReturn_) {
        for (double& t : tauReturn_) {
            t = omega * t;
            t = clampBounds(t, tauLB_, tauUB_);
        }
    }
}

void PheromoneModelTBASMove::onIterationEnd(int iterationIndex,
                                            const EvaluatedSolution& iterBest,
                                            const EvaluatedSolution& globalBest,
                                            bool globalImproved)
{
    std::vector<EvaluatedSolution> ranked;
    ranked.reserve(1);
    ranked.push_back(iterBest);
    onIterationEndRanked(iterationIndex, ranked, globalBest, globalImproved);
}

void PheromoneModelTBASMove::onIterationEndRanked(int iterationIndex,
                                                  const std::vector<EvaluatedSolution>& ranked,
                                                  const EvaluatedSolution& globalBest,
                                                  bool globalImproved)
{
    // Set bounds from first global best
    if (!boundsKnown_ && std::isfinite(globalBest.cost) && globalBest.cost > 0.0) {
        ensureTauBoundsKnownFromGlobalBest(globalBest.cost);
    }

    // Archive update
    if (archiveOpt_.enabled)
    {
        tryAddToArchive_(globalBest);
        for (const auto& es : ranked)
            tryAddToArchive_(es);
        trimArchive_();
    }

    // Update stagnation counter
    if (globalImproved) noImprove_ = 0;
    else               noImprove_++;

    // TBAS: Update Q variable: Q_i = Q_{i-1} / (1 - ρ)
    const double rho = params_.rhoCars;
    if (rho > 0.0 && rho < 1.0) {
        Q_ = Q_ / (1.0 - rho);
    }

    // TBAS: NO EVAPORATION (unlike MMAS)
    // We only deposit and clip

    // Deposit solutions
    const bool useGb = false;  // TBAS doesn't use gbPeriod like MMAS
    const int W = (int)ranked.size();

    const double bonusNormal = (W <= 6) ? 0.5 : 0.35;

    auto isValidSolution = [this](const Solution& s) -> bool {
        return (int)s.node.size() == N_ && (int)s.car.size() == N_;
    };

    if (globalImproved)
    {
        if (isValidSolution(globalBest.sol)) {
            depositSolution(globalBest.sol, globalBest.cost, 1.0);
        }

        const double bonusOnGlobalImproved = 0.2;
        const double eps = 1e-12;
        int r0 = 0;

        if (W > 0 && std::isfinite(ranked[0].cost) && std::isfinite(globalBest.cost))
        {
            const double scale = std::max(1.0, std::abs(globalBest.cost));
            if (std::abs(ranked[0].cost - globalBest.cost) <= eps * scale)
                r0 = 1;
        }

        for (int r = r0; r < W; ++r)
        {
            if (r < (int)ranked.size() && isValidSolution(ranked[r].sol)) {
                const double w = bonusOnGlobalImproved * (double)(W - r) / (double)W;
                depositSolution(ranked[r].sol, ranked[r].cost, w);
            }
        }
    }
    else
    {
        // Ranked update
        if (!ranked.empty() && isValidSolution(ranked[0].sol)) {
            depositSolution(ranked[0].sol, ranked[0].cost, 1.0);
        }

        for (int r = 1; r < W; ++r)
        {
            if (r < (int)ranked.size() && isValidSolution(ranked[r].sol)) {
                const double w = bonusNormal * (double)(W - r) / (double)W;
                depositSolution(ranked[r].sol, ranked[r].cost, w);
            }
        }
    }

    // Update bounds if global improved
    if (globalImproved) {
        ensureTauBoundsKnownFromGlobalBest(globalBest.cost);
    }

    // TBAS: Perform clipping if needed
    performClipping();

    // Archive memory deposit
    if (archiveOpt_.enabled && boundsKnown_) 
    {
        depositArchiveMemory_();
    }
}

void PheromoneModelTBASMove::reducePheromonesOnTour(const Solution& solution, double gamma)
{
    if ((int)solution.node.size() != N_ || (int)solution.car.size() != N_) return;
    
    if (!boundsKnown_) return;
    
    const double tauBase = tauLB_;  // TBAS uses tauLB_ as base
    const double oneMinusGamma = 1.0 - gamma;
    const double gammaTauBase = gamma * tauBase;
    
    for (int k = 0; k < N_; ++k)
    {
        const int i = solution.node[k];
        const int j = solution.node[(k + 1) % N_];
        const int c = solution.car[k];
        
        if (i < 0 || i >= N_ || j < 0 || j >= N_ || c < 0 || c >= C_) continue;
        
        const size_t im = idxMove(c, i, j);
        if (im < tau_.size())
        {
            double& t = tau_[im];
            t = oneMinusGamma * t + gammaTauBase;
            t = clampOne(t);
        }
    }
    
    if (hasReturn_ && tauReturn_.size() == (size_t)C_ * (size_t)N_ * (size_t)N_)
    {
        int currentCar = solution.car[0];
        int lastCarNode = solution.node[0];
        
        for (int i = 0; i < N_; ++i)
        {
            const int currentNode = solution.node[i];
            const int lastCar = (i > 0) ? solution.car[i - 1] : currentCar;
            currentCar = solution.car[i];
            
            if (currentCar != lastCar)
            {
                const size_t ir = idxReturn(lastCar, currentNode, lastCarNode);
                if (ir < tauReturn_.size())
                {
                    double& tr = tauReturn_[ir];
                    tr = oneMinusGamma * tr + gammaTauBase;
                    tr = clampOne(tr);
                }
                lastCarNode = currentNode;
            }
        }
        
        const int startNode = solution.node[0];
        const size_t ir = idxReturn(currentCar, startNode, lastCarNode);
        if (ir < tauReturn_.size())
        {
            double& tr = tauReturn_[ir];
            tr = oneMinusGamma * tr + gammaTauBase;
            tr = clampOne(tr);
        }
    }
}

void PheromoneModelTBASMove::addToArchive(const EvaluatedSolution& solution)
{
    tryAddToArchive_(solution);
    trimArchive_();
}

double PheromoneModelTBASMove::getBestArchivedCost() const
{
    if (!archiveOpt_.enabled || archive_.empty())
        return std::numeric_limits<double>::infinity();
    
    double bestCost = std::numeric_limits<double>::infinity();
    for (const auto& entry : archive_)
    {
        if ((int)entry.s.sol.node.size() != N_ || (int)entry.s.sol.car.size() != N_)
            continue;
        if (std::isfinite(entry.s.cost) && entry.s.cost < bestCost)
            bestCost = entry.s.cost;
    }
    
    return bestCost;
}

bool PheromoneModelTBASMove::getBestArchivedSolution(EvaluatedSolution& out) const
{
    if (!archiveOpt_.enabled || archive_.empty())
        return false;
    
    const ArchiveEntry* bestEntry = nullptr;
    double bestCost = std::numeric_limits<double>::infinity();
    
    for (const auto& entry : archive_)
    {
        // Skip entries with inconsistent size (guards against corruption/races).
        if ((int)entry.s.sol.node.size() != N_ || (int)entry.s.sol.car.size() != N_)
            continue;
        if (std::isfinite(entry.s.cost) && entry.s.cost < bestCost)
        {
            bestCost = entry.s.cost;
            bestEntry = &entry;
        }
    }

    if (bestEntry)
    {
        out = bestEntry->s;
        return true;
    }

    return false;
}

} // namespace aco
