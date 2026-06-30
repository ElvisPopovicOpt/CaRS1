#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <colony.hpp>
#include <limits>
#include <chrono>
#include <params_cli.hpp>
#include <parser.hpp>
#include <pheromone_model.hpp>
#include <runRecorder.hpp>
#include <CachedFixedTourCarAssigner.hpp>
#include <logger.hpp>
#include <sstream>
#include <intensifier.hpp>
#include <research/intensifier_binomial_logger.hpp>
#include <research/surrogate_correlation_logger.hpp>
#include "LkLiteLocalSearch.hpp"
#include "ThreeOptLiteLocalSearch.hpp"

namespace aco 
{

// Estimates the Weibull shape parameter k from ant costs (moment method: CV = sigma/mu).
// Returns NaN if there are too few samples or no positive costs.
static double weibullKFromCosts(const std::vector<EvaluatedSolution>& ranked)
{
    std::vector<double> costs;
    costs.reserve(ranked.size());
    for (const auto& es : ranked)
        if (std::isfinite(es.cost) && es.cost > 0.0)
            costs.push_back(es.cost);
    const size_t n = costs.size();
    if (n < 3u) return std::numeric_limits<double>::quiet_NaN();

    double sum = 0.0;
    for (double x : costs) sum += x;
    const double mu = sum / static_cast<double>(n);
    if (mu <= 0.0) return std::numeric_limits<double>::quiet_NaN();

    double sumSq = 0.0;
    for (double x : costs) sumSq += (x - mu) * (x - mu);
    const double var = sumSq / static_cast<double>(n);
    if (var <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    const double cvSample = std::sqrt(var) / mu;  // coefficient of variation
    const double cvSq = cvSample * cvSample;

    // Weibull: CV^2 = Gamma(1+2/k)/Gamma(1+1/k)^2 - 1. Solve for k.
    // cv_theory_sq(k) is monotonically decreasing in k.
    auto cvTheorySq = [](double k) -> double {
        if (k <= 0.0) return std::numeric_limits<double>::infinity();
        const double g1 = std::tgamma(1.0 + 2.0 / k);
        const double g2 = std::tgamma(1.0 + 1.0 / k);
        if (g2 <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        const double ratio = g1 / (g2 * g2);
        return ratio - 1.0;
    };

    double kLo = 0.1, kHi = 100.0;
    const double fLo = cvTheorySq(kLo);
    const double fHi = cvTheorySq(kHi);
    if (!std::isfinite(fLo) || !std::isfinite(fHi)) return std::numeric_limits<double>::quiet_NaN();
    if (cvSq <= fHi) return kHi;  // very homogeneous sample -> large k
    if (cvSq >= fLo) return kLo;  // very heterogeneous sample -> small k
    for (int step = 0; step < 50; ++step)
    {
        const double kMid = 0.5 * (kLo + kHi);
        const double fMid = cvTheorySq(kMid);
        if (!std::isfinite(fMid)) break;
        if (std::abs(fMid - cvSq) < 1e-12) return kMid;
        if (fMid > cvSq)
            kLo = kMid;
        else
            kHi = kMid;
    }
    return 0.5 * (kLo + kHi);
}

static inline void doubleBridgeKick(std::vector<int>& tour, aco::Rng& rng)
{
    const int N = (int)tour.size();
    if (N < 8) return;
    if (tour[0] != 0) return; // LS requires node[0]==0

    for (int tries = 0; tries < 50; ++tries)
    {
        // cut points in [1..N-2]
        int a = rng.uniformInt(1, N - 2);
        int b = rng.uniformInt(1, N - 2);
        int c = rng.uniformInt(1, N - 2);
        int d = rng.uniformInt(1, N - 2);

        std::array<int,4> cut{a,b,c,d};
        std::sort(cut.begin(), cut.end());
        a = cut[0]; 
        b = cut[1]; 
        c = cut[2]; 
        d = cut[3];

        // non-empty segments
        if (b <= a + 1) continue;
        if (c <= b + 1) continue;
        if (d <= c + 1) continue;
        if (d >= N - 1) continue;

        std::vector<int> out;
        out.reserve((std::size_t)N);
        out.push_back(0);

        auto append = [&](int l, int rInclusive){
            for (int i = l; i <= rInclusive; ++i)
                out.push_back(tour[(std::size_t)i]);
        };

        // S1=[1..a], S2=[a+1..b], S3=[b+1..c], S4=[c+1..d], S5=[d+1..N-1]
        // 0 + S1 + S4 + S3 + S2 + S5
        append(1, a);
        append(c + 1, d);      // S4
        append(b + 1, c);      // S3
        append(a + 1, b);      // S2
        append(d + 1, N - 1);  // S5

        tour.swap(out);
        return;
    }
}

static bool preferThreeOptPolish(const cars_tsplib::Instance& inst)
{
    const bool isATSP = (inst.type == cars_tsplib::ProblemType::ATSP);

    const bool ewAsym =
        inst.edgeWeightIsSymmetric.has_value() && !*inst.edgeWeightIsSymmetric;

    const bool rrAsym =
        inst.returnRateIsAsymmetric.has_value() && *inst.returnRateIsAsymmetric;

    return isATSP || ewAsym || rrAsym;
}


Colony::Colony(std::shared_ptr<const cars_tsplib::Instance> inst,
               const aco_cli::ParamsData& params,
               int runIndex,
               uint64_t seed,
               std::shared_ptr<const aco::IFixedTourCarAssigner> dp,
               std::shared_ptr<const IAntPolicy> antPolicy,
               std::shared_ptr<IPheromoneModel> pherModel,
               std::shared_ptr<const ICostModel> costModel,
               std::shared_ptr<const ILocalSearch> localSearch,
               std::shared_ptr<const ILocalSearch> lkLite,
               std::shared_ptr<const ILocalSearch> threeOptLite,
               std::shared_ptr<const IIntensifier> intensifier)
    : inst_(std::move(inst))
    , params_(params)
    , runIndex_(runIndex)
    , dp_(dp)
    , seed_(seed)
    , rng_(seed)
    , antPolicy_(std::move(antPolicy))
    , pherModel_(std::move(pherModel))
    , costModel_(std::move(costModel))
    , localSearch_(std::move(localSearch))
    , lkLite_(std::move(lkLite))  
    , threeOptLite_(std::move(threeOptLite))
    , intensifier_(std::move(intensifier))
{
    // Auto-create the intensifier if none was passed in, it's enabled in params, and LS+DP
    // are available; otherwise it stays nullptr (so an explicit disable via -intf 0 is respected).
    if (!intensifier_ && params_.intensifierEnabled && localSearch_ && dp_)
    {
        intensifier_ = std::make_shared<Intensifier>(localSearch_, dp_, inst_);
    }

    // Precompute min travel cost cache for surrogate evaluation
    precomputeMinTravel_();
}

RunResult Colony::run()
{
    RunResult rr;

    auto ensureCars = [&](aco::Solution& s, double& cost)
    {
        if (!dp_) return;
        const int N = inst_->n();
        if ((int)s.car.size() != N) s.car.assign((size_t)N, 0);
        const double c = dp_->reassignCars(s.node, s.car);
        cost = std::isfinite(c) ? c : std::numeric_limits<double>::infinity();
    };

    auto ensureCostOnly = [&](const aco::Solution& s) -> double
    {

        if (dp_) return dp_->evaluateCost(s.node);
        if (costModel_) return costModel_->evaluate(s);
        return std::numeric_limits<double>::infinity();
    };


    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();

    rr.runIndex = runIndex_;
    rr.seed = seed_;
    rr.best.cost = std::numeric_limits<double>::infinity();
    rr.iterBestCost.reserve((size_t)params_.iterations);

    // Tracks whether best was reset to infinity (for archived-best checks)
    bool bestWasResetToInfinity = false;

    AntContext ctx;
    ctx.inst   = inst_;
    ctx.params = &params_;
    ctx.pher   = pherModel_.get();
    ctx.rng    = &rng_;
    ctx.adaptiveQ0 = -1.0; // -1.0 means "use fixed q0"

    if (recorder_) recorder_->setSeed(runIndex_, seed_);
    auto vprint = [&](const std::string& msg)
    {
        if (!verbose_) return;
        if (logger_) 
        {
            logger_->log(msg);
        }
        else if (printMx_)
        {
            std::lock_guard<std::mutex> mlock(*printMx_);
            std::cout << msg << "\n";
        }
        else
        {
            std::cout << msg << "\n";
        }
    };


    const int K = std::clamp(params_.eliteKAnts, 1, params_.antsN);

    // Stagnation counter for kick/LK/3-opt (independent of pheromone_model).
    // Reset after kick/LK/3-opt activates so they don't repeat every iteration.
    int itersSinceGlobalImprovement = 0;

    // Waiting-time Weibull: absolute best cost ever seen (never reset to infinity)
    double allTimeBestCost = std::numeric_limits<double>::infinity();

    const double improveEps = 1e-9;
    // Number of solutions to locally improve per iteration.
    // 0 = dynamic (starts at 2, grows only during stagnation); >0 = fixed count (default: 2)
    const int LSW_FIXED = (params_.lsTopW > 0) ? params_.lsTopW : 0;
    // Current LSW for dynamic mode (starts at 2, grows during stagnation, resets to 2 on exit)
    int currentLSW = 2;


    auto tryInsertElite = [&](std::vector<EvaluatedSolution >& elite, Solution&& s, double cost)
    {
        if (!std::isfinite(cost)) return;

        if ((int)elite.size() < K) 
        {
            elite.push_back({std::move(s), cost});
            return;
        }
        // find the worst in elite and replace it if this one is better
        int worstIdx = 0;
        double worstCost = elite[0].cost;
        for (int i = 1; i < K; ++i) 
        {
            if (elite[i].cost > worstCost) 
            { 
                worstCost = elite[i].cost; worstIdx = i; 
            }
        }
        if (cost < worstCost) 
        {
            elite[worstIdx] = {std::move(s), cost};
        }
    };

    int stagnationEvents = 0;

    // Intensifier research logging (binomial analysis): cost at activation start, at deactivation, improved
    int prevActivationCount = 0;
    double lastActivationStartCost = std::numeric_limits<double>::infinity();

    for (int it = 0; it < params_.iterations; ++it)
    {
        // Report progress (tracked by runner.cpp)
        if (progressCallback_) {
            progressCallback_(it);
        }

        // Dynamic q0: linear ramp over iterations.
        // Early iterations favor exploration (lower q0); later ones favor exploitation (higher q0).
        double q0Value = params_.policyQ0; // default: fixed q0
        if (params_.adaptiveQ0) {
            const double progress = (double)it / (double)params_.iterations;
            const double q0Start = params_.q0Start;
            const double q0End = params_.q0End;
            ctx.adaptiveQ0 = q0Start + (q0End - q0Start) * progress;
            ctx.adaptiveQ0 = std::clamp(ctx.adaptiveQ0, 0.0, 1.0);
            q0Value = ctx.adaptiveQ0;
        } else {
            ctx.adaptiveQ0 = -1.0; // use fixed q0
        }

        // Dynamic LSW: grows during stagnation (more ants get LS-improved) and is restored to
        // the baseline on global improvement. currentLSW carries over from the previous
        // iteration and is updated at the end of this one.
        int LSW;
        if (LSW_FIXED > 0) {
            LSW = LSW_FIXED;
        } else {
            LSW = currentLSW;
        }

        // Pull noImprove_ from pheromone_model for checks (value from before this iteration's update)
        int noImprove = 0;
        if (pherModel_) {
            noImprove = pherModel_->getStagnationCounter();
        }
        
        std::vector<EvaluatedSolution > elite;
        elite.reserve((size_t)K);

        if (params_.ablationIntensifierOnlyNoAnts)
        {
            // Ablation: no ants — one candidate tour per iteration (iter 0: initial tour; later: copy of global best)
            const int Nn = inst_->n();
            auto buildInitialTour = [&](bool randomTour) -> EvaluatedSolution
            {
                EvaluatedSolution out;
                out.sol.node.resize((size_t)Nn);
                if (randomTour)
                {
                    out.sol.node[0] = 0;
                    for (int i = 1; i < Nn; ++i)
                        out.sol.node[(size_t)i] = i;
                    rng_.shuffle(out.sol.node.begin() + 1, out.sol.node.end());
                }
                else
                {
                    for (int i = 0; i < Nn; ++i)
                        out.sol.node[(size_t)i] = i;
                }
                out.sol.car.assign((size_t)Nn, 0);
                if (dp_)
                {
                    out.cost = dp_->evaluateCost(out.sol.node);
                    if (std::isfinite(out.cost))
                    {
                        const double c2 = dp_->reassignCars(out.sol.node, out.sol.car);
                        if (std::isfinite(c2))
                            out.cost = c2;
                    }
                }
                else if (costModel_)
                    out.cost = costModel_->evaluate(out.sol);
                else
                    out.cost = 0.0;
                return out;
            };

            EvaluatedSolution es;
            if (it == 0)
            {
                constexpr int kMaxAttempts = 48;
                for (int att = 0; att < kMaxAttempts; ++att)
                {
                    es = buildInitialTour(params_.ablationIntensifierOnlyRandomInitialTour);
                    if (std::isfinite(es.cost))
                        break;
                    // next attempt: always shuffle (even if -airit 0 made the first attempt the identity tour)
                    es = buildInitialTour(true);
                    if (std::isfinite(es.cost))
                        break;
                }
                if (!std::isfinite(es.cost))
                    es = buildInitialTour(false);
                vprint("run: " + std::to_string(runIndex_) + ", ablation intensifier-only (-aion): initial tour, cost="
                       + std::to_string(es.cost));
            }
            else
            {
                es.sol = rr.best.sol;
                es.cost = rr.best.cost;
                // After archive+reset, cost can be inf even though the tour is still in rr.best.sol — LS expects a consistent cost
                if ((int)es.sol.node.size() == Nn && dp_ && !std::isfinite(es.cost))
                {
                    if ((int)es.sol.car.size() != Nn)
                        es.sol.car.assign((size_t)Nn, 0);
                    es.cost = dp_->evaluateCost(es.sol.node);
                    if (std::isfinite(es.cost))
                    {
                        const double c2 = dp_->reassignCars(es.sol.node, es.sol.car);
                        if (std::isfinite(c2))
                            es.cost = c2;
                    }
                }
                if ((int)es.sol.node.size() != Nn || !std::isfinite(es.cost))
                {
                    es = buildInitialTour(true);
                    if (!std::isfinite(es.cost))
                        es = buildInitialTour(false);
                }
            }
            elite.push_back(std::move(es));
        }
        else
        {
        // 1) Construction + two-phase ranking:
        //    - cheap surrogate cost for all ants
        //    - DP evaluateCost only for the shortlisted candidates
        struct Candidate
        {
            aco::Solution sol;
            double surrogate = std::numeric_limits<double>::infinity();
        };

        // Two-phase only makes sense if we have both dp_ and costModel_; otherwise fall back to the old behavior.
        const bool useTwoPhase = (dp_ != nullptr) && (costModel_ != nullptr);

        if (!useTwoPhase)
        {
            for (int a = 0; a < params_.antsN; ++a)
            {
                Solution s = antPolicy_->construct(ctx);
                double cost = 0.0;
                if (dp_) cost = dp_->evaluateCost(s.node);
                else if (costModel_) cost = costModel_->evaluate(s);
                else cost = 0.0;
                tryInsertElite(elite, std::move(s), cost);
            }
        }
        else // two-phase construction path
        {
            std::vector<Candidate> cands;
            cands.reserve((size_t)params_.antsN);

            // (A) Construct + surrogate cost for all ants.
            // Uses min travel cost per edge (cheaper, and a better predictor of DP cost).
            const int N = inst_->n();
            const int C = inst_->cars();
            for (int a = 0; a < params_.antsN; ++a)
            {
                Candidate c;
                c.sol = antPolicy_->construct(ctx);

                // Surrogate: min travel cost per edge (ignores car switching and return cost).
                // Cheaper and a better DP-cost predictor since it doesn't depend on a poor car[] assignment.
                // Uses the precomputed minTravelCache_ instead of recomputing per loop iteration.
                double surrogate = 0.0;
                for (int i = 0; i < N; ++i)
                {
                    const int u = c.sol.node[i];
                    const int v = c.sol.node[(i + 1) % N];
                    const double minTravel = minTravelCache_[u][v];
                    if (std::isfinite(minTravel))
                        surrogate += minTravel;
                    else
                    {
                        surrogate = std::numeric_limits<double>::infinity();
                        break;
                    }
                }
                // Surrogate with estimated return costs (improves correlation on larger instances)
                if (params_.surrogateReturnGamma > 0 && avgReturnCost_ > 0)
                {
                    const int estReturns = std::max(1, N / 5);  // heuristic estimate of number of returns
                    surrogate += params_.surrogateReturnGamma * static_cast<double>(estReturns) * avgReturnCost_;
                }
                c.surrogate = surrogate;
                cands.push_back(std::move(c));
            }

            // Research logging: surrogate vs full DP cost (for correlation analysis)
            if (researchSurrogateLogger_ && dp_)
            {
                for (size_t i = 0; i < cands.size(); ++i)
                {
                    const double dpCost = dp_->evaluateCost(cands[i].sol.node);
                    researchSurrogateLogger_->log(runIndex_, it, static_cast<int>(i),
                                                  cands[i].surrogate, dpCost);
                }
            }

            // (B) Shortlist size heuristic: typically 4*K, at least 20, never more than antsN.
            const int shortlistM = std::clamp(std::max(20, 4 * K), K, params_.antsN);

            // (C) Wildcard count, so the surrogate doesn't kill diversity
            const int wildcards = std::clamp(
                int(0.1 * params_.antsN),
                2,  // minimum 2
                std::min(30, std::max(2, params_.antsN / 10))  // max 30 or 10% of antsN
            );

            // (D) Select top-shortlistM by surrogate (O(n))
            auto compSur = [](const Candidate& a, const Candidate& b) { return a.surrogate < b.surrogate; };
            if ((int)cands.size() > shortlistM)
                std::nth_element(cands.begin(), cands.begin() + shortlistM, cands.end(), compSur);

            const int baseCount = std::min(shortlistM, (int)cands.size());

            std::vector<int> picked;
            picked.reserve((size_t)(baseCount + wildcards + 16));

            std::vector<uint8_t> pickedFlag(cands.size(), 0);
            auto pickIdx = [&](int idx)
            {
                if (idx < 0 || idx >= (int)cands.size()) return;
                if (pickedFlag[(size_t)idx]) return;
                pickedFlag[(size_t)idx] = 1;
                picked.push_back(idx);
            };

            // (E) DP evaluation only for the shortlist
            for (int i = 0; i < baseCount; ++i)
                pickIdx(i);


            // (F) Wildcards from the remainder (random)
            for (int w = 0; w < wildcards; ++w)
            {
                const int span = (int)cands.size() - baseCount;
                if (span <= 0) break;

                bool pickedOne = false;
                for (int tries = 0; tries < 30 && !pickedOne; ++tries)
                {
                    const int r = baseCount + rng_.uniformInt(span);
                    if (pickedFlag[(size_t)r]) continue;
                    pickIdx(r);
                    pickedOne = true;
                }
            }

            // (G) If elite isn't full yet (e.g. INFs), extend DP evaluation to the rest by surrogate order
            if ((int)picked.size() < K)
            {
                std::sort(cands.begin() + baseCount, cands.end(), compSur);
                for (int i = baseCount; i < (int)cands.size() && (int)picked.size() < K; ++i)
                    pickIdx(i);
            }

            for (int idx : picked)
            {
                // DP exact cost
                const double cost = dp_->evaluateCost(cands[idx].sol.node);
                tryInsertElite(elite, std::move(cands[idx].sol), cost);
            }


        }

        } // !ablationIntensifierOnlyNoAnts


        // If elite is somehow empty (should not happen)
        if (elite.empty()) 
        {
            const auto now = clock::now();
            const double elapsedMs =
            std::chrono::duration<double, std::milli>(now - t0).count();
            rr.iterBestCost.push_back(rr.best.cost);
            if (recorder_)
                recorder_->recordIteration2(runIndex_, it, rr.best.cost, rr.best.cost, rr.best.cost, elapsedMs);
            continue;
        }

        // 2) DP reassign for top-K (fill car[]) is intentionally skipped here: elite[i].cost is
        // already DP-optimal from construction (dp_->evaluateCost(s.node)), and LS moves are
        // evaluated through dp_->evaluateCostView(...), so no "initial" car[] is needed.

        // 2.5) sort elite before LS so LS only runs on the best
        std::sort(elite.begin(), elite.end(),
            [](const EvaluatedSolution & a, const EvaluatedSolution & b) { return a.cost < b.cost; });

        // 3) Local search on the top-LSW elite ants (if LS is configured).
        // LSW is dynamic (starts at 2, grows only during stagnation) or fixed if lsTopW > 0.
        // LS runs every iteration and performs an expensive chain (2opt->reloc->2opt).
        if (localSearch_) 
        {
            const int limit = std::min(LSW, (int)elite.size());
            for (int i = 0; i < limit; ++i) 
            {
                try 
                {
                    localSearch_->improve(elite[i].sol, elite[i].cost);
                    // (LS already calls reassignCars at the end, but double-check anyway)
                    if (dp_ && std::isfinite(elite[i].cost))
                    {
                        const int N = inst_->n();
                        if ((int)elite[i].sol.car.size() != N)
                            elite[i].sol.car.assign((size_t)N, 0);
                        const double newCost = dp_->reassignCars(elite[i].sol.node, elite[i].sol.car);
                        if (std::isfinite(newCost))
                            elite[i].cost = newCost;
                    }
                } catch (const std::exception& ex) 
                {
                    std::cerr << "[LS] " << ex.what() << "\n";
                }
            }
        }
        // 4) Rank elite after LS
        std::sort(elite.begin(), elite.end(),
                [](const EvaluatedSolution & a, const EvaluatedSolution & b) { return a.cost < b.cost; });

        const int W = std::clamp(params_.eliteKAnts, 1, (int)elite.size());

        std::vector<aco::EvaluatedSolution> ranked;
        ranked.reserve((size_t)W);

        for (int r = 0; r < W; ++r)
        {
            aco::EvaluatedSolution es;
            es.cost = elite[r].cost;
            es.sol  = std::move(elite[r].sol);
            ranked.push_back(std::move(es));
        }

        // 4.1) Final consistency: ensure a DP-optimal car[] and cost for the top-W solutions.
        // Cheap since W is small, and it prevents pheromones from learning off inconsistent car assignments.
        if (dp_)
        {
            const int N = inst_->n();
            for (auto& es : ranked)
            {
                if ((int)es.sol.car.size() != N)
                    es.sol.car.assign((size_t)N, 0);

                const double newCost = dp_->reassignCars(es.sol.node, es.sol.car);
                es.cost = std::isfinite(newCost) ? newCost : std::numeric_limits<double>::infinity();
            }

            // cost can change after reassignCars, so re-sort ranked
            std::sort(ranked.begin(), ranked.end(),
                      [](const aco::EvaluatedSolution& a, const aco::EvaluatedSolution& b) 
                      {
                          return a.cost < b.cost;
                      });
        }

        // iterBest is now truly ranked[0] after DP consolidation
        const double iterBestCost = ranked[0].cost;
        rr.iterBestCost.push_back(iterBestCost);

        // remember best before this iteration (for the final globalImproved check)
        const double prevBestCost = rr.best.cost;

        // per-iteration improvement (EPS)
        const bool improvedByIter = (iterBestCost < rr.best.cost - improveEps);

        // update rr.best immediately (before stagnation handling), but as a COPY so ranked stays intact for pheromones
        if (improvedByIter)
            rr.best = ranked[0];   // copy, not move

        // update the local stagnation counter for kick/LK/3-opt
        if (improvedByIter) {
            itersSinceGlobalImprovement = 0;
        } else {
            ++itersSinceGlobalImprovement;
        }

        // iterBestPost is this iteration's final/best cost after all actions (construction, LS on
        // top-LSW elite ants, intensification, polish, kick) — starts at iterBestCost and may
        // improve further. It is the per-iteration best, NOT the cumulative minimum (that's
        // rr.best.cost), and is printed/recorded as "it best" but is never used to update
        // pheromones (those use the ranked top-W solutions instead).
        double iterBestPost = iterBestCost;

        // Remember the cost from intensification improvement (if any) for later use
        double intensificationCost = std::numeric_limits<double>::infinity();

        // 4.5) Intensification: an extra LS cycle applied only to the global best solution.
        // Activates when best was just improved (improvedByIter), or during the first few
        // iterations after stagnation, allowing cheap deeper exploitation since it only touches
        // the best solution. intensify() checks the activation conditions and applies LS + DP to
        // best; on improvement it updates rr.best directly (pass-by-reference) and Colony just
        // checks whether the cost changed.
        bool intensifierActiveThisIteration = false;
        bool disablePheromoneLearningOnActivationThisIteration = false;
        if (intensifier_ && std::isfinite(rr.best.cost))
        {
            intensifierActiveThisIteration =
                improvedByIter ||
                (itersSinceGlobalImprovement > 0 &&
                 itersSinceGlobalImprovement <= params_.intensifierMaxStagnationIterations);
            // remember cost before intensification (to check whether it improved)
            const double costBefore = rr.best.cost;

            // call the intensifier - updates rr.best directly on improvement
            const double newCost = intensifier_->intensify(rr.best, improvedByIter,
                                                          itersSinceGlobalImprovement, improveEps);

            // intensify() already updated rr.best if it improved; just detect whether it happened
            // (check rr.best.cost since it's pass-by-reference and may have been updated)
            if (std::isfinite(newCost) && rr.best.cost < costBefore - improveEps)
            {
                intensificationCost = rr.best.cost;  // remember cost for tracking (iterBestPost)

                // this counts as an improvement, so reset stagnation
                if (itersSinceGlobalImprovement > 0) {
                    itersSinceGlobalImprovement = 0;
                }
            }

            // Research logging: detect the start of a new activation (counter increased)
            if (intensifier_->getActivationCount() > prevActivationCount)
            {
                if (params_.ablationResetPheromonesOnIntensifierActivation)
                {
                    disablePheromoneLearningOnActivationThisIteration = true;
                    if (pherModel_)
                    {
                        pherModel_->reset();
                        vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) +
                               ", ablation: pheromones reset on intensifier activation");
                    }
                }
                lastActivationStartCost = costBefore;
                prevActivationCount = intensifier_->getActivationCount();
                if (researchLogger_)
                    researchLogger_->resetActivationDepth(runIndex_);
            }

            // Research logging: on deactivation, record (run, activation_index, cost_start, cost_end, improved)
            if (intensifier_->wasJustDeactivated() && researchLogger_)
            {
                const int activationCount = intensifier_->getActivationCount();
                const bool improvedThisActivation = (std::isfinite(rr.best.cost) && std::isfinite(lastActivationStartCost) &&
                                                    rr.best.cost < lastActivationStartCost - improveEps);
                researchLogger_->onDeactivation(runIndex_, activationCount, lastActivationStartCost, rr.best.cost, improvedThisActivation);
            }
            
            // Log every time the intensifier turns off (was active, now isn't)
            if (intensifier_->wasJustDeactivated())
            {
                const int activationCount = intensifier_->getActivationCount();
                if (recorder_) {
                    recorder_->recordIntensifierDeactivationEvent(runIndex_, it, activationCount, rr.best.cost);
                }
                vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                       ", intensifier deactivated (activation count: " + std::to_string(activationCount) + ")");
            }
            
            // Pheromone reset after the first deactivation, but only if intensification did NOT
            // improve the best: archive the global best, reset it to infinity, then selectively
            // reduce pheromones on its tour (see steps below).
            if (intensifier_->shouldResetPheromones() && pherModel_)
            {
                // If best wasn't improved by intensification, reset pheromones (selective + smoothing)
                const bool bestNotImproved = (rr.best.cost >= costBefore - improveEps);

                if (bestNotImproved)
                {
                    // 1. Archive the global best (if enabled), before resetting it.
                    // Uses rr.best.sol and costBefore since rr.best.cost is about to be reset to infinity.
                    if (params_.intensifierArchiveAndResetBest && std::isfinite(costBefore) &&
                        (int)rr.best.sol.node.size() == inst_->n())
                    {
                        EvaluatedSolution bestToArchive = rr.best;
                        bestToArchive.cost = costBefore;  // the real cost before resetting
                        pherModel_->addToArchive(bestToArchive);
                        vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) +
                               ", archived global best (cost: " + std::to_string(costBefore) + ")");

                        // Reset global best to infinity (allows fresh exploration)
                        rr.best.cost = std::numeric_limits<double>::infinity();
                        bestWasResetToInfinity = true;
                        vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) +
                               ", reset global best to infinity (allowing fresh exploration)");
                    }

                    // 2. Selective pheromone reset (smoothing on the best tour): instead of a full
                    // reset, reduce pheromones only on the best tour, using gamma as the smoothing
                    // factor (0.0 = full reset, 1.0 = no change).
                    if (std::isfinite(costBefore))
                    {
                        // Use rr.best.sol if still valid (before the cost reset above),
                        // otherwise fall back to ranked[0] if available
                        const Solution* tourToReset = nullptr;
                        if ((int)rr.best.sol.node.size() == inst_->n())
                        {
                            tourToReset = &rr.best.sol;
                        }
                        else if (!ranked.empty() && std::isfinite(ranked[0].cost) && 
                                 (int)ranked[0].sol.node.size() == inst_->n())
                        {
                            tourToReset = &ranked[0].sol;
                        }
                        
                        if (tourToReset)
                        {
                            const double gamma = params_.intensifierPheromoneReductionGamma;
                            pherModel_->reducePheromonesOnTour(*tourToReset, gamma);
                            vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                                   ", selectively reduced pheromones on best tour (gamma: " + std::to_string(gamma) + ")");
                        }
                        else
                        {
                            // Fallback: full reset if we don't have a valid tour
                            pherModel_->reset();
                            vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                                   ", full pheromone reset (no valid tour for selective reset)");
                        }
                    }
                    else
                    {
                        // Fallback: full reset if we don't have costBefore
                        pherModel_->reset();
                        vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                               ", full pheromone reset (no costBefore available)");
                    }
                }
                else
                {
                    // Best WAS improved, do not reset pheromones
                    vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                           ", intensifier deactivated but best improved (no pheromone reset)");
                }
            }
            
            // Update the intensifier's internal state for the next call
            // (after checking deactivation and pheromone reset)
            intensifier_->updateState();
        }
        // When the intensifier is disabled (--intensifier 0), there's no extra LS exclusively on
        // the global best here; LS on elite ants in the main pipeline above still applies.

        // Update iterBestPost if intensification improved this iteration
        if (std::isfinite(intensificationCost) && intensificationCost < iterBestPost)
        {
            iterBestPost = intensificationCost;
        }
        bool skipKick = false;

        // Stagnation handling: LK, 3-opt or kick, based on the local itersSinceGlobalImprovement.
        // --- Kick (double-bridge) on stagnation ---
        if (params_.stagnation > 0 &&
            itersSinceGlobalImprovement >= params_.stagnation &&
            std::isfinite(rr.best.cost))
        {
            ++stagnationEvents;
            if (recorder_) recorder_->incStagnationEvent(runIndex_);

            const bool prefer3 = preferThreeOptPolish(*inst_);

            const aco::ILocalSearch* polish = nullptr;
            uint64_t polishTag = 0;

            enum class PolishKind { None, Lk, ThreeOpt };
            PolishKind kind = PolishKind::None;

            if (prefer3) {
                if (threeOptLite_) { polish = threeOptLite_.get(); polishTag = aco::TAG_THREEOPT; kind = PolishKind::ThreeOpt; }
                else if (lkLite_)  { polish = lkLite_.get();      polishTag = aco::TAG_LKLITE;  kind = PolishKind::Lk; }
            } else {
                if (lkLite_)       { polish = lkLite_.get();      polishTag = aco::TAG_LKLITE;  kind = PolishKind::Lk; }
                else if (threeOptLite_) { polish = threeOptLite_.get(); polishTag = aco::TAG_THREEOPT; kind = PolishKind::ThreeOpt; }
            }

            bool improvedByStag = false;

            // 1) Polish on top elite ants (not just rr.best). Same approach as LS: apply polish to
            // several top elite ants so more than just one elite solution can improve.
            if (polish)
            {
                if (recorder_) 
                {
                    if (kind == PolishKind::ThreeOpt) recorder_->incThreeOptCall(runIndex_);
                    else if (kind == PolishKind::Lk)  recorder_->incLkCall(runIndex_);
                }

                if (kind == PolishKind::ThreeOpt)
                    vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + ", 3-opt call.");
                else if (kind == PolishKind::Lk)
                    vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + ", LK call.");

                const uint64_t seedPolish =
                    aco::mixSeed4(seed_, polishTag, (uint64_t)it, (uint64_t)stagnationEvents);
                polish->resetSeed(seedPolish);

                // Adaptive elite-ant count for polish: top-3 to top-5 based on instance size.
                // Polish is expensive (6 attempts x DP verify per ant), so we focus on the
                // top ants where improvement is most likely.
                const int N = inst_->n();
                int targetPolishCount;
                if (N <= 100) {
                    targetPolishCount = 3;  // small instances: top-3
                } else if (N <= 200) {
                    targetPolishCount = 4;  // medium instances: top-4
                } else {
                    targetPolishCount = 5;  // large instances: top-5
                }
                const int polishLimit = std::min({targetPolishCount, K, (int)ranked.size()});

                for (int i = 0; i < polishLimit; ++i)
                {
                    // Make a copy so we don't mutate the original ranked entry (used for pheromones)
                    aco::EvaluatedSolution cand = ranked[i];

                    try {
                        polish->improve(cand.sol, cand.cost);
                    } catch (const std::exception& ex) {
                        std::cerr << "[Stagnation-Polish] " << ex.what() << "\n";
                        cand.cost = std::numeric_limits<double>::infinity();
                        continue;
                    }

                    if (std::isfinite(cand.cost))
                    {
                        // DP/cost consolidation
                        if (dp_) {
                            const int N = inst_->n();
                            if ((int)cand.sol.car.size() != N) cand.sol.car.assign((size_t)N, 0);
                            const double c2 = dp_->reassignCars(cand.sol.node, cand.sol.car);
                            cand.cost = std::isfinite(c2) ? c2 : std::numeric_limits<double>::infinity();
                        } else if (costModel_) {
                            cand.cost = costModel_->evaluate(cand.sol);
                        }

                        if (std::isfinite(cand.cost))
                        {
                            if (cand.cost < iterBestPost) iterBestPost = cand.cost;

                            // update rr.best on global improvement (check global first)
                            bool globalImproved = false;
                            if (cand.cost < rr.best.cost - improveEps)
                            {
                                rr.best = cand;  // copy, not move (still needed for ranked)
                                improvedByStag = true;
                                globalImproved = true;

                                if (recorder_) 
                                {
                                    if (kind == PolishKind::ThreeOpt) recorder_->incThreeOptSuccess(runIndex_);
                                    else if (kind == PolishKind::Lk)  recorder_->incLkSuccess(runIndex_);
                                }

                                if (kind == PolishKind::ThreeOpt)
                                    vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + ", 3-opt success.");
                                else if (kind == PolishKind::Lk)
                                    vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + ", LK success.");
                            }

                            // update ranked[i] on improvement (local or global)
                            if (cand.cost < ranked[i].cost - improveEps)
                            {
                                ranked[i] = std::move(cand);
                            }
                        }
                    }
                }

                // Sort ranked after all polish improvements (so pheromones learn from sorted solutions)
                if (polishLimit > 0)
                {
                    std::sort(ranked.begin(), ranked.end(),
                        [](const aco::EvaluatedSolution& a, const aco::EvaluatedSolution& b)
                        {
                            return a.cost < b.cost;
                        });
                    // Update iterBestPost after sorting
                    if (ranked[0].cost < iterBestPost) iterBestPost = ranked[0].cost;
                }
            }

            // 2) kick only if polish failed
            // Apply kick to top-2/top-3 ants for greater diversification.
            // Kick is cheaper than polish (perturbation + LS), so it can be applied to more ants.
            if (!improvedByStag)
            {
                // Adaptive ant count for kick: top-2 to top-3 based on instance size
                const int N = inst_->n();
                int targetKickCount;
                if (N <= 100) {
                    targetKickCount = 2;  // small instances: top-2
                } else {
                    targetKickCount = 3;  // large instances: top-3
                }
                const int kickLimit = std::min({targetKickCount, (int)ranked.size()});

                for (int i = 0; i < kickLimit; ++i)
                {
                    if (recorder_ && i == 0) recorder_->incKickCall(runIndex_);
                    if (i == 0) {
                        vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + ", kick call (on " + std::to_string(kickLimit) + " solutions).");
                    }

                    aco::EvaluatedSolution cand = ranked[i];
                    // Generate a distinct seed per ant by mixing in i along with other arguments
                    const uint64_t kickSeed =
                        aco::mixSeed4(seed_, aco::TAG_KICK, (uint64_t)it, (uint64_t)(stagnationEvents * 1000 + i));
                    aco::Rng kickRng(kickSeed);

                    doubleBridgeKick(cand.sol.node, kickRng);

                    if (dp_) cand.cost = dp_->evaluateCost(cand.sol.node);
                    else if (costModel_) cand.cost = costModel_->evaluate(cand.sol);

                    if (localSearch_ && std::isfinite(cand.cost)) {
                        try { localSearch_->improve(cand.sol, cand.cost); }
                        catch (const std::exception& ex) { std::cerr << "[Kick-LS] " << ex.what() << "\n"; }
                    }

                    if (dp_ && std::isfinite(cand.cost)) ensureCars(cand.sol, cand.cost);
                    else if (costModel_ && std::isfinite(cand.cost)) cand.cost = costModel_->evaluate(cand.sol);

                    if (std::isfinite(cand.cost) && cand.cost < iterBestPost)
                        iterBestPost = cand.cost;

                    bool kickImproved = false;
                    if (std::isfinite(cand.cost) && cand.cost < rr.best.cost - improveEps)
                    {
                        rr.best = std::move(cand);
                        kickImproved = true;
                        if (recorder_ && i == 0) recorder_->incKickSuccess(runIndex_);
                        if (i == 0) {
                            vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + ", kick success.");
                        }
                    }

                    // update ranked[i] on improvement
                    if (std::isfinite(cand.cost) && cand.cost < ranked[i].cost - improveEps)
                    {
                        ranked[i] = std::move(cand);
                    }
                }

                // Sort ranked after the kick (so pheromones learn from sorted solutions)
                if (kickLimit > 0)
                {
                    std::sort(ranked.begin(), ranked.end(),
                        [](const aco::EvaluatedSolution& a, const aco::EvaluatedSolution& b)
                        {
                            return a.cost < b.cost;
                        });
                    // Update iterBestPost after sorting
                    if (ranked[0].cost < iterBestPost) iterBestPost = ranked[0].cost;
                }
            }

            // Reset the local stagnation counter after kick/LK/3-opt activation (cooldown).
            // noImprove_ in pheromone_model is left untouched (used for smoothing).
            itersSinceGlobalImprovement = 0;
        }


        // final pheromone update
        // FINAL globalImproved: whether rr.best is better than at the start of the iteration
        // (already checked and logged above when iterBestPost is updated)
        const bool globalImprovedFinal = (rr.best.cost < prevBestCost - improveEps);

        // Waiting-time Weibull: records the true improvement (insensitive to archive+reset)
        const bool trueImprovementThisIter =
            std::isfinite(iterBestPost) && iterBestPost < allTimeBestCost - improveEps;
        if (trueImprovementThisIter) {
            allTimeBestCost = iterBestPost;
            if (recorder_) recorder_->recordGlobalBestImprovement(runIndex_, it);
        }
        
        // Check the archived best after a reset: if best was reset to infinity and a new best was
        // just found, see whether the new best beats the archived one and replace it if so.
        // If it does, reset best to infinity and the pheromones again for continued exploration.
        if (bestWasResetToInfinity && globalImprovedFinal && std::isfinite(rr.best.cost) && pherModel_)
        {
            const double bestArchivedCost = pherModel_->getBestArchivedCost();

            // If the new best beats the archived one, replace it in the archive and reset again
            if (std::isfinite(bestArchivedCost) && rr.best.cost < bestArchivedCost - improveEps)
            {
                // new best is better - add to archive (auto-replaces the old one since archive is sorted)
                pherModel_->addToArchive(rr.best);
                vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) +
                       ", new best (" + std::to_string(rr.best.cost) + ") better than archived (" +
                       std::to_string(bestArchivedCost) + ") - updated archive");

                // reset best to infinity and pheromones again for continued exploration
                if (params_.intensifierArchiveAndResetBest)
                {
                    // remember the tour before resetting, for selective pheromone reset
                    const Solution* tourToReset = nullptr;
                    if ((int)rr.best.sol.node.size() == inst_->n())
                    {
                        tourToReset = &rr.best.sol;
                    }
                    else if (!ranked.empty() && std::isfinite(ranked[0].cost) && 
                             (int)ranked[0].sol.node.size() == inst_->n())
                    {
                        tourToReset = &ranked[0].sol;
                    }
                    
                    // reset best to infinity (allows fully fresh exploration)
                    rr.best.cost = std::numeric_limits<double>::infinity();
                    bestWasResetToInfinity = true;  // keep the flag since we'll reset again
                    vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) +
                           ", reset global best to infinity again (continuous exploration)");

                    // reset pheromones (selectively or fully)
                    if (tourToReset)
                    {
                        const double gamma = params_.intensifierPheromoneReductionGamma;
                        pherModel_->reducePheromonesOnTour(*tourToReset, gamma);
                        vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) +
                               ", selectively reduced pheromones on best tour again (gamma: " + std::to_string(gamma) + ")");
                    }
                    else
                    {
                        // fallback: full reset if there's no valid tour
                        pherModel_->reset();
                        vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) +
                               ", full pheromone reset again (no valid tour for selective reset)");
                    }
                }
            }
            else if (std::isfinite(bestArchivedCost))
            {
                // new best isn't better - archived stays untouched, don't reset
                vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) +
                       ", new best (" + std::to_string(rr.best.cost) + ") not better than archived (" +
                       std::to_string(bestArchivedCost) + ") - keeping archived");

                // reset the flag since we didn't reset (best stays untouched)
                bestWasResetToInfinity = false;
            }
            else
            {
                // no archived best yet (first time), reset the flag
                bestWasResetToInfinity = false;
            }
        }

        // Pheromone update: uses ranked (top-W solutions) and rr.best, NOT iterBestPost.
        // ranked holds the top-W solutions after LS and DP consolidation (before stagnation actions).
        // iterBestPost is the iteration's final cost (after all actions) but isn't used for
        // pheromones, since they learn from the already DP-optimal, sorted top-W solutions —
        // which gives better learning of solution structure.
        bool smoothingActive = false;
        if (pherModel_) {
            // check noImprove_ before the call (to see whether it got reset)
            const int noImproveBefore = pherModel_->getStagnationCounter();

            const bool disablePheromoneLearningDuringIntensification =
                params_.ablationDisableAntDepositsDuringIntensification &&
                intensifier_ &&
                intensifierActiveThisIteration;
            const bool disablePheromoneLearningThisIteration =
                disablePheromoneLearningDuringIntensification ||
                disablePheromoneLearningOnActivationThisIteration ||
                params_.ablationIntensifierOnlyNoAnts;

            if (disablePheromoneLearningThisIteration)
            {
                vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) +
                       ", ablation: skipped pheromone update (no pheromone learning)");
            }
            else
            {
                // don't use iterBestPost here - pheromones learn from top-W solutions, not the iteration's final cost
                pherModel_->onIterationEndRanked(it, ranked, rr.best, globalImprovedFinal);
            }

            // check noImprove_ after the call (to see whether it got reset)
            const int noImproveAfter = pherModel_->getStagnationCounter();

            // log if noImprove_ was reset (decreased)
            if (noImproveBefore > 0 && noImproveAfter == 0) {
                if (globalImprovedFinal) {
                    if (recorder_) {
                        recorder_->recordStagnationResetGlobalImprovementEvent(runIndex_, it, noImproveBefore, rr.best.cost);
                    }
                    if (researchLogger_ && intensifierActiveThisIteration) {
                        researchLogger_->noteGlobalStagnationResetDuringIntensifier(runIndex_);
                    }
                    vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + ", stagnation counter reset: " + std::to_string(noImproveBefore) + " -> 0 (global improvement)");
                } else {
                    vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + ", stagnation counter reset: " + std::to_string(noImproveBefore) + " -> 0 (smoothing activation)");
                }
            }
            
            // check whether smoothing was applied this iteration (from pheromone_model)
            bool smoothingFromModel = pherModel_->wasSmoothingApplied();

            // Smoothing also ramps up as we approach stagnation (same idea as LSW). It activates at
            // 2x the stagnation threshold in pheromone_model; we reuse noImprove_ from there as the
            // single source of truth for the stagnation counter, and start ramping 10 iterations early.
            const int smoothingStagPrag = (params_.stagnation > 0) ? (2 * params_.stagnation) : 0;
            const int noImprovePher = pherModel_->getStagnationCounter();
            if (smoothingStagPrag > 0 && noImprovePher >= smoothingStagPrag - 10) {
                smoothingActive = true;
            } else if (smoothingFromModel) {
                smoothingActive = true;
            }
        }

        // elapsed time at the end of the loop, after stagnation actions (so it includes LK/kick/LS)
        const auto now = clock::now();
        const double elapsedMs =
            std::chrono::duration<double, std::milli>(now - t0).count();

        // Update LSW at the end of the iteration (for the next iteration and for logging).
        // Hybrid approach: LSW grows gradually across all iterations, plus extra growth during stagnation.
        // Note: if LSW_FIXED > 0, LSW was already set to LSW_FIXED earlier, so it isn't changed here.
        if (LSW_FIXED == 0) {
            const int LSW_START = 2;
            const int LSW_END = std::min(K, 7); // raised from 5 to 7 for more LS during stagnation

            // 1) Base LSW: gradual growth across all iterations (80% of max growth, for faster convergence)
            const double progress = (double)it / (double)params_.iterations;
            const double baseLSW = LSW_START + (LSW_END - LSW_START) * 0.8 * progress;

            // 2) Stagnation bonus: extra growth during stagnation (remaining 20% of max growth)
            double stagnationBonus = 0.0;
            if (params_.stagnation > 0 && itersSinceGlobalImprovement >= params_.stagnation) {
                const int stagnationProgress = itersSinceGlobalImprovement - params_.stagnation;
                const int maxStagnationGrowth = 10;
                const double stagnationFactor = std::min(1.0, (double)stagnationProgress / (double)maxStagnationGrowth);
                stagnationBonus = (LSW_END - LSW_START) * 0.2 * stagnationFactor;
            } else if (params_.stagnation > 0 && itersSinceGlobalImprovement >= params_.stagnation - 10) {
                const int approachProgress = itersSinceGlobalImprovement - (params_.stagnation - 10);
                const int approachRange = 10;
                const double approachFactor = (double)approachProgress / (double)approachRange;
                stagnationBonus = (LSW_END - LSW_START) * 0.2 * 0.5 * approachFactor;
            }
            
            // 3) Combine base LSW and stagnation bonus
            double targetLSW = baseLSW + stagnationBonus;
            currentLSW = (int)std::round(targetLSW);
            currentLSW = std::clamp(currentLSW, LSW_START, LSW_END);

            // Update LSW for logging and recording
            LSW = currentLSW;
        }

        // Final iteration cost check, done at the end of the iteration.
        // iterBestPost is the best cost achieved this iteration (after all actions: construction,
        // LS, intensification, polish, kick) — NOT the cumulative minimum (that's rr.best.cost).
        // It's the true/final iteration cost used for: CLI/log output ("it best"), recorder
        // (iterBestPostCost). It is NOT used for the pheromone update (which uses ranked top-W).
        if (!ranked.empty() && std::isfinite(ranked[0].cost))
        {
            iterBestPost = std::min(iterBestPost, ranked[0].cost);
        }

        // If iterBestPost beats rr.best.cost, update rr.best (the whole solution, not just cost).
        // rr.best.cost must remain the cumulative minimum across all iterations. rr.best may
        // already have been updated earlier in the iteration (intensification/polish/kick) or
        // via the iterBestCost check above; here we confirm iterBestPost (this iteration's best,
        // after all actions) beats the current rr.best.cost.
        EvaluatedSolution bestPostSolution;  // solution corresponding to iterBestPost
        bool hasBestPostSolution = false;

        if (std::isfinite(iterBestPost) && iterBestPost < rr.best.cost - improveEps)
        {
            const double oldBestCost = rr.best.cost;
            // If iterBestPost came from ranked[0], use it directly; otherwise rr.best was already
            // updated earlier in the iteration (intensification/polish/kick or the iterBestCost check).
            if (!ranked.empty() && std::isfinite(ranked[0].cost) &&
                std::abs(iterBestPost - ranked[0].cost) < improveEps)
            {
                // iterBestPost comes from ranked[0]
                rr.best = ranked[0];  // copy of the whole solution (sol + cost), not a move
                bestPostSolution = ranked[0];
                hasBestPostSolution = true;
            }
            else
            {
                // iterBestPost == rr.best.cost; rr.best was already updated earlier
                bestPostSolution = rr.best;
                hasBestPostSolution = true;
            }

            // Log every time a new global best (cumulative minimum) is reached
            vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) +
                   ", NEW GLOBAL BEST: " + std::to_string(oldBestCost) +
                   " -> " + std::to_string(rr.best.cost));
        }
        else if (std::isfinite(iterBestPost))
        {
            // iterBestPost doesn't beat rr.best, but remember the solution for the archived-best check
            if (!ranked.empty() && std::isfinite(ranked[0].cost) &&
                std::abs(iterBestPost - ranked[0].cost) < improveEps)
            {
                bestPostSolution = ranked[0];
                hasBestPostSolution = true;
            }
            else if (std::isfinite(rr.best.cost) && std::abs(iterBestPost - rr.best.cost) < improveEps)
            {
                bestPostSolution = rr.best;
                hasBestPostSolution = true;
            }
        }

        // Archived-best check after iterBestPost: if best was reset to infinity and iterBestPost
        // is better than the archived best, this catches the case where iterBestPost beats an
        // rr.best that was archived earlier.
        if (bestWasResetToInfinity && hasBestPostSolution && std::isfinite(bestPostSolution.cost) && pherModel_)
        {
            const double bestArchivedCost = pherModel_->getBestArchivedCost();

            // if best post beats the archived one, replace it in the archive
            if (std::isfinite(bestArchivedCost) && bestPostSolution.cost < bestArchivedCost - improveEps)
            {
                // best post is better - add to archive (auto-replaces the old one since archive is sorted)
                pherModel_->addToArchive(bestPostSolution);
                vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) +
                       ", iterBestPost (" + std::to_string(bestPostSolution.cost) + ") better than archived (" +
                       std::to_string(bestArchivedCost) + ") - updated archive with best post");

                // reset best to infinity and pheromones again for continued exploration
                if (params_.intensifierArchiveAndResetBest)
                {
                    // remember the tour before resetting, for selective pheromone reset
                    const Solution* tourToReset = nullptr;
                    if ((int)bestPostSolution.sol.node.size() == inst_->n())
                    {
                        tourToReset = &bestPostSolution.sol;
                    }
                    else if (!ranked.empty() && std::isfinite(ranked[0].cost) &&
                             (int)ranked[0].sol.node.size() == inst_->n())
                    {
                        tourToReset = &ranked[0].sol;
                    }

                    // reset best to infinity (allows fully fresh exploration)
                    rr.best.cost = std::numeric_limits<double>::infinity();
                    bestWasResetToInfinity = true;  // keep the flag since we'll reset again
                    vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) +
                           ", reset global best to infinity again (best post better than archived)");

                    // reset pheromones (selectively or fully)
                    if (tourToReset)
                    {
                        const double gamma = params_.intensifierPheromoneReductionGamma;
                        pherModel_->reducePheromonesOnTour(*tourToReset, gamma);
                        vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) +
                               ", selectively reduced pheromones on best post tour (gamma: " + std::to_string(gamma) + ")");
                    }
                    else
                    {
                        // fallback: full reset if there's no valid tour
                        pherModel_->reset();
                        vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) +
                               ", full pheromone reset again (no valid tour for selective reset)");
                    }
                }
            }
        }
        
        if (recorder_) {
            // Weibull k from this iteration's ant costs (ranked)
            const double weibullK = weibullKFromCosts(ranked);
            // Recorder fields:
            // - iterBestCost = best cost before stagnation actions (after construction + LS)
            // - iterBestPost = best cost this iteration after all actions (final iteration cost)
            // - rr.best.cost = cumulative minimum (already up to date since it's only updated on improvement)
            recorder_->recordIteration2(runIndex_, it, iterBestCost, iterBestPost, rr.best.cost, elapsedMs, weibullK);
            // record q0, LSW and smoothing for this iteration
            recorder_->recordIterationParams(runIndex_, it, q0Value, LSW, smoothingActive);

            // Diversity stats for the H5/H6 hypotheses (from the same costs as weibullKFromCosts)
            {
                double dMean = 0.0, dStd = 0.0, dRange = 0.0, dCV = 0.0;
                double qMin = std::numeric_limits<double>::quiet_NaN();
                double q1 = std::numeric_limits<double>::quiet_NaN();
                double q2 = std::numeric_limits<double>::quiet_NaN();
                double q3 = std::numeric_limits<double>::quiet_NaN();
                double qMax = std::numeric_limits<double>::quiet_NaN();
                std::vector<double> costs;
                costs.reserve(ranked.size());
                for (const auto& es : ranked)
                    if (std::isfinite(es.cost) && es.cost > 0.0)
                        costs.push_back(es.cost);
                const size_t n = costs.size();
                if (n >= 1) {
                    double sum = 0.0;
                    for (double x : costs) sum += x;
                    dMean = sum / static_cast<double>(n);
                }
                if (n >= 2) {
                    double sumSq = 0.0;
                    for (double x : costs) sumSq += (x - dMean) * (x - dMean);
                    dStd = std::sqrt(sumSq / static_cast<double>(n));
                    dRange = *std::max_element(costs.begin(), costs.end())
                           - *std::min_element(costs.begin(), costs.end());
                    if (dMean > 0.0) dCV = dStd / dMean;
                } else if (n == 1) {
                    dStd = 0.0;
                    dRange = 0.0;
                    dCV = 0.0;
                }
                if (n >= 1) {
                    std::sort(costs.begin(), costs.end());
                    auto quantile = [&](double p) -> double {
                        if (costs.size() == 1) return costs[0];
                        const double pos = (static_cast<double>(costs.size()) - 1.0) * p;
                        const size_t i = static_cast<size_t>(pos);
                        const double frac = pos - static_cast<double>(i);
                        if (i + 1 >= costs.size()) return costs.back();
                        return (1.0 - frac) * costs[i] + frac * costs[i + 1];
                    };
                    qMin = costs.front();
                    q1 = quantile(0.25);
                    q2 = quantile(0.50);
                    q3 = quantile(0.75);
                    qMax = costs.back();
                }
                recorder_->recordDiversity(runIndex_, it, dStd, dRange, dMean, dCV,
                                           trueImprovementThisIter);
                recorder_->recordAntQuantiles(runIndex_, it, qMin, q1, q2, q3, qMax, dMean, dStd,
                                              static_cast<int>(n));
            }
        }
        
        // Format: run, iteration, stagCount (kick/LK), q0, LSW, noImprovCount (smoothing), smoothing, it best, best
        // stagCount = local itersSinceGlobalImprovement (for kick/LK, resets on activation)
        // noImprovCount = noImprove_ from pheromone_model (for smoothing, resets on global improvement or smoothing)
        int noImprovCount = 0;
        if (pherModel_) {
            noImprovCount = pherModel_->getStagnationCounter();
        }
        
        // Compute median/mean global-best for this iteration across all runs
        // (only runs that have already computed their value; slower ones are ignored)
        double medianGlobalBest = StatsRow::nan();
        double meanGlobalBest = StatsRow::nan();
        if (recorder_)
        {
            medianGlobalBest = recorder_->getMedianGlobalBestForIteration(it);
            meanGlobalBest = recorder_->getMeanGlobalBestForIteration(it);
        }
        
        // CLI/log output:
        // - "it best" = iterBestPost (final iteration cost: best cost after all actions)
        // - "best" = rr.best.cost (cumulative minimum across all iterations)
        // - "median"/"mean" = median/mean global-best for this iteration across all runs (computed so far)
        std::ostringstream oss;
        oss << "run: " << runIndex_
            << ", iteration: " << it
            << ", stagCount: " << itersSinceGlobalImprovement
            << ", q0: " << std::fixed << std::setprecision(3) << q0Value
            << ", LSW: " << LSW
            << ", noImprovCount: " << noImprovCount
            << ", smoothing: " << (smoothingActive ? 1 : 0)
            << ", it best: " << std::setprecision(3) << iterBestPost  // final iteration cost
            << ", best: " << std::setprecision(3) << rr.best.cost;     // cumulative minimum
        
        if (std::isfinite(medianGlobalBest))
        {
            oss << ", median: " << std::setprecision(3) << medianGlobalBest;
        }
        
        if (std::isfinite(meanGlobalBest))
        {
            oss << ", mean: " << std::setprecision(3) << meanGlobalBest;
        }
        
        vprint(oss.str());

        // Cache statistics: printed every 50 iterations or on the last iteration,
        // only for the verbose run (the currently selected run).
        if (verbose_ && dp_ && (it % 50 == 0 || it == params_.iterations - 1))
        {
            // check whether dp_ is a CachedFixedTourCarAssigner
            const aco::CachedFixedTourCarAssigner* cached = 
                dynamic_cast<const aco::CachedFixedTourCarAssigner*>(dp_.get());
            
            if (cached)
            {
                const auto stats = cached->stats();
                const uint64_t evalTotal = stats.evalHits + stats.evalMiss;
                const uint64_t reasTotal = stats.reasHits + stats.reasMiss;
                
                double evalHitRate = 0.0;
                double reasHitRate = 0.0;
                
                if (evalTotal > 0)
                    evalHitRate = 100.0 * (double)stats.evalHits / (double)evalTotal;
                if (reasTotal > 0)
                    reasHitRate = 100.0 * (double)stats.reasHits / (double)reasTotal;
                
                std::ostringstream oss;
                oss << "[Cache] run=" << runIndex_ << ", iter=" << it
                    << ", evalHitRate=" << std::fixed << std::setprecision(2) << evalHitRate << "%"
                    << " (" << stats.evalHits << "/" << evalTotal << ")"
                    << ", reasHitRate=" << reasHitRate << "%"
                    << " (" << stats.reasHits << "/" << reasTotal << ")"
                    << ", inserts=" << stats.inserts;
                
                if (logger_)
                {
                    logger_->log(oss.str());
                }
                else if (printMx_)
                {
                    std::lock_guard<std::mutex> mlock(*printMx_);
                    std::cout << oss.str() << "\n";
                }
                else
                {
                    std::cout << oss.str() << "\n";
                }
            }
        }

    }    

    // final consolidation, just in case (global best vs DP)
    if (dp_ && std::isfinite(rr.best.cost))
    {
        const int N = inst_->n();
        try {
            if ((int)rr.best.sol.node.size() != N || (int)rr.best.sol.car.size() != N) {
                if (verbose_)
                    vprint("run: " + std::to_string(runIndex_) + ", WARN final consolidation: rr.best size mismatch "
                           "(node=" + std::to_string(rr.best.sol.node.size()) + ", car=" + std::to_string(rr.best.sol.car.size()) + ", expected N=" + std::to_string(N) + ")");
            } else {
                if ((int)rr.best.sol.car.size() != N) rr.best.sol.car.assign((size_t)N, 0);
                const double c = dp_->reassignCars(rr.best.sol.node, rr.best.sol.car);
                rr.best.cost = std::isfinite(c) ? c : std::numeric_limits<double>::infinity();
            }
        } catch (const std::exception& e) {
            if (verbose_)
                vprint("run: " + std::to_string(runIndex_) + ", WARN final consolidation failed: " + std::string(e.what())
                       + " (node.size=" + std::to_string(rr.best.sol.node.size()) + ", inst->n()=" + std::to_string(N) + ")");
        }
    }
    
    // FINAL ARCHIVED-BEST CHECK (global best vs archived/cache):
    // At the end of the run, check whether the archived best beats the current rr.best, and if
    // so, use it as the final solution. This matters because best can be reset to infinity during
    // the run, but the archived best is preserved and may be better than the final rr.best.
    if (params_.intensifierArchiveAndResetBest && pherModel_)
    {
        EvaluatedSolution archivedBest;
        if (pherModel_->getBestArchivedSolution(archivedBest))
        {
            // Check whether the archived best beats the current rr.best
            // (using the same epsilon as the rest of the code, 1e-12)
            const double improveEps = 1e-12;
            const bool archivedIsBetter = std::isfinite(archivedBest.cost) && 
                                         (rr.best.cost >= archivedBest.cost - improveEps || !std::isfinite(rr.best.cost));
            
            if (archivedIsBetter)
            {
                const int N = inst_->n();
                if ((int)archivedBest.sol.node.size() != N || (int)archivedBest.sol.car.size() != N) {
                    if (verbose_)
                        vprint("run: " + std::to_string(runIndex_) + ", WARN archived best size mismatch - skipping "
                               "(node=" + std::to_string(archivedBest.sol.node.size()) + ", car=" + std::to_string(archivedBest.sol.car.size()) + ", expected N=" + std::to_string(N) + ")");
                } else {
                    // Archived best is better - use it as rr.best
                    rr.best = archivedBest;

                    // Final consolidation for the archived best (ensure a DP-optimal car[])
                    if (dp_ && std::isfinite(rr.best.cost))
                    {
                        try {
                            if ((int)rr.best.sol.car.size() != N) rr.best.sol.car.assign((size_t)N, 0);
                            const double c = dp_->reassignCars(rr.best.sol.node, rr.best.sol.car);
                            rr.best.cost = std::isfinite(c) ? c : std::numeric_limits<double>::infinity();
                            if (verbose_)
                                vprint("run: " + std::to_string(runIndex_) + ", final: using archived best (cost: " + 
                                       std::to_string(rr.best.cost) + ") as final solution");
                        } catch (const std::exception& e) {
                            if (verbose_)
                                vprint("run: " + std::to_string(runIndex_) + ", WARN archived best reassignCars failed: " + std::string(e.what())
                                       + " (node.size=" + std::to_string(rr.best.sol.node.size()) + ") - keeping archived cost");
                        }
                    } else if (verbose_) {
                        vprint("run: " + std::to_string(runIndex_) + ", final: using archived best (cost: " + 
                               std::to_string(rr.best.cost) + ") as final solution");
                    }
                }
            }
        }
    }

    // Final cache statistics at the end of the run.
    // Recorded for all runs (not just the verbose one).
    if (dp_)
    {
        const aco::CachedFixedTourCarAssigner* cached = 
            dynamic_cast<const aco::CachedFixedTourCarAssigner*>(dp_.get());
        
        if (cached)
        {
            const auto stats = cached->stats();
            const uint64_t evalTotal = stats.evalHits + stats.evalMiss;
            const uint64_t reasTotal = stats.reasHits + stats.reasMiss;
            
            double evalHitRate = 0.0;
            double reasHitRate = 0.0;
            
            if (evalTotal > 0)
                evalHitRate = 100.0 * (double)stats.evalHits / (double)evalTotal;
            if (reasTotal > 0)
                reasHitRate = 100.0 * (double)stats.reasHits / (double)reasTotal;
            
            // record into the recorder for CSV output
            if (recorder_)
            {
                recorder_->setCacheStats(runIndex_, evalHitRate, reasHitRate);
            }

            // only print for the verbose run
            if (verbose_)
            {
                std::ostringstream oss;
                oss << "[Cache FINAL] run=" << runIndex_
                    << ", evalHitRate=" << std::fixed << std::setprecision(2) << evalHitRate << "%"
                    << " (" << stats.evalHits << "/" << evalTotal << ")"
                    << ", reasHitRate=" << reasHitRate << "%"
                    << " (" << stats.reasHits << "/" << reasTotal << ")"
                    << ", inserts=" << stats.inserts;
                
                if (logger_)
                {
                    logger_->log(oss.str());
                }
                else if (printMx_)
                {
                    std::lock_guard<std::mutex> mlock(*printMx_);
                    std::cout << oss.str() << "\n";
                }
                else
                {
                    std::cout << oss.str() << "\n";
                }
            }
        }
    }

    return rr;
}

void Colony::precomputeMinTravel_()
{
    const int N = inst_->n();
    const int C = inst_->cars();
    
    if (N <= 0 || C <= 0) return;
    
    minTravelCache_.assign(N, std::vector<double>(N, std::numeric_limits<double>::infinity()));
    
    for (int u = 0; u < N; ++u)
    {
        for (int v = 0; v < N; ++v)
        {
            if (u == v) continue; // skip diagonal
            
            double minT = std::numeric_limits<double>::infinity();
            for (int c = 0; c < C; ++c)
            {
                const double tc = inst_->travelCost(c, u, v);
                if (std::isfinite(tc))
                    minT = std::min(minT, tc);
            }
            minTravelCache_[u][v] = minT;
        }
    }

    // Average return cost (for the surrogate's return term — improves correlation on larger instances)
    avgReturnCost_ = 0.0;
    if (inst_->hasReturnCosts())
    {
        double sumRet = 0.0;
        int countRet = 0;
        for (int c = 0; c < C; ++c)
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    if (i != j)
                    {
                        const double rc = inst_->returnCost(c, i, j);
                        if (std::isfinite(rc))
                        {
                            sumRet += rc;
                            ++countRet;
                        }
                    }
        if (countRet > 0)
            avgReturnCost_ = sumRet / countRet;
    }
}

} // namespace aco
