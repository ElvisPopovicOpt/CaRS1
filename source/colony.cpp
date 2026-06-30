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

// Procjena Weibull shape parametra k iz costova mrava (moment method: CV = sigma/mu).
// Vraća NaN ako premalo uzoraka ili nema pozitivnih costova.
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

    // Weibull: CV^2 = Gamma(1+2/k)/Gamma(1+1/k)^2 - 1. Tražimo k.
    // cv_theory_sq(k) = Gamma(1+2/k)/Gamma(1+1/k)^2 - 1 je monotono padajuća u k.
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
    if (cvSq <= fHi) return kHi;  // uzorak vrlo homogen -> veliki k
    if (cvSq >= fLo) return kLo;  // uzorak vrlo heterogen -> mali k
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
    if (tour[0] != 0) return; // jer tvoj LS zahtijeva node<a href="" class="citation-link" target="_blank" style="vertical-align: super; font-size: 0.8em; margin-left: 3px;">[0]</a>==0

    for (int tries = 0; tries < 50; ++tries)
    {
        // rezovi u [1..N-2]
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

        // ne-prazni segmenti
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
    // Ovo prilagodi tvojim stvarnim getterima/fieldovima:
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
    // AUTOMATSKO KREIRANJE INTENZIFIKATORA:
    // Ako intenzifikator nije proslijeđen i u parametrima je uključen, kreiraj ga ako imamo LS+DP.
    // Ako je params_.intensifierEnabled == false (npr. -intf 0), nullptr ostaje nullptr — nema
    // auto-kreiranja, inače bi se intenzifikator palio protiv eksplicitnog isključivanja.
    if (!intensifier_ && params_.intensifierEnabled && localSearch_ && dp_)
    {
        intensifier_ = std::make_shared<Intensifier>(localSearch_, dp_, inst_);
    }
    
    // Precompute min travel cost cache za surrogate evaluaciju
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
    
    // Pratimo je li best resetiran na infinity (za provjeru arhiviranog besta)
    bool bestWasResetToInfinity = false;

    AntContext ctx;
    ctx.inst   = inst_;
    ctx.params = &params_;
    ctx.pher   = pherModel_.get();
    ctx.rng    = &rng_;
    ctx.adaptiveQ0 = -1.0; // inicijaliziraj na -1.0 (koristi fiksni q0)
    ctx.adaptiveQ0 = -1.0; // inicijaliziraj na -1.0 (koristi fiksni q0)

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


    const int K = std::clamp(params_.eliteKAnts, 1, params_.antsN); // dodaj u params; ili hardcode 5

    // stagnacija brojač za kick/LK/3-opt (neovisno od pheromone_model)
    // Resetuje se nakon aktivacije kick/LK/3-opt da se ne ponavljaju u svakoj iteraciji
    int itersSinceGlobalImprovement = 0;

    // Waiting-time Weibull: apsolutni minimum costa (nikad se ne resetira na infinity)
    double allTimeBestCost = std::numeric_limits<double>::infinity();
    
    // EPS
    const double improveEps = 1e-9;
    // koliko rješenja lokalno popravljamo
    // 0 = auto (dinamički: počinje sa 2, raste samo tokom stagnacije)
    // >0 = fiksni broj (default: 2)
    const int LSW_FIXED = (params_.lsTopW > 0) ? params_.lsTopW : 0; // 0 = dinamički, >0 = fiksni
    // Trenutni LSW za dinamički mod (počinje sa 2, raste tokom stagnacije, vraća se na 2 kada se izađe)
    int currentLSW = 2;


    auto tryInsertElite = [&](std::vector<EvaluatedSolution >& elite, Solution&& s, double cost)
    {
        if (!std::isfinite(cost)) return;

        if ((int)elite.size() < K) 
        {
            elite.push_back({std::move(s), cost});
            return;
        }
        // nađi najgoreg u elite i zamijeni ako je bolji
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

    int stagnationEvents = 0; // NOVO

    // Istraživačko logiranje intenzifikatora (binomna analiza): cost na početku aktivacije, na deaktivaciji, improved
    int prevActivationCount = 0;
    double lastActivationStartCost = std::numeric_limits<double>::infinity();

    for (int it = 0; it < params_.iterations; ++it)
    {
        // Ažuriraj napredak (za praćenje u runner.cpp)
        if (progressCallback_) {
            progressCallback_(it);
        }
        
        // Dinamički q0: linearno povećanje tijekom iteracija
        // U ranim iteracijama više eksploracije (niži q0), u kasnim više eksploatacije (viši q0)
        double q0Value = params_.policyQ0; // default: fiksni q0
        if (params_.adaptiveQ0) {
            const double progress = (double)it / (double)params_.iterations;
            const double q0Start = params_.q0Start;
            const double q0End = params_.q0End;
            ctx.adaptiveQ0 = q0Start + (q0End - q0Start) * progress;
            ctx.adaptiveQ0 = std::clamp(ctx.adaptiveQ0, 0.0, 1.0);
            q0Value = ctx.adaptiveQ0; // dinamički q0
        } else {
            ctx.adaptiveQ0 = -1.0; // koristi fiksni q0
        }

        // Dinamički LSW: raste samo tokom stagnacije, vraća se na 2 kada se izađe
        // Kada algoritam stagnira, povećavamo LSW da bi se više mrava poboljšalo LS-om
        // Kada se izađe iz stagnacije (global improvement), vraćamo LSW na početnu vrednost
        // Koristimo currentLSW iz prethodne iteracije (računa se na kraju prethodne iteracije)
        // Na kraju ove iteracije ćemo ažurirati currentLSW za sljedeću iteraciju
        int LSW;
        if (LSW_FIXED > 0) {
            LSW = LSW_FIXED; // fiksni
        } else {
            LSW = currentLSW; // koristimo currentLSW iz prethodne iteracije
        }
        
        // Dobij noImprove_ iz pheromone_model za provjere (vrijednost iz prethodne iteracije, prije ažuriranja)
        int noImprove = 0;
        if (pherModel_) {
            noImprove = pherModel_->getStagnationCounter();
        }
        
        std::vector<EvaluatedSolution > elite;
        elite.reserve((size_t)K);

        if (params_.ablationIntensifierOnlyNoAnts)
        {
            // Ablation: bez mrava — jedna kandidat-tura po iteraciji (iter 0: početna; kasnije: kopija global besta)
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
                    // sljedeći pokušaj: uvijek promiješaj (i ako je -airit 0, prvi pokušaj je identitet)
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
                // Nakon archive+reset besta cost može biti inf iako je tura još u rr.best.sol — LS očekuje konzistentan cost
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
        // 1) konstrukcija + DVOFAZNO rangiranje:
        //    - surrogate za sve mrave (jeftino)
        //    - DP evaluateCost samo za shortlist kandidata
        struct Candidate
        {
            aco::Solution sol;
            double surrogate = std::numeric_limits<double>::infinity();
        };

        // Ako nemamo i dp_ i costModel_, nema smisla raditi dvofazno => fallback na staro ponašanje.
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
        else // two phase else blok
        {
            std::vector<Candidate> cands;
            cands.reserve((size_t)params_.antsN);

            // (A) Construct + surrogate za sve
            // Koristimo min travel cost po bridu (jeftinije i bolji prediktor DP cost-a)
            const int N = inst_->n();
            const int C = inst_->cars();
            for (int a = 0; a < params_.antsN; ++a)
            {
                Candidate c;
                c.sol = antPolicy_->construct(ctx);
                
                // Surrogate: min travel cost po bridu (bez car switching-a i return cost-a)
                // Ovo je jeftinije i bolji prediktor DP cost-a jer ne ovisi o lošem car[] assignmentu
                // OPTIMIZACIJA: Koristi precomputed minTravelCache_ umjesto računanja u petlji
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
                // Surogat s procjenom return costova (poboljšanje korelacije na većim instancama)
                if (params_.surrogateReturnGamma > 0 && avgReturnCost_ > 0)
                {
                    const int estReturns = std::max(1, N / 5);  // heuristička procjena broja returnova
                    surrogate += params_.surrogateReturnGamma * static_cast<double>(estReturns) * avgReturnCost_;
                }
                c.surrogate = surrogate;
                cands.push_back(std::move(c));
            }

            // Istraživačko logiranje: surogat vs puna DP cijena (za analizu korelacije)
            if (researchSurrogateLogger_ && dp_)
            {
                for (size_t i = 0; i < cands.size(); ++i)
                {
                    const double dpCost = dp_->evaluateCost(cands[i].sol.node);
                    researchSurrogateLogger_->log(runIndex_, it, static_cast<int>(i),
                                                  cands[i].surrogate, dpCost);
                }
            }

            // (B) shortlist veličina (bez novih params: heuristika)
            // Tipično 4*K, ali barem 20 i ne više od antsN.
            const int shortlistM = std::clamp(std::max(20, 4 * K), K, params_.antsN);  

            // (C) wildcard count (da surrogate ne "ubije" diverzitet)
            //const int wildcards = std::min(2, std::max(0, params_.antsN - shortlistM));
            const int wildcards = std::clamp(
                int(0.1 * params_.antsN),
                2,  // minimum 2 (umjesto 5 za manje instance)
                std::min(30, std::max(2, params_.antsN / 10))  // max 30 ili 10% od antsN
            );

            // (D) Izaberi top-shortlistM po surrogate (O(n))
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

            // (E) DP evaluacija samo za shortlist
            for (int i = 0; i < baseCount; ++i)
                pickIdx(i);


            // (F) Wildcards iz ostatka (slučajno)
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

            // (G) Ako elite nije puna (INF-ovi ili sl.), proširi DP evaluaciju na preostale po surrogate
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


        // Ako iz nekog razloga nema elite (ne bi se smjelo dogoditi)
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

        // 2) Za top-K: napravi DP reassign (popuni car[]) da pheromoni/LS budu konzistentni
        // Zasto je zakomentirano:
        // elite[i].cost je već DP-optimalan jer u konstrukciji radiš dp_->evaluateCost(s.node).
        // LS poteze evaluira kroz dp_->evaluateCostView(...) (DP-view), znači ne treba “inicijalni” car[].
        /*
        if (dp_) 
        {
            const int N = inst_->n();
            for (auto& e : elite) 
            {
                if ((int)e.sol.car.size() != N)
                    e.sol.car.assign((size_t)N, 0);
                const double newCost = dp_->reassignCars(e.sol.node, e.sol.car);
                e.cost = std::isfinite(newCost) ? newCost : std::numeric_limits<double>::infinity();
            }
        }
        */
        // 2.5) sortiraj elite prije LS da LS ide samo na najbolje
        std::sort(elite.begin(), elite.end(),
            [](const EvaluatedSolution & a, const EvaluatedSolution & b) { return a.cost < b.cost; });

        // 3) Local search na top-LSW elite mrave (ako imaš LS)
        // LSW je dinamički (počinje sa 2, raste samo tokom stagnacije) ili fiksni ako je lsTopW > 0
        // LS se poziva svaku iteraciju i radi skupi chain (2opt->reloc->2opt)
        if (localSearch_) 
        {
            const int limit = std::min(LSW, (int)elite.size());
            for (int i = 0; i < limit; ++i) 
            {
                try 
                {
                    localSearch_->improve(elite[i].sol, elite[i].cost);
                    // (LS poziva reassignCars na kraju, ali provjerimo za svaki slučaj)
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
        // 4) Rankaj elite nakon LS
        std::sort(elite.begin(), elite.end(),
                [](const EvaluatedSolution & a, const EvaluatedSolution & b) { return a.cost < b.cost; });

        const int W = std::clamp(params_.eliteKAnts, 1, (int)elite.size());
        // Ako eliteKAnts zna biti 8–10, a update postane previše “razmazan”, privremeno u Colony staviti
        //const int W = std::min(std::clamp(params_.eliteKAnts, 1, (int)elite.size()), 5);

        std::vector<aco::EvaluatedSolution> ranked;
        ranked.reserve((size_t)W);

        for (int r = 0; r < W; ++r)
        {
            aco::EvaluatedSolution es;
            es.cost = elite[r].cost;
            es.sol  = std::move(elite[r].sol);
            ranked.push_back(std::move(es));
        }

        // 4.1) FINAL CONSISTENCY: za top-W rješenja osiguraj DP-optimalan car[] i cost
        // (Ovo je jeftino jer je W mali, a sprječava da feromoni uče iz nekonzistentnih car-assignments)
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

            // Pošto se cost može promijeniti nakon reassignCars, ponovno sortiraj ranked
            std::sort(ranked.begin(), ranked.end(),
                      [](const aco::EvaluatedSolution& a, const aco::EvaluatedSolution& b) 
                      {
                          return a.cost < b.cost;
                      });
        }

        // Sada je iterBest stvarno ranked[0] nakon DP konsolidacije
        const double iterBestCost = ranked[0].cost;
        rr.iterBestCost.push_back(iterBestCost);

        // zapamti best prije ove iteracije (za finalni globalImproved)
        const double prevBestCost = rr.best.cost;

        // iteracijsko poboljšanje (EPS)
        const bool improvedByIter = (iterBestCost < rr.best.cost - improveEps);

        // update rr.best ODMAH (prije stagnation), ali KOPIJOM da ranked ostane netaknut za feromone
        if (improvedByIter)
            rr.best = ranked[0];   // <-- copy (nemoj std::move)
        
        // Ažuriraj lokalni brojač stagnacije za kick/LK/3-opt
        if (improvedByIter) {
            itersSinceGlobalImprovement = 0;
        } else {
            ++itersSinceGlobalImprovement;
        }

        // iterBestPost = najbolji cost u ovoj iteraciji (nakon svih akcija)
        // Počinje s iterBestCost, ali se može poboljšati kroz intensification, polish, kick
        // VAŽNO: iterBestPost je najbolji cost U TOJ ITERACIJI, ne kumulativni minimum
        // 
        // FINALNA CIJENA ITERACIJE:
        // iterBestPost je prava/finalna cijena iteracije - najbolji cost nakon svih akcija:
        // - Konstrukcija mrava
        // - Local Search na top-LSW elite mrave
        // - Intensification (dodatni LS na best)
        // - Polish (LK/3-opt na stagnaciju)
        // - Kick (double-bridge + LS na stagnaciju)
        // 
        // KORIŠTENJE:
        // - Ispisuje se u CLI/log kao "it best" (iteration best)
        // - Snima se u recorder kao iterBestPostCost
        // - NE koristi se za ažuriranje feromona (feromoni koriste ranked top-W rješenja)
        double iterBestPost = iterBestCost;
        
        // Zapamti cost intensification poboljšanja (ako postoji) za kasnije
        double intensificationCost = std::numeric_limits<double>::infinity();
        
        // 4.5) INTENSIFIKACIJA: Dodatni LS ciklus samo za global best rješenje
        // 
        // INTENSIFIKATOR:
        // - Komponenta koja primjenjuje dodatni Local Search ciklus na global best rješenje
        // - Aktivira se kada je best upravo promijenjen (improvedByIter) ili u prvih 3 iteracije nakon stagnacije
        // - Omogućava dublju eksploataciju bez prevelikog utroška vremena (radi samo na best rješenju)
        // 
        // KAKO SE PRIKLJUČUJE:
        // - Kreira se u main.cpp factory funkciji i prosljeđuje u Colony konstruktor
        // - Colony automatski kreira intenzifikator ako nije proslijeđen ali postoje potrebne komponente
        // - Koristi se ovdje nakon što se ažurira rr.best s iterBestCost
        // 
        // LOGIKA:
        // - intensify() provjerava uvjete aktivacije i primjenjuje LS + DP na best
        // - Ako je poboljšanje, intensify() direktno ažurira rr.best (pass-by-reference)
        // - Colony samo provjerava je li se cost promijenio i ažurira tracking varijable
        bool intensifierActiveThisIteration = false;
        bool disablePheromoneLearningOnActivationThisIteration = false;
        if (intensifier_ && std::isfinite(rr.best.cost))
        {
            intensifierActiveThisIteration =
                improvedByIter ||
                (itersSinceGlobalImprovement > 0 &&
                 itersSinceGlobalImprovement <= params_.intensifierMaxStagnationIterations);
            // Zapamti cost prije intenzifikacije (za provjeru je li se dogodilo poboljšanje)
            const double costBefore = rr.best.cost;
            
            // Pozovi intenzifikator - ažurira rr.best direktno ako je poboljšanje
            const double newCost = intensifier_->intensify(rr.best, improvedByIter, 
                                                          itersSinceGlobalImprovement, improveEps);
            
            // intensify() već ažurira rr.best ako je poboljšanje, samo provjerimo je li se dogodilo
            // (provjeravamo rr.best.cost jer je pass-by-reference i može biti ažuriran)
            if (std::isfinite(newCost) && rr.best.cost < costBefore - improveEps)
            {
                intensificationCost = rr.best.cost;  // Zapamti cost za tracking (iterBestPost)
                
                // Ako je ovo poboljšanje, resetiraj stagnaciju (izlaz iz stagnacije)
                if (itersSinceGlobalImprovement > 0) {
                    itersSinceGlobalImprovement = 0;
                }
            }

            // Istraživačko logiranje: detektiraj početak nove aktivacije (brojač se povećao)
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

            // Istraživačko logiranje: na deaktivaciju zapiši (run, activation_index, cost_start, cost_end, improved)
            if (intensifier_->wasJustDeactivated() && researchLogger_)
            {
                const int activationCount = intensifier_->getActivationCount();
                const bool improvedThisActivation = (std::isfinite(rr.best.cost) && std::isfinite(lastActivationStartCost) &&
                                                    rr.best.cost < lastActivationStartCost - improveEps);
                researchLogger_->onDeactivation(runIndex_, activationCount, lastActivationStartCost, rr.best.cost, improvedThisActivation);
            }
            
            // ISPIS KADA SE INTENZIFIKATOR DEAKTIVIRA:
            // Ispiši svaki put kada se intenzifikator isključi (bio aktivan, sada nije)
            if (intensifier_->wasJustDeactivated())
            {
                const int activationCount = intensifier_->getActivationCount();
                if (recorder_) {
                    recorder_->recordIntensifierDeactivationEvent(runIndex_, it, activationCount, rr.best.cost);
                }
                vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                       ", intensifier deactivated (activation count: " + std::to_string(activationCount) + ")");
            }
            
            // POBOLJŠANO RESETIRANJE FEROMONA NAKON PRVE DEAKTIVACIJE:
            // Kombinacija preporučenih pristupa:
            // 1. Selektivno resetiranje (samo feromoni na best tour-u) + smoothing
            // 2. Resetiranje samo ako best NIJE poboljšan nakon intenzifikacije
            // 3. Arhiviranje global besta i resetiranje global besta na infinity
            if (intensifier_->shouldResetPheromones() && pherModel_)
            {
                // Provjeri je li best poboljšan nakon intenzifikacije
                // Ako nije poboljšan, resetiraj feromone (selektivno + smoothing)
                const bool bestNotImproved = (rr.best.cost >= costBefore - improveEps);
                
                if (bestNotImproved)
                {
                    // 1. ARHIVIRANJE GLOBAL BESTA (ako je opcija omogućena):
                    // Arhiviraj trenutni global best prije resetiranja
                    // Koristimo rr.best.sol i costBefore jer ćemo resetirati rr.best.cost na infinity
                    if (params_.intensifierArchiveAndResetBest && std::isfinite(costBefore) && 
                        (int)rr.best.sol.node.size() == inst_->n())
                    {
                        EvaluatedSolution bestToArchive = rr.best;
                        bestToArchive.cost = costBefore;  // Koristimo costBefore jer je to stvarni cost prije resetiranja
                        pherModel_->addToArchive(bestToArchive);
                        vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                               ", archived global best (cost: " + std::to_string(costBefore) + ")");
                        
                        // Resetiraj global best na infinity (omogućava potpuno novo istraživanje)
                        rr.best.cost = std::numeric_limits<double>::infinity();
                        bestWasResetToInfinity = true;  // Označi da je best resetiran
                        vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                               ", reset global best to infinity (allowing fresh exploration)");
                    }
                    
                    // 2. SELEKTIVNO RESETIRANJE FEROMONA (smoothing na best tour-u):
                    // Umjesto potpunog resetiranja, smanjujemo feromone samo na best tour-u
                    // Koristimo smoothing faktor gamma (0.0 = potpuno resetiranje, 1.0 = bez promjene)
                    if (std::isfinite(costBefore))
                    {
                        // Koristimo rr.best.sol ako je još valjan (prije resetiranja cost-a), 
                        // inače koristimo ranked[0] ako postoji
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
                            // Fallback: potpuno resetiranje ako nemamo valjan tour
                            pherModel_->reset();
                            vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                                   ", full pheromone reset (no valid tour for selective reset)");
                        }
                    }
                    else
                    {
                        // Fallback: potpuno resetiranje ako nemamo costBefore
                        pherModel_->reset();
                        vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                               ", full pheromone reset (no costBefore available)");
                    }
                }
                else
                {
                    // Best JE poboljšan, ne resetiraj feromone
                    vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                           ", intensifier deactivated but best improved (no pheromone reset)");
                }
            }
            
            // Ažuriraj interno stanje intenzifikatora za sljedeći poziv
            // (nakon što se provjeri deaktivacija i resetiranje feromona)
            intensifier_->updateState();
        }
        // Kad je intenzifikator isključen (--intensifier 0), nema dodatnog LS-a isključivo na
        // global best ovdje; ostaje LS na elitnim mravima u glavnom pipelineu iznad.

        // Ažuriraj iterBestPost ako je intensification poboljšao u ovoj iteraciji
        if (std::isfinite(intensificationCost) && intensificationCost < iterBestPost)
        {
            iterBestPost = intensificationCost;
        }
        bool skipKick = false;

        // STAGNATION - LK, 3-opt or KICK
        // (kick/strong LS kasnije ovdje, temeljeno na lokalnom itersSinceGlobalImprovement)
        // --- Kick (double-bridge) na stagnaciju ---
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

            // 1) polish na top elite mrave (ne samo rr.best)
            // Koristimo isti pristup kao za LS: primenjujemo polish na nekoliko top elite mrava
            // Da bi polish mogao poboljšavati više mrava, ne samo jednog elitnog
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

                // Adaptivni broj elite mrava za polish: top-3 do top-5 ovisno o veličini instance
                // Polish je skup (6 pokušaja × DP verify po mravu), pa fokusiramo se na top mrave
                // gdje je najveća vjerojatnost poboljšanja
                const int N = inst_->n();
                int targetPolishCount;
                if (N <= 100) {
                    targetPolishCount = 3;  // Manje instance: top-3
                } else if (N <= 200) {
                    targetPolishCount = 4;  // Srednje instance: top-4
                } else {
                    targetPolishCount = 5;  // Veće instance: top-5
                }
                const int polishLimit = std::min({targetPolishCount, K, (int)ranked.size()});

                for (int i = 0; i < polishLimit; ++i)
                {
                    // Napravi kopiju da ne menjamo originalni ranked (koristi se za feromone)
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

                            // Ažuriraj rr.best ako je globalno poboljšanje (prvo proveri globalno)
                            bool globalImproved = false;
                            if (cand.cost < rr.best.cost - improveEps)
                            {
                                rr.best = cand;  // kopija, ne move (treba i za ranked)
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

                            // Ažuriraj ranked[i] ako je poboljšanje (lokalno ili globalno)
                            if (cand.cost < ranked[i].cost - improveEps)
                            {
                                ranked[i] = std::move(cand);
                            }
                        }
                    }
                }

                // Sortiraj ranked nakon svih polish poboljšanja (da feromoni uče iz sortiranih rešenja)
                if (polishLimit > 0)
                {
                    std::sort(ranked.begin(), ranked.end(),
                        [](const aco::EvaluatedSolution& a, const aco::EvaluatedSolution& b) 
                        {
                            return a.cost < b.cost;
                        });
                    // Ažuriraj iterBestPost nakon sortiranja
                    if (ranked[0].cost < iterBestPost) iterBestPost = ranked[0].cost;
                }
            }

            // 2) kick only if polish failed
            // Primjenjujemo kick na top-2 do top-3 mrava za veću diverzifikaciju
            // Kick je jeftiniji od polisha (perturbacija + LS), pa možemo primijeniti na više mrava
            if (!improvedByStag)
            {
                // Adaptivni broj mrava za kick: top-2 do top-3 ovisno o veličini instance
                const int N = inst_->n();
                int targetKickCount;
                if (N <= 100) {
                    targetKickCount = 2;  // Manje instance: top-2
                } else {
                    targetKickCount = 3;  // Veće instance: top-3
                }
                const int kickLimit = std::min({targetKickCount, (int)ranked.size()});

                for (int i = 0; i < kickLimit; ++i)
                {
                    if (recorder_ && i == 0) recorder_->incKickCall(runIndex_);
                    if (i == 0) {
                        vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + ", kick call (on " + std::to_string(kickLimit) + " solutions).");
                    }

                    aco::EvaluatedSolution cand = ranked[i];
                    // Generiraj različit seed za svaki mrav kombinacijom i s drugim argumentima
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

                    // Ažuriraj ranked[i] ako je poboljšanje
                    if (std::isfinite(cand.cost) && cand.cost < ranked[i].cost - improveEps)
                    {
                        ranked[i] = std::move(cand);
                    }
                }

                // Sortiraj ranked nakon kick-a (da feromoni uče iz sortiranih rešenja)
                if (kickLimit > 0)
                {
                    std::sort(ranked.begin(), ranked.end(),
                        [](const aco::EvaluatedSolution& a, const aco::EvaluatedSolution& b) 
                        {
                            return a.cost < b.cost;
                        });
                    // Ažuriraj iterBestPost nakon sortiranja
                    if (ranked[0].cost < iterBestPost) iterBestPost = ranked[0].cost;
                }
            }
            
            // Resetuj lokalni brojač stagnacije nakon aktivacije kick/LK/3-opt (cooldown)
            // noImprove_ u pheromone_model ostaje netaknut (koristi se za smoothing)
            itersSinceGlobalImprovement = 0;
        }


        // finalni update feromona
        // FINAL globalImproved: je li rr.best bolji nego na početku iteracije
        // (Ovo je već provjereno i ispisano gore kada se iterBestPost ažurira)
        const bool globalImprovedFinal = (rr.best.cost < prevBestCost - improveEps);

        // Waiting-time Weibull: bilježi pravo poboljšanje (neosjetljivo na archive+reset)
        const bool trueImprovementThisIter =
            std::isfinite(iterBestPost) && iterBestPost < allTimeBestCost - improveEps;
        if (trueImprovementThisIter) {
            allTimeBestCost = iterBestPost;
            if (recorder_) recorder_->recordGlobalBestImprovement(runIndex_, it);
        }
        
        // PROVJERA ARHIVIRANOG BESTA NAKON RESETIRANJA:
        // Ako je best bio resetiran na infinity i sada je pronađen novi best,
        // provjeri je li novi best bolji od arhiviranog i zamijeni ako je bolji.
        // Ako je novi best bolji, ponovno resetiraj best na infinity i feromone za kontinuirano istraživanje.
        if (bestWasResetToInfinity && globalImprovedFinal && std::isfinite(rr.best.cost) && pherModel_)
        {
            const double bestArchivedCost = pherModel_->getBestArchivedCost();
            
            // Ako je novi best bolji od arhiviranog, zamijeni u archive-u i ponovno resetiraj
            if (std::isfinite(bestArchivedCost) && rr.best.cost < bestArchivedCost - improveEps)
            {
                // Novi best je bolji - dodaj u archive (automatski će zamijeniti stari jer je archive sortiran)
                pherModel_->addToArchive(rr.best);
                vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                       ", new best (" + std::to_string(rr.best.cost) + ") better than archived (" + 
                       std::to_string(bestArchivedCost) + ") - updated archive");
                
                // PONOVNO RESETIRANJE: Resetiraj best na infinity i feromone za kontinuirano istraživanje
                if (params_.intensifierArchiveAndResetBest)
                {
                    // Zapamti tour prije resetiranja za selektivno resetiranje feromona
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
                    
                    // Resetiraj best na infinity (omogućava potpuno novo istraživanje)
                    rr.best.cost = std::numeric_limits<double>::infinity();
                    bestWasResetToInfinity = true;  // Zadrži flag jer ćemo ponovno resetirati
                    vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                           ", reset global best to infinity again (continuous exploration)");
                    
                    // Resetiraj feromone (selektivno ili potpuno)
                    if (tourToReset)
                    {
                        const double gamma = params_.intensifierPheromoneReductionGamma;
                        pherModel_->reducePheromonesOnTour(*tourToReset, gamma);
                        vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                               ", selectively reduced pheromones on best tour again (gamma: " + std::to_string(gamma) + ")");
                    }
                    else
                    {
                        // Fallback: potpuno resetiranje ako nemamo valjan tour
                        pherModel_->reset();
                        vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                               ", full pheromone reset again (no valid tour for selective reset)");
                    }
                }
            }
            else if (std::isfinite(bestArchivedCost))
            {
                // Novi best nije bolji - arhivirani ostaje netaknut, ne resetiraj
                vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                       ", new best (" + std::to_string(rr.best.cost) + ") not better than archived (" + 
                       std::to_string(bestArchivedCost) + ") - keeping archived");
                
                // Resetiraj flag jer nismo resetirali (best ostaje netaknut)
                bestWasResetToInfinity = false;
            }
            else
            {
                // Nema arhiviranog besta (prvi put), resetiraj flag
                bestWasResetToInfinity = false;
            }
        }
        
        // AŽURIRANJE FEROMONA:
        // Feromoni se ažuriraju s ranked (top-W rješenja) i rr.best, NE s iterBestPost
        // ranked sadrži top-W rješenja nakon LS-a i DP konsolidacije (prije stagnation akcija)
        // iterBestPost je finalna cijena iteracije (nakon svih akcija), ali se ne koristi za feromone
        // 
        // Razlog: feromoni uče iz top-W rješenja koja su već DP-optimalna i sortirana,
        // što omogućava bolje učenje strukture rješenja
        bool smoothingActive = false;
        if (pherModel_) {
            // Provjeri noImprove_ prije poziva (da vidimo je li se resetovao)
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
                // NE koristi iterBestPost - feromoni uče iz top-W rješenja, ne iz finalne cijene iteracije
                pherModel_->onIterationEndRanked(it, ranked, rr.best, globalImprovedFinal);
            }
            
            // Provjeri noImprove_ poslije poziva (da vidimo je li se resetovao)
            const int noImproveAfter = pherModel_->getStagnationCounter();
            
            // Ispis ako se noImprove_ resetovao (smanjio)
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
            
            // Check if smoothing was applied in this iteration (iz pheromone_model)
            bool smoothingFromModel = pherModel_->wasSmoothingApplied();
            
            // Smoothing također raste kako se približavamo stagnaciji (isto kao LSW)
            // Smoothing se aktivira na 2× stagnation prag u pheromone_model,
            // koristimo noImprove_ iz pheromone_model (jedinstveni brojač stagnacije)
            // Smoothing počinje rasti kada se približavamo 2× stagnation pragu
            const int smoothingStagPrag = (params_.stagnation > 0) ? (2 * params_.stagnation) : 0;
            const int noImprovePher = pherModel_->getStagnationCounter();
            if (smoothingStagPrag > 0 && noImprovePher >= smoothingStagPrag - 10) {
                // Smoothing je aktivan kada se približavamo 2× stagnation pragu (10 iteracija prije)
                smoothingActive = true;
            } else if (smoothingFromModel) {
                // Smoothing je aktivan iz pheromone_model (na 2× stagnation pragu)
                smoothingActive = true;
            }
        }

        // elapsed isto na kraju petlje
        // elapsed nakon stagnation akcija (da uključi LK/kick/LS):
        const auto now = clock::now();
        const double elapsedMs =
            std::chrono::duration<double, std::milli>(now - t0).count();

        // Ažuriraj LSW na kraju iteracije (za sljedeću iteraciju i ispis)
        // Hibridni pristup: LSW raste postepeno kroz sve iteracije + dodatni rast u stagnaciji
        // Napomena: Ako je LSW_FIXED > 0, LSW je već postavljen na LSW_FIXED u liniji 241, tako da se ovdje ne mijenja
        if (LSW_FIXED == 0) {
            const int LSW_START = 2;
            const int LSW_END = std::min(K, 7); // Povećano s 5 na 7 za više LS-a u stagnaciji
            
            // 1) Base LSW: postepeni rast kroz sve iteracije (80% od maksimalnog rasta za bržu konvergenciju)
            const double progress = (double)it / (double)params_.iterations;
            const double baseLSW = LSW_START + (LSW_END - LSW_START) * 0.8 * progress;
            
            // 2) Stagnation bonus: dodatni rast u stagnaciji (20% od maksimalnog rasta, jer base LSW je već 80%)
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
            
            // 3) Kombiniraj base LSW i stagnation bonus
            double targetLSW = baseLSW + stagnationBonus;
            currentLSW = (int)std::round(targetLSW);
            currentLSW = std::clamp(currentLSW, LSW_START, LSW_END);
            
            // Ažuriraj LSW za ispis i recording
            LSW = currentLSW;
        }
        
        // FINALNA CIJENA ITERACIJE - provjera na kraju iteracije
        // iterBestPost = najbolji cost u ovoj iteraciji (nakon svih akcija: konstrukcija + LS + intensification + polish + kick)
        // iterBestPost se već ažurira tijekom iteracije (intensification, polish, kick)
        // Na kraju iteracije, provjerimo je li ranked[0].cost bolji od trenutnog iterBestPost
        // VAŽNO: iterBestPost je najbolji cost U TOJ ITERACIJI, ne kumulativni minimum
        // Ne koristimo rr.best.cost jer je to kumulativni minimum kroz sve iteracije
        // 
        // Ovo je prava/finalna cijena iteracije koja se koristi za:
        // - Ispis u CLI/log (kao "it best")
        // - Snimanje u recorder (kao iterBestPostCost)
        // - NE koristi se za ažuriranje feromona (feromoni koriste ranked top-W rješenja)
        if (!ranked.empty() && std::isfinite(ranked[0].cost))
        {
            iterBestPost = std::min(iterBestPost, ranked[0].cost);
        }
        
        // Ako je iterBestPost bolji od rr.best.cost, ažuriraj rr.best (cijelo rješenje, ne samo cost)
        // rr.best.cost mora biti kumulativni minimum (najbolji cost kroz sve iteracije)
        // VAŽNO: rr.best se može ažurirati na liniji 497-498 (ako je iterBestCost < rr.best.cost)
        // ili tijekom iteracije (intensification/polish/kick), ali ovdje provjeravamo je li iterBestPost
        // (najbolji cost u ovoj iteraciji nakon svih akcija) bolji od trenutnog rr.best.cost
        EvaluatedSolution bestPostSolution;  // Rješenje koje odgovara iterBestPost
        bool hasBestPostSolution = false;
        
        if (std::isfinite(iterBestPost) && iterBestPost < rr.best.cost - improveEps)
        {
            const double oldBestCost = rr.best.cost;
            // Ako je iterBestPost iz ranked[0], koristi ranked[0]
            // Inače, rr.best je već ažuriran tijekom iteracije (intensification/polish/kick)
            // ili na liniji 497-498 (ako je iterBestCost < rr.best.cost)
            if (!ranked.empty() && std::isfinite(ranked[0].cost) && 
                std::abs(iterBestPost - ranked[0].cost) < improveEps)
            {
                // iterBestPost dolazi iz ranked[0]
                rr.best = ranked[0];  // kopija cijelog rješenja (sol + cost), ne move
                bestPostSolution = ranked[0];
                hasBestPostSolution = true;
            }
            else
            {
                // iterBestPost == rr.best.cost, rr.best je već ažuriran (intensification/polish/kick ili na liniji 497-498)
                bestPostSolution = rr.best;
                hasBestPostSolution = true;
            }
            
            // ISPIS KADA SE POSTIGNE NOVA GLOBAL BEST:
            // Ispiši svaki put kada se postigne nova global best (poboljšanje kumulativnog minimuma)
            vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                   ", NEW GLOBAL BEST: " + std::to_string(oldBestCost) + 
                   " -> " + std::to_string(rr.best.cost));
        }
        else if (std::isfinite(iterBestPost))
        {
            // iterBestPost nije bolji od rr.best, ali zapamti rješenje za provjeru arhiviranog besta
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
        
        // PROVJERA ARHIVIRANOG BESTA NAKON ITERBESTPOST:
        // Ako je best bio resetiran na infinity i iterBestPost je bolji od arhiviranog,
        // provjeri je li iterBestPost (best post) bolji od arhiviranog besta.
        // Ovo je važno jer iterBestPost može biti bolji od rr.best koji je arhiviran ranije.
        if (bestWasResetToInfinity && hasBestPostSolution && std::isfinite(bestPostSolution.cost) && pherModel_)
        {
            const double bestArchivedCost = pherModel_->getBestArchivedCost();
            
            // Ako je iterBestPost (best post) bolji od arhiviranog, zamijeni u archive-u
            if (std::isfinite(bestArchivedCost) && bestPostSolution.cost < bestArchivedCost - improveEps)
            {
                // iterBestPost je bolji - dodaj u archive (automatski će zamijeniti stari jer je archive sortiran)
                pherModel_->addToArchive(bestPostSolution);
                vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                       ", iterBestPost (" + std::to_string(bestPostSolution.cost) + ") better than archived (" + 
                       std::to_string(bestArchivedCost) + ") - updated archive with best post");
                
                // PONOVNO RESETIRANJE: Resetiraj best na infinity i feromone za kontinuirano istraživanje
                if (params_.intensifierArchiveAndResetBest)
                {
                    // Zapamti tour prije resetiranja za selektivno resetiranje feromona
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
                    
                    // Resetiraj best na infinity (omogućava potpuno novo istraživanje)
                    rr.best.cost = std::numeric_limits<double>::infinity();
                    bestWasResetToInfinity = true;  // Zadrži flag jer ćemo ponovno resetirati
                    vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                           ", reset global best to infinity again (best post better than archived)");
                    
                    // Resetiraj feromone (selektivno ili potpuno)
                    if (tourToReset)
                    {
                        const double gamma = params_.intensifierPheromoneReductionGamma;
                        pherModel_->reducePheromonesOnTour(*tourToReset, gamma);
                        vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                               ", selectively reduced pheromones on best post tour (gamma: " + std::to_string(gamma) + ")");
                    }
                    else
                    {
                        // Fallback: potpuno resetiranje ako nemamo valjan tour
                        pherModel_->reset();
                        vprint("run: " + std::to_string(runIndex_) + ", iteration: " + std::to_string(it) + 
                               ", full pheromone reset again (no valid tour for selective reset)");
                    }
                }
            }
        }
        
        if (recorder_) {
            // Weibull k iz costova mrava (ranked) za ovu iteraciju
            const double weibullK = weibullKFromCosts(ranked);
            // SNIMANJE U RECORDER:
            // - iterBestCost = najbolji cost prije stagnation akcija (nakon konstrukcije + LS)
            // - iterBestPost = najbolji cost u ovoj iteraciji (nakon svih akcija) - FINALNA CIJENA ITERACIJE
            // - rr.best.cost = kumulativni minimum (najbolji cost kroz sve iteracije)
            //   rr.best.cost je već kumulativni minimum jer se ažurira samo kada se nađe bolje rješenje
            // 
            // iterBestPost je prava/finalna cijena iteracije koja se snima za analizu
            recorder_->recordIteration2(runIndex_, it, iterBestCost, iterBestPost, rr.best.cost, elapsedMs, weibullK);
            // Snimi q0, LSW i smoothing za ovu iteraciju
            recorder_->recordIterationParams(runIndex_, it, q0Value, LSW, smoothingActive);

            // Diversity stats za H5/H6 hipoteze (iz istih costova kao weibullKFromCosts)
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
        // stagCount = lokalni itersSinceGlobalImprovement (za kick/LK, resetuje se nakon aktivacije)
        // noImprovCount = noImprove_ iz pheromone_model (za smoothing, resetuje se na global improvement ili smoothing)
        int noImprovCount = 0;
        if (pherModel_) {
            noImprovCount = pherModel_->getStagnationCounter();
        }
        
        // Izračunaj medijan i srednju vrijednost global best-a za ovu iteraciju iz svih runova
        // Samo runovi koji su već izračunali svoju vrijednost (sporiji se ignoriraju)
        double medianGlobalBest = StatsRow::nan();
        double meanGlobalBest = StatsRow::nan();
        if (recorder_)
        {
            medianGlobalBest = recorder_->getMedianGlobalBestForIteration(it);
            meanGlobalBest = recorder_->getMeanGlobalBestForIteration(it);
        }
        
        // ISPIS U CLI/LOG:
        // - "it best" = iterBestPost (finalna cijena iteracije - najbolji cost nakon svih akcija)
        // - "best" = rr.best.cost (kumulativni minimum - najbolji cost kroz sve iteracije)
        // - "median" = medijan global best-a za ovu iteraciju iz svih runova (samo izračunati)
        // - "mean" = srednja vrijednost global best-a za ovu iteraciju iz svih runova (samo izračunati)
        std::ostringstream oss;
        oss << "run: " << runIndex_
            << ", iteration: " << it
            << ", stagCount: " << itersSinceGlobalImprovement
            << ", q0: " << std::fixed << std::setprecision(3) << q0Value
            << ", LSW: " << LSW
            << ", noImprovCount: " << noImprovCount
            << ", smoothing: " << (smoothingActive ? 1 : 0)
            << ", it best: " << std::setprecision(3) << iterBestPost  // FINALNA CIJENA ITERACIJE
            << ", best: " << std::setprecision(3) << rr.best.cost;     // KUMULATIVNI MINIMUM
        
        if (std::isfinite(medianGlobalBest))
        {
            oss << ", median: " << std::setprecision(3) << medianGlobalBest;
        }
        
        if (std::isfinite(meanGlobalBest))
        {
            oss << ", mean: " << std::setprecision(3) << meanGlobalBest;
        }
        
        vprint(oss.str());

        // Cache statistics: ispis svakih 50 iteracija ili na zadnjoj iteraciji
        // Ispisujemo samo za verbose run (trenutno odabrani run)
        if (verbose_ && dp_ && (it % 50 == 0 || it == params_.iterations - 1))
        {
            // Provjeri je li dp_ CachedFixedTourCarAssigner
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
                
                // Ispis samo za verbose run
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

    // za svaki slucaj finalna konsolidacija (global best vs DP)
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
    
    // FINALNA PROVJERA ARHIVIRANOG BESTA (global best vs archived/cache):
    // Na kraju run-a provjeri je li arhivirani best bolji od trenutnog rr.best.
    // Ako je bolji, zamijeni rr.best s arhiviranim bestom kao finalno rješenje.
    // Ovo je važno jer se best može resetirati na infinity tijekom run-a, ali arhivirani best
    // ostaje sačuvan i može biti bolji od finalnog rr.best.
    if (params_.intensifierArchiveAndResetBest && pherModel_)
    {
        EvaluatedSolution archivedBest;
        if (pherModel_->getBestArchivedSolution(archivedBest))
        {
            // Provjeri je li arhivirani best bolji od trenutnog rr.best
            // Koristimo isti epsilon kao u ostatku koda (1e-12)
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
                    // Arhivirani best je bolji - zamijeni rr.best s arhiviranim bestom
                    rr.best = archivedBest;
                    
                    // Finalna konsolidacija za arhivirani best (osiguraj DP-optimalan car[])
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

    // Final cache statistics na kraju runa
    // Snimimo u recorder za sve runove (ne samo verbose)
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
            
            // Snimi u recorder za CSV ispis
            if (recorder_)
            {
                recorder_->setCacheStats(runIndex_, evalHitRate, reasHitRate);
            }
            
            // Ispis samo za verbose run
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
            if (u == v) continue; // skip dijagonala
            
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

    // Prosječni return cost (za surogat s return članom — poboljšanje korelacije na većim instancama)
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
