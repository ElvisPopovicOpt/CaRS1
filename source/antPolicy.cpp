#include <cmath>
#include <limits>
#include <stdexcept>
#include <parser.hpp>
#include <rng.hpp>
#include <interfaces.hpp>
#include <params_cli.hpp>
#include <candidate_list.hpp>
#include <antPolicy.hpp>

namespace aco 
{

AntPolicyCandidateListRoulette::AntPolicyCandidateListRoulette(
    std::shared_ptr<const CandidateListCache> cl,
    double epsilonEta,
    double returnFactor,
    bool sampleCars,
    bool enforceNoRerent,
    double q0)
    : cl_(std::move(cl))
    , eps_(epsilonEta)
    , returnFactor_(returnFactor)
    , sampleCars_(sampleCars)
    , enforceNoRerent_(enforceNoRerent)
    , q0_(std::max(0.0, std::min(1.0, q0))) // NEW
{}


static inline double eta(double eps, double cost) 
{
    return 1.0 / (eps + cost);
}

Solution AntPolicyCandidateListRoulette::construct(const AntContext& ctx) const
{
    if (!ctx.inst || !ctx.params || !ctx.pher || !ctx.rng)
        throw std::runtime_error("AntPolicyCandidateListRoulette: ctx.inst/params/pher/rng must be set.");

    auto& rng = *ctx.rng;
    const auto& inst = *ctx.inst;
    const int N = inst.n();
    const int C = inst.cars();

    std::vector<uint8_t> carLocked;
    if (enforceNoRerent_) 
    {
        carLocked.assign((size_t)C, 0);
    }

    Solution s;
    s.node.resize((size_t)N);
    s.car.resize((size_t)N);

    std::vector<uint8_t> visited((size_t)N, 0);
    s.node[0] = 0;
    visited[0] = 1;

    // Kretanje po turi (odabir sljedećeg čvora): alpha/beta za čvorove
    const double alphaNodes = ctx.params->alphaNodes;
    const double betaNodes  = ctx.params->betaNodes;
    // Izbor automobila na bridu: alpha/beta za automobile
    const double alphaCars = ctx.params->alphaCars;
    const double betaCars   = ctx.params->betaCars;

    // Return komponenta je aktivna samo ako instanca ima return troškove
    // i pheromone accessor podržava return feromone.
    const bool useReturn =
        inst.hasReturnCosts() && ctx.pher->hasReturnPheromones();
    // Return komponenta u POLICYJU (konstrukciji) se koristi samo ako:
    //  - instanca ima return troškove,
    //  - pheromone accessor ima tauReturn,
    //  - i returnFactor_ > 0.0.
    // Kad je returnFactor_ == 0.0, kompletan return dio se preskače (nema overhead-a).
    const bool useReturnPolicy = useReturn && (returnFactor_ > 0.0);

    // Stanje potrebno za switch trošak (ekvivalentno CostModelCars::evaluate):
    int prevCar    = -1; // nepoznato na pos=0 => nema switch penala
    int rentalNode = 0;  // lastCarNode u cost modelu; start je 0

    auto choose_from_pairs = [&](int currentNode,
                                 const std::vector<int>& nextCandidates,
                                 int forcedNextOrNeg1,
                                 int prevCarState,
                                 int rentalNodeState,
                                 bool isLastStepToStart) -> std::pair<int,int>
    {
        struct Item { int j; int c; double w; };
        std::vector<Item> items;
        items.reserve((size_t)nextCandidates.size() * (size_t)std::max(1, C));

        double sum = 0.0;
        // Precompute za return komponentu u policyju (koristi se samo ako useReturnPolicy==true).
        double alphaR = 0.0;
        double betaR  = 0.0;
        if (useReturnPolicy) 
        {
            alphaR = alphaCars * returnFactor_;
            betaR  = betaCars * returnFactor_;
        }

        if (!sampleCars_) 
        {
            struct NItem { int j; double w; };
            std::vector<NItem> nitems;
            nitems.reserve(nextCandidates.size());

            double nsum = 0.0;
            for (int j : nextCandidates) 
            {
                if (forcedNextOrNeg1 >= 0 && j != forcedNextOrNeg1) continue;

                // Agregacija za čvor (brid): tauNode, etaNode za kretanje po turi (alphaNodes, betaNodes)
                double tauNode = 0.0;
                double minCost = std::numeric_limits<double>::max();
                for (int c = 0; c < C; ++c) 
                {
                    if (enforceNoRerent_ && c != prevCarState && carLocked[(size_t)c]) continue;
                    tauNode += ctx.pher->tauMove(c, currentNode, j);
                    minCost = std::min(minCost, inst.travelCost(c, currentNode, j));
                }
                const double etaNode = eta(eps_, minCost);
                const double W_node = std::pow(tauNode, alphaNodes) * std::pow(etaNode, betaNodes);
                if (!(W_node > 0.0) || !std::isfinite(W_node)) continue;

                // Suma po autima (izbor auta: alphaCars, betaCars); DP će kasnije dodijeliti auto
                double wj = 0.0;
                for (int c = 0; c < C; ++c) 
                {
                    if (enforceNoRerent_ && c != prevCarState && carLocked[(size_t)c]) continue;

                    const double tMove = ctx.pher->tauMove(c, currentNode, j);
                    const double hMove = eta(eps_, inst.travelCost(c, currentNode, j));
                    double wc = std::pow(tMove, alphaCars) * std::pow(hMove, betaCars);
                    if (!(wc > 0.0) || !std::isfinite(wc)) continue;

                    if (useReturn && returnFactor_ > 0.0) 
                    {
                        if (prevCarState >= 0 && c != prevCarState) 
                        {
                            const double rc = inst.returnCost(prevCarState, currentNode, rentalNodeState);
                            const double tR = ctx.pher->tauReturn(prevCarState, currentNode, rentalNodeState);
                            const double hR = eta(eps_, rc);
                            wc *= std::pow(tR, alphaR) * std::pow(hR, betaR);
                        }

                        if (isLastStepToStart) 
                        {
                            const int rentalAfter = (prevCarState >= 0 && c != prevCarState)
                                                    ? currentNode
                                                    : rentalNodeState;
                            const double rcFinal = inst.returnCost(c, 0, rentalAfter);
                            const double tFinal  = ctx.pher->tauReturn(c, 0, rentalAfter);
                            const double hFinal  = eta(eps_, rcFinal);
                            wc *= std::pow(tFinal, alphaR) * std::pow(hFinal, betaR);
                        }
                    }

                    if (!(wc > 0.0) || !std::isfinite(wc)) continue;
                    wj += wc;
                }

                const double w = W_node * wj;
                if (w > 0.0 && std::isfinite(w)) 
                {
                    nitems.push_back({j, w});
                    nsum += w;
                }
            }

            if (nitems.empty() || !(nsum > 0.0) || !std::isfinite(nsum)) 
            {
                int j = (forcedNextOrNeg1 >= 0)
                    ? forcedNextOrNeg1
                    : nextCandidates[rng.uniformInt((int)nextCandidates.size())];
                return {j, 0}; // dummy car
            }

            // exploitation jump (argmax) ---
            // OPTIMIZACIJA: Koristi dinamički q0 iz ctx.adaptiveQ0 ako je dostupan
            const double q0 = (ctx.adaptiveQ0 >= 0.0) ? ctx.adaptiveQ0 : q0_;
            if (q0 > 0.0 && rng.uniform01() < q0) 
            {
                size_t bestIdx = 0;
                double bestW = nitems[0].w;
                for (size_t t = 1; t < nitems.size(); ++t) {
                    if (nitems[t].w > bestW) { bestW = nitems[t].w; bestIdx = t; }
                }
                return {nitems[bestIdx].j, 0};
            }

            double r = rng.uniform01() * nsum;
            double acc = 0.0;
            for (const auto& it : nitems) 
            {
                acc += it.w;
                if (acc >= r) return {it.j, 0};
            }
            return {nitems.back().j, 0};
        }


        for (int j : nextCandidates) 
        {
            if (forcedNextOrNeg1 >= 0 && j != forcedNextOrNeg1) continue;

            // Težina za čvor j (kretanje po turi): alphaNodes, betaNodes
            double tauNode = 0.0;
            double minCost = std::numeric_limits<double>::max();
            for (int c = 0; c < C; ++c) 
            {
                if (enforceNoRerent_ && c != prevCarState && carLocked[(size_t)c]) continue;
                tauNode += ctx.pher->tauMove(c, currentNode, j);
                minCost = std::min(minCost, inst.travelCost(c, currentNode, j));
            }
            const double etaNode = eta(eps_, minCost);
            const double W_node = std::pow(tauNode, alphaNodes) * std::pow(etaNode, betaNodes);
            if (!(W_node > 0.0) || !std::isfinite(W_node)) continue;

            for (int c = 0; c < C; ++c) 
            {
                if (enforceNoRerent_) 
                {
                    if (c != prevCarState && carLocked[(size_t)c]) continue;
                }
                // Težina za auto na bridu (i->j): alphaCars, betaCars
                const double tMove = ctx.pher->tauMove(c, currentNode, j);
                const double hMove = eta(eps_, inst.travelCost(c, currentNode, j));
                double w = W_node * (std::pow(tMove, alphaCars) * std::pow(hMove, betaCars));

                if (!(w > 0.0) || !std::isfinite(w)) continue;

                if (useReturnPolicy)
                {
                    if (prevCarState >= 0 && c != prevCarState) 
                    {
                        const double rc = inst.returnCost(prevCarState, currentNode, rentalNodeState);
                        const double tR = ctx.pher->tauReturn(prevCarState, currentNode, rentalNodeState);
                        const double hR = eta(eps_, rc);

                        // multiplicativno (minimalno invazivno): pojača/oslabi switch
                        w *= std::pow(tR, alphaR) * std::pow(hR, betaR);
                    }

                    // (2) Ako je ovo zadnji korak (lastNode -> 0), uključi i final return
                    if (isLastStepToStart) 
                    {
                        // rentalNode nakon mogućeg switcha:
                        const int rentalAfter = (prevCarState >= 0 && c != prevCarState)
                                                ? currentNode
                                                : rentalNodeState;

                        const double rcFinal = inst.returnCost(c, /*start*/ 0, rentalAfter);
                        const double tFinal  = ctx.pher->tauReturn(c, /*start*/ 0, rentalAfter);
                        const double hFinal  = eta(eps_, rcFinal);

                        w *= std::pow(tFinal, alphaR) * std::pow(hFinal, betaR);
                    }
                }

                if (!(w > 0.0) || !std::isfinite(w)) continue;

                items.push_back({j, c, w});
                sum += w;
            }
        }

        // Fallback: uniformno
        if (items.empty() || !(sum > 0.0) || !std::isfinite(sum)) 
        {
            int j = (forcedNextOrNeg1 >= 0)
                ? forcedNextOrNeg1
                : nextCandidates[rng.uniformInt((int)nextCandidates.size())];
            int c = rng.uniformInt(C);
            return {j, c};
        }

        // exploitation jump (argmax over (j,c)) ---
        // OPTIMIZACIJA: Koristi dinamički q0 iz ctx.adaptiveQ0 ako je dostupan
        const double q0 = (ctx.adaptiveQ0 >= 0.0) ? ctx.adaptiveQ0 : q0_;
        if (q0 > 0.0 && rng.uniform01() < q0) {
            size_t bestIdx = 0;
            double bestW = items[0].w;
            for (size_t t = 1; t < items.size(); ++t) {
                if (items[t].w > bestW) { bestW = items[t].w; bestIdx = t; }
            }
            return {items[bestIdx].j, items[bestIdx].c};
        }

        double r = rng.uniform01() * sum;
        double acc = 0.0;
        for (const auto& it : items) 
        {
            acc += it.w;
            if (acc >= r) return {it.j, it.c};
        }
        return {items.back().j, items.back().c};
    }; // choose_from_pairs
    // Napomena:
    // - sampleCars_==true: policy uzorkuje (nextNode, car) parove (default i najbolji za Brasil30n CaRS).
    // - sampleCars_==false: node-only mod (ablacija); DP kasnije dodjeljuje aute.


    // Glavna konstrukcija ture (pos = indeks brida/odluke)
    const int fav = ctx.params->favorites;
    const int maxCand = (fav > 0) ? std::min(fav, N - 1) : 0;

    for (int pos = 0; pos < N - 1; ++pos) 
    {
        const int i = s.node[(size_t)pos];

        std::vector<int> cand;
        cand.reserve((size_t)(maxCand > 0 ? maxCand : N - 1));

        if (cl_ && cl_->N() == N) 
        {
            for (int j : cl_->candidates(i)) 
            {
                if (!visited[(size_t)j]) cand.push_back(j);
                if (maxCand > 0 && (int)cand.size() >= maxCand) break;
            }
        }
        if (cand.empty()) 
        {
            for (int j = 0; j < N; ++j)
                if (!visited[(size_t)j]) cand.push_back(j);
        }

        auto [nextNode, chosenCar] =
            choose_from_pairs(i, cand, -1, prevCar, rentalNode, /*isLastStepToStart*/ false);

        // Postavi auto za brid i->nextNode
        s.car[(size_t)pos] = chosenCar;

        // Ažuriraj stanje switcha (ekvivalent evaluate() logike)
        // "no re-rent" je neovisan od returnFactor-a (možeš ga htjeti i kad returnFactor=0),
        // ali nema smisla dirati carLocked ako feature nije uključen.
        if (enforceNoRerent_ && prevCar >= 0 && chosenCar != prevCar) {
            carLocked[(size_t)prevCar] = 1;
        }

        // rentalNode je potreban samo ako policy stvarno koristi return logiku.
        if (useReturnPolicy && prevCar >= 0 && chosenCar != prevCar) {
            rentalNode = i;
        }

prevCar = chosenCar;



        s.node[(size_t)pos + 1] = nextNode;
        visited[(size_t)nextNode] = 1;
    }

    // Zadnji brid: lastNode -> start (0), uz moguće switch + final return u težini
    const int lastNode = s.node[(size_t)N - 1];
    {
        std::vector<int> onlyStart = {0};
        auto [_, chosenCar] =
            choose_from_pairs(lastNode, onlyStart, 0, prevCar, rentalNode, /*isLastStepToStart*/ true);
        (void)_;
        s.car[(size_t)N - 1] = chosenCar;

        if (enforceNoRerent_ && prevCar >= 0 && chosenCar != prevCar) 
        {
            carLocked[(size_t)prevCar] = 1;
        }
        if (useReturnPolicy && prevCar >= 0 && chosenCar != prevCar) 
        {
            rentalNode = lastNode;
        }

    prevCar = chosenCar;


    }

    return s;
}

} // namespace aco
