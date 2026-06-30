#pragma once
#include <solution.hpp>
#include <rng.hpp>  
#include <mutex>
#include <memory>
#include <limits>

namespace cars_tsplib { struct Instance; }
namespace aco_cli { struct ParamsData; }

namespace aco 
{

struct IPheromoneAccessor
{
    virtual ~IPheromoneAccessor() = default;

    // MOVE pheromoni (postoje uvijek)
    virtual double tauMove(int car, int i, int j) const = 0;

    // RETURN pheromoni (opcionalno).
    // Default: nema return feromona -> neutralno 1.0 (ne utječe na težine).
    virtual double tauReturn(int /*car*/, int /*from*/, int /*to*/) const { return 1.0; }

    // Pomoćna informacija (opcionalno, ali korisno za ant-policy):
    virtual bool hasReturnPheromones() const { return false; }
};

class UniformPheromoneAccessor final : public IPheromoneAccessor
{
public:
    explicit UniformPheromoneAccessor(double tau = 1.0) : tau_(tau) {}
    double tauMove(int, int, int) const override { return tau_; }

    // uniform accessor nema return feromone
    double tauReturn(int, int, int) const override { return 1.0; }
    bool hasReturnPheromones() const override { return false; }

private:
    double tau_;
};


struct IPheromoneModel : public IPheromoneAccessor
{
    virtual ~IPheromoneModel() = default;

    virtual void reset() = 0;

    virtual void onIterationEnd(int iterationIndex,
                                const EvaluatedSolution& iterBest,
                                const EvaluatedSolution& globalBest,
                                bool globalImproved) = 0;

    // NOVO: ranked (top-W) update. Default: koristi samo ranked<a href="" class="citation-link" target="_blank" style="vertical-align: super; font-size: 0.8em; margin-left: 3px;">[0]</a> kao iterBest.
    // Colony će W zasad vezati uz eliteKAnts (kasnije se može dodati poseban parametar).
    virtual void onIterationEndRanked(int iterationIndex,
                                      const std::vector<EvaluatedSolution>& ranked,
                                      const EvaluatedSolution& globalBest,
                                      bool globalImproved)
    {
        if (ranked.empty()) return;
        onIterationEnd(iterationIndex, ranked.front(), globalBest, globalImproved);
    }

    // Check if smoothing was applied in the last iteration (default: false)
    virtual bool wasSmoothingApplied() const { return false; }
    
    // Get stagnation counter (noImprove_) for smoothing calculation (default: 0)
    virtual int getStagnationCounter() const { return 0; }
    
    /**
     * Selektivno smanjuje feromone na određenom tour-u (smoothing).
     * Umjesto potpunog resetiranja, smanjuje feromone na edges u tour-u s faktorom gamma.
     * Formula: tau = (1 - gamma) * tau + gamma * tauBase
     * 
     * @param solution Rješenje čiji tour treba smanjiti feromone
     * @param gamma Faktor smanjenja (0.0 = potpuno resetiranje na tauBase, 1.0 = bez promjene)
     */
    virtual void reducePheromonesOnTour(const Solution& solution, double gamma) { (void)solution; (void)gamma; }
    
    /**
     * Dodaje rješenje u archive (ako je archive omogućen).
     * 
     * @param solution Rješenje koje treba dodati u archive
     */
    virtual void addToArchive(const EvaluatedSolution& solution) { (void)solution; }
    
    /**
     * Vraća najbolji arhivirani best (ili infinity ako archive nije omogućen ili je prazan).
     * 
     * @return Najbolji cost iz archive-a, ili infinity ako nema arhiviranih rješenja
     */
    virtual double getBestArchivedCost() const { return std::numeric_limits<double>::infinity(); }
    
    /**
     * Vraća najbolje arhivirano rješenje (ili prazno rješenje ako archive nije omogućen ili je prazan).
     * 
     * @param out Output parametar za najbolje arhivirano rješenje
     * @return true ako postoji arhivirano rješenje, false inače
     */
    virtual bool getBestArchivedSolution(EvaluatedSolution& out) const { (void)out; return false; }
    
};


struct ICostModel 
{
    virtual ~ICostModel() = default;
    virtual double evaluate(const Solution& s) const = 0;
};

struct ILocalSearch 
{
    virtual ~ILocalSearch() = default;
    virtual void improve(Solution& s, double& cost) const = 0;
    // opcionalno seedanje (default: ništa)
    virtual void resetSeed(uint64_t /*seed*/) const {}
};



struct AntContext 
{
    std::shared_ptr<const cars_tsplib::Instance> inst;
    const aco_cli::ParamsData* params = nullptr;
    const IPheromoneAccessor* pher = nullptr;

    // per-run RNG 
    Rng* rng = nullptr;
    
    // Dinamički q0 override: -1.0 = koristi fiksni q0 iz ant policy-ja, >=0.0 = override q0
    double adaptiveQ0 = -1.0;
};

struct IAntPolicy 
{
    virtual ~IAntPolicy() = default;
    virtual Solution construct(const AntContext& ctx) const = 0;
};

struct IFixedTourCarAssigner
{
    virtual ~IFixedTourCarAssigner() = default;

    // za zadani Hamiltonov ciklus (nodes<a href="" class="citation-link" target="_blank" style="vertical-align: super; font-size: 0.8em; margin-left: 3px;">[0]</a>==0, size N)
    // popuni carPerEdgeOut (size N) i vrati minimalni cost
    virtual double reassignCars(const std::vector<int>& nodes,
                                std::vector<int>& carPerEdgeOut) const = 0;

    // cost-only (opcionalno, ali korisno)
    virtual double evaluateCost(const std::vector<int>& nodes) const = 0;
};

/**
 * Interface za intenzifikator komponentu.
 * 
 * Intenzifikator primjenjuje dodatni Local Search ciklus na global best rješenje
 * za dublju eksploataciju bez prevelikog utroška vremena.
 * 
 * Aktivira se kada:
 * - Best je upravo promijenjen u ovoj iteraciji (improvedByIter == true)
 * - ILI smo u prvih nekoliko iteracija nakon stagnacije
 */
struct IIntensifier
{
    virtual ~IIntensifier() = default;
    
    /**
     * Primjenjuje intenzifikaciju na best rješenje.
     * 
     * @param best Global best rješenje (može biti ažurirano ako postoji poboljšanje - pass-by-reference)
     * @param improvedByIter Je li best upravo promijenjen u ovoj iteraciji
     * @param itersSinceGlobalImprovement Broj iteracija od zadnjeg globalnog poboljšanja
     * @param improveEps Epsilon za provjeru poboljšanja (numerička preciznost)
     * @return Cost poboljšanja (ili infinity ako nije bilo poboljšanja)
     */
    virtual double intensify(EvaluatedSolution& best,
                            bool improvedByIter,
                            int itersSinceGlobalImprovement,
                            double improveEps) const = 0;
    
    /**
     * Provjerava treba li resetirati feromone.
     * Vraća true ako je intenzifikator upravo deaktiviran (bio aktivan, sada nije),
     * i brojač aktivacija <= parametru (npr. parametar 1 = reset nakon 1. gašenja, 2 = nakon 1. i 2.).
     * Parametar 0 ili negativno = nikad resetirati.
     * 
     * @return True ako treba resetirati feromone, false inače
     */
    virtual bool shouldResetPheromones() const = 0;
    
    /**
     * Provjerava je li intenzifikator upravo deaktiviran (bio aktivan, sada nije).
     * 
     * @return True ako je intenzifikator upravo deaktiviran, false inače
     */
    virtual bool wasJustDeactivated() const = 0;
    
    /**
     * Vraća broj puta koliko je intenzifikator bio uključen.
     * 
     * @return Broj aktivacija intenzifikatora
     */
    virtual int getActivationCount() const = 0;
    
    /**
     * Ažurira interno stanje za sljedeći poziv.
     * Treba se pozvati nakon što se provjeri deaktivacija i resetiranje feromona.
     */
    virtual void updateState() const = 0;
};

} // namespace aco
