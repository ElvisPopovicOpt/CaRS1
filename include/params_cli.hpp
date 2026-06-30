#pragma once

#include <stdexcept>
#include <string>
#include <iosfwd>

namespace aco_cli
{

// Defaults (as requested)
#define DEFAULT_FILENAME "default/Afe4ns.car"
#define MAX_NRUNS 101
#define DEFAULT_NRUNS 1
#define DEFAULT_ITERATIONS 1000
#define DEFAULT_ANT_NUMBER 100
#define DEFAULT_ALPHA 1.5
#define DEFAULT_BETA 1.1
#define DEFAULT_RHO 0.1
#define DEFAULT_MAXMIN 0 // if zero, no maxmin
#define DEFAULT_STAGNATION 30
#define DEFAULT_ELITE_K_ANTS 1   // top-1 = iterBest only
#define DEFAULT_FAVORITES 20  // Recommended: 0.5N to 0.67N for good exploration/exploitation balance
#define DEFAULT_GB_PERIOD 0 //0 auto, >0 fixed, -1 off

// New defaults (recommended)
#define DEFAULT_SMOOTHING_GAMMA 0.7  // Povećano s 0.5 na 0.7 za jači reset feromona u kasnoj stagnaciji
#define DEFAULT_EXPLORE_ENABLED 0
#define DEFAULT_EXPLORE_ITERS 20
#define DEFAULT_EXPLORE_RHO_MULT 2.0
#define DEFAULT_EPSILON_ETA 1e-6

#define DEFAULT_POLICY_Q0 0.1
#define DEFAULT_ADAPTIVE_Q0 1  // 0 = disabled, 1 = enabled (default: enabled za bržu konvergenciju)
#define DEFAULT_Q0_START 0.2    // početna vrijednost q0 (više eksploatacije od početka)
#define DEFAULT_Q0_END 0.9      // konačna vrijednost q0 (eksploatacija)

// New defaults (requested)
#define DEFAULT_RETURN_FACTOR_CARS 0.0
#define DEFAULT_SAMPLE_CARS 1
#define DEFAULT_ENFORCE_NO_RERENT 1 // zabranjeno ponovno iznajmljivanje

// Local search and DP tuning
#define DEFAULT_LS_TOP_W 2  // number of elite solutions to improve with LS per iteration
#define DEFAULT_DP_CACHE_CAPACITY 5000
#define DEFAULT_DP_MAX_SEGMENT_LEN_OFFSET 5  // offset from ceil(N/C)
#define DEFAULT_USE_FULL_POLISH 0  // 0 = lite (brži), 1 = full (sporiji ali može pronaći bolja rješenja)
#define DEFAULT_LS_QUALITY_MODE 3  // 1=low, 2=medium, 3=high, 4=extreme

// Intensifier tuning
#define DEFAULT_INTENSIFIER_MAX_STAGNATION_ITERATIONS 3  // Broj iteracija nakon stagnacije kada se još aktivira intenzifikator
#define DEFAULT_INTENSIFIER_RESET_PHEROMONES_COUNT 1  // Broj deaktivacija nakon kojih se resetiraju feromoni (0 ili negativno = nikad, 1 = nakon 1. gašenja, 2 = nakon 1. i 2., itd.)
#define DEFAULT_INTENSIFIER_PHEROMONE_REDUCTION_GAMMA 0.3  // Faktor smanjenja feromona na best tour-u (0.0=potpuno resetiranje, 1.0=bez promjene, default: 0.3)
#define DEFAULT_INTENSIFIER_ARCHIVE_AND_RESET_BEST 1  // Arhiviraj global best i resetiraj ga nakon prve deaktivacije (0=disabled, 1=enabled, default: enabled) 
#define DEFAULT_ABLATION_DISABLE_ANT_DEPOSITS_DURING_INTENSIFICATION 0
#define DEFAULT_ABLATION_RESET_PHEROMONES_ON_INTENSIFIER_ACTIVATION 0
#define DEFAULT_ABLATION_INTENSIFIER_ONLY_NO_ANTS 0
#define DEFAULT_ABLATION_INTENSIFIER_ONLY_RANDOM_INITIAL_TOUR 1



// Optional (if you add MMAS tauMin from pBest
// Bude obicno 0.05 a 0 je iskljuceno koristi se maxMin)
#define DEFAULT_PBEST 0.0

enum class RestartTarget
{
    Mid,    // (tauMin + tauMax)/2  (trenutno ponašanje)
    TauMax  // smoothing prema tauMax (često bolji MMAS restart)
    // TauMin (možeš dodati ako ikad zatreba)
};

static inline const char* toString(RestartTarget x)
{
    switch (x) {
        case RestartTarget::Mid:   return "mid";
        case RestartTarget::TauMax:return "tauMax";
    }
    return "mid";
}

struct ParamsData
{
    std::string filename;
    long long seed = 0;

    int nRuns = DEFAULT_NRUNS;

    double maxMinCars = DEFAULT_MAXMIN;
    double maxMinNodes = DEFAULT_MAXMIN;

    int iterations = DEFAULT_ITERATIONS;
    int antsN = DEFAULT_ANT_NUMBER;

    double rhoNodes = DEFAULT_RHO;
    double rhoCars = DEFAULT_RHO;
    double rhoPass = DEFAULT_RHO;

    double alphaNodes = DEFAULT_ALPHA;
    double alphaCars = DEFAULT_ALPHA;
    double alphaPass = DEFAULT_ALPHA;

    double betaNodes = DEFAULT_BETA;
    double betaCars = DEFAULT_BETA;
    double betaPass = DEFAULT_BETA;

    int stagnation = DEFAULT_STAGNATION;

    int eliteKAnts = DEFAULT_ELITE_K_ANTS;
    int favorites = DEFAULT_FAVORITES;
    int gbPeriod = DEFAULT_GB_PERIOD;

    // --- Cars: stagnation / restart tuning
    double smoothingGammaCars = DEFAULT_SMOOTHING_GAMMA;
    RestartTarget restartTargetCars = RestartTarget::Mid;

    // --- Cars: exploration phase after stagnation
    bool exploreEnabledCars = (DEFAULT_EXPLORE_ENABLED != 0);
    int exploreItersCars = DEFAULT_EXPLORE_ITERS;
    double exploreRhoMultiplierCars = DEFAULT_EXPLORE_RHO_MULT;

    double policyQ0 = DEFAULT_POLICY_Q0;
    
    // Dinamički q0 parametri
    bool adaptiveQ0 = (DEFAULT_ADAPTIVE_Q0 != 0);  // omogući dinamički q0
    double q0Start = DEFAULT_Q0_START;              // početna vrijednost (eksploracija)
    double q0End = DEFAULT_Q0_END;                  // konačna vrijednost (eksploatacija)

    // --- Cars: heuristic eta=1/(eps+cost)
    double epsilonEtaCars = DEFAULT_EPSILON_ETA;

    // --- Optional: MMAS p_best (if you later implement tauMin via pBest)
    double pBestCars = DEFAULT_PBEST;

        // --- NEW: return policy influence (cars)
    double returnFactorCars = DEFAULT_RETURN_FACTOR_CARS; // 0 disables return logic in policy
    // --- NEW: sampling mode
    bool sampleCars = (DEFAULT_SAMPLE_CARS != 0); // true: sample (nextNode, car), false: only node permutation
    // --- NEW: no re-rent constraint during construction
    bool enforceNoRerent = (DEFAULT_ENFORCE_NO_RERENT != 0);

    // --- Local search and DP tuning
    int lsTopW = DEFAULT_LS_TOP_W;  // number of elite solutions to improve with LS per iteration
    int dpCacheCapacity = DEFAULT_DP_CACHE_CAPACITY;  // DP cache size (0 = disabled)
    int dpMaxSegmentLenOffset = DEFAULT_DP_MAX_SEGMENT_LEN_OFFSET;  // DP maxSegmentLen = ceil(N/C) + offset
    bool useFullPolish = (DEFAULT_USE_FULL_POLISH != 0);  // 0 = lite (brži), 1 = full (sporiji ali može pronaći bolja rješenja)
    
    // --- LS Quality parameters
    int lsQualityMode = DEFAULT_LS_QUALITY_MODE;  // 1=low, 2=medium, 3=high, 4=extreme
    int lsTwoOpt1Passes = -1;  // -1 = use mode, >0 = override
    int lsRelocPasses = -1;    // -1 = use mode, >0 = override
    int lsTwoOpt2Passes = -1;  // -1 = use mode, >0 = override
    int lsDpRetries = -1;      // -1 = use mode, >0 = override

    // --- Reinsertion Local Search (single-node relocation in LS chain: 2-opt -> Reloc -> Reinsert -> 2-opt)
    bool useReinsertion = false;   // 0 = isključeno (default), 1 = uključeno
    int reinsertionMaxPasses = 5;  // max pass-ova za Reinsertion kad je uključen
    
    // --- Intensifier parameters
    bool intensifierEnabled = true;  // 0 = isključi intenzifikator u potpunosti (ne upliće se u proračun ture)
    int intensifierMaxStagnationIterations = DEFAULT_INTENSIFIER_MAX_STAGNATION_ITERATIONS;  // Broj iteracija nakon stagnacije kada se još aktivira intenzifikator
    int intensifierResetPheromonesCount = DEFAULT_INTENSIFIER_RESET_PHEROMONES_COUNT;  // Broj deaktivacija nakon kojih se resetiraju feromoni (0 ili negativno = nikad, 1 = nakon 1. gašenja, 2 = nakon 1. i 2., itd.)
    double intensifierPheromoneReductionGamma = DEFAULT_INTENSIFIER_PHEROMONE_REDUCTION_GAMMA;  // Faktor smanjenja feromona na best tour-u (selektivno resetiranje)
    bool intensifierArchiveAndResetBest = (DEFAULT_INTENSIFIER_ARCHIVE_AND_RESET_BEST != 0);  // Arhiviraj global best i resetiraj ga nakon prve deaktivacije
    
    // --- Ablation: tijekom aktivne intenzifikacije ugasi depoziciju mrava (ranked top-W)
    bool ablationDisableAntDepositsDuringIntensification =
        (DEFAULT_ABLATION_DISABLE_ANT_DEPOSITS_DURING_INTENSIFICATION != 0);
    // --- Ablation: na svaku aktivaciju intenzifikatora resetiraj feromone
    bool ablationResetPheromonesOnIntensifierActivation =
        (DEFAULT_ABLATION_RESET_PHEROMONES_ON_INTENSIFIER_ACTIVATION != 0);
    // --- Ablation: bez mrava — samo početna tura (random ili 0,1,..,N-1) + isti ostatak petlje (LS, intenzifikator, …)
    bool ablationIntensifierOnlyNoAnts = (DEFAULT_ABLATION_INTENSIFIER_ONLY_NO_ANTS != 0);
    bool ablationIntensifierOnlyRandomInitialTour =
        (DEFAULT_ABLATION_INTENSIFIER_ONLY_RANDOM_INITIAL_TOUR != 0);

    // --- Surrogate poboljšanje (korelacija s punom DP cijenom)
    double surrogateReturnGamma = 0.5;  // 0 = isključeno; >0 = dodaj procjenu return costova u surogat (gamma * estReturns * avgReturnCost). Preporuka 0.3–0.5.

    // --- Research / statistical logging (binomna analiza intenzifikatora, ablacije)
    bool researchLog = false;  // 1 = snimi intensifier_binomial.csv (run, activation_index, cost_start, cost_end, improved, depth)
    bool researchSurrogateLog = false;  // 1 = snimi surrogate_correlation.csv (run, iter, ant_index, surrogate, dp_cost) za analizu korelacije surogata s punim DP costom

    // --- Pheromone Model Selection
    bool useTBAS = false;  // 0 = MMAS (default), 1 = TBAS (Three Bounds Ant System)
};



struct ParseError : std::runtime_error
{
    using std::runtime_error::runtime_error;
};

struct ParserOptions
{
    bool strictUnknownParams = true; // unknown -> error (true) or ignore (false)
    bool printHelpOnEmpty = true;    // if argc==1 => behave like --help
};

class Parser
{
public:
    explicit Parser(ParserOptions opt = {});

    // Parses argc/argv and returns filled ParamsData (with defaults).
    // Throws ParseError on invalid input.
    // Throws ParseError("HELP_REQUESTED") if help requested or (opt.printHelpOnEmpty && argc==1).
    ParamsData parse(int argc, const char* const* argv) const;

    void printUsage(std::ostream& os, const std::string& exeName = "program") const;

private:
    ParserOptions opt_;
};

} // namespace aco_cli
