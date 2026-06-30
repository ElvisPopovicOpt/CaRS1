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
#define DEFAULT_SMOOTHING_GAMMA 0.7  // stronger pheromone reset in late stagnation
#define DEFAULT_EXPLORE_ENABLED 0
#define DEFAULT_EXPLORE_ITERS 20
#define DEFAULT_EXPLORE_RHO_MULT 2.0
#define DEFAULT_EPSILON_ETA 1e-6

#define DEFAULT_POLICY_Q0 0.1
#define DEFAULT_ADAPTIVE_Q0 1  // 0 = disabled, 1 = enabled (default: enabled, for faster convergence)
#define DEFAULT_Q0_START 0.2    // starting q0 value (more exploitation from the start)
#define DEFAULT_Q0_END 0.9      // final q0 value (exploitation)

// New defaults (requested)
#define DEFAULT_RETURN_FACTOR_CARS 0.0
#define DEFAULT_SAMPLE_CARS 1
#define DEFAULT_ENFORCE_NO_RERENT 1 // re-renting after return is disallowed

// Local search and DP tuning
#define DEFAULT_LS_TOP_W 2  // number of elite solutions to improve with LS per iteration
#define DEFAULT_DP_CACHE_CAPACITY 5000
#define DEFAULT_DP_MAX_SEGMENT_LEN_OFFSET 5  // offset from ceil(N/C)
#define DEFAULT_USE_FULL_POLISH 0  // 0 = lite (faster), 1 = full (slower but may find better solutions)
#define DEFAULT_LS_QUALITY_MODE 3  // 1=low, 2=medium, 3=high, 4=extreme

// Intensifier tuning
#define DEFAULT_INTENSIFIER_MAX_STAGNATION_ITERATIONS 3  // iterations after stagnation during which the intensifier stays active
#define DEFAULT_INTENSIFIER_RESET_PHEROMONES_COUNT 1  // number of deactivations that trigger a pheromone reset (0 or negative = never, 1 = after 1st, 2 = after 1st and 2nd, etc.)
#define DEFAULT_INTENSIFIER_PHEROMONE_REDUCTION_GAMMA 0.3  // pheromone reduction factor on the best tour (0.0 = full reset, 1.0 = no change)
#define DEFAULT_INTENSIFIER_ARCHIVE_AND_RESET_BEST 1  // archive global best and reset it after the first deactivation (0=disabled, 1=enabled)
#define DEFAULT_ABLATION_DISABLE_ANT_DEPOSITS_DURING_INTENSIFICATION 0
#define DEFAULT_ABLATION_RESET_PHEROMONES_ON_INTENSIFIER_ACTIVATION 0
#define DEFAULT_ABLATION_INTENSIFIER_ONLY_NO_ANTS 0
#define DEFAULT_ABLATION_INTENSIFIER_ONLY_RANDOM_INITIAL_TOUR 1



// Optional MMAS tauMin via pBest. Typically 0.05; 0 disables it and maxMin is used instead.
#define DEFAULT_PBEST 0.0

enum class RestartTarget
{
    Mid,    // (tauMin + tauMax)/2 (current behavior)
    TauMax  // smoothing toward tauMax (often a better MMAS restart)
    // TauMin could be added later if needed
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
    
    // Adaptive q0 parameters
    bool adaptiveQ0 = (DEFAULT_ADAPTIVE_Q0 != 0);  // enables dynamic q0
    double q0Start = DEFAULT_Q0_START;              // starting value (exploration)
    double q0End = DEFAULT_Q0_END;                  // final value (exploitation)

    // --- Cars: heuristic eta=1/(eps+cost)
    double epsilonEtaCars = DEFAULT_EPSILON_ETA;

    // --- Optional: MMAS p_best (if you later implement tauMin via pBest)
    double pBestCars = DEFAULT_PBEST;

    // --- Return policy influence (cars)
    double returnFactorCars = DEFAULT_RETURN_FACTOR_CARS; // 0 disables return logic in policy
    // --- Sampling mode
    bool sampleCars = (DEFAULT_SAMPLE_CARS != 0); // true: sample (nextNode, car), false: only node permutation
    // --- No re-rent constraint during construction
    bool enforceNoRerent = (DEFAULT_ENFORCE_NO_RERENT != 0);

    // --- Local search and DP tuning
    int lsTopW = DEFAULT_LS_TOP_W;  // number of elite solutions to improve with LS per iteration
    int dpCacheCapacity = DEFAULT_DP_CACHE_CAPACITY;  // DP cache size (0 = disabled)
    int dpMaxSegmentLenOffset = DEFAULT_DP_MAX_SEGMENT_LEN_OFFSET;  // DP maxSegmentLen = ceil(N/C) + offset
    bool useFullPolish = (DEFAULT_USE_FULL_POLISH != 0);  // 0 = lite (faster), 1 = full (slower but may find better solutions)

    // --- LS Quality parameters
    int lsQualityMode = DEFAULT_LS_QUALITY_MODE;  // 1=low, 2=medium, 3=high, 4=extreme
    int lsTwoOpt1Passes = -1;  // -1 = use mode, >0 = override
    int lsRelocPasses = -1;    // -1 = use mode, >0 = override
    int lsTwoOpt2Passes = -1;  // -1 = use mode, >0 = override
    int lsDpRetries = -1;      // -1 = use mode, >0 = override

    // --- Reinsertion local search (single-node relocation in LS chain: 2-opt -> Reloc -> Reinsert -> 2-opt)
    bool useReinsertion = false;   // 0 = disabled (default), 1 = enabled
    int reinsertionMaxPasses = 5;  // max passes for Reinsertion when enabled

    // --- Intensifier parameters
    bool intensifierEnabled = true;  // 0 = fully disable the intensifier (no effect on tour computation)
    int intensifierMaxStagnationIterations = DEFAULT_INTENSIFIER_MAX_STAGNATION_ITERATIONS;  // iterations after stagnation during which the intensifier stays active
    int intensifierResetPheromonesCount = DEFAULT_INTENSIFIER_RESET_PHEROMONES_COUNT;  // number of deactivations that trigger a pheromone reset (0 or negative = never, 1 = after 1st, 2 = after 1st and 2nd, etc.)
    double intensifierPheromoneReductionGamma = DEFAULT_INTENSIFIER_PHEROMONE_REDUCTION_GAMMA;  // selective pheromone reduction factor on the best tour
    bool intensifierArchiveAndResetBest = (DEFAULT_INTENSIFIER_ARCHIVE_AND_RESET_BEST != 0);  // archive global best and reset it after the first deactivation

    // --- Ablation: disable ant deposits (ranked top-W) while the intensifier is active
    bool ablationDisableAntDepositsDuringIntensification =
        (DEFAULT_ABLATION_DISABLE_ANT_DEPOSITS_DURING_INTENSIFICATION != 0);
    // --- Ablation: reset pheromones on every intensifier activation
    bool ablationResetPheromonesOnIntensifierActivation =
        (DEFAULT_ABLATION_RESET_PHEROMONES_ON_INTENSIFIER_ACTIVATION != 0);
    // --- Ablation: no ants — only an initial tour (random or 0,1,..,N-1) plus the rest of the loop (LS, intensifier, ...)
    bool ablationIntensifierOnlyNoAnts = (DEFAULT_ABLATION_INTENSIFIER_ONLY_NO_ANTS != 0);
    bool ablationIntensifierOnlyRandomInitialTour =
        (DEFAULT_ABLATION_INTENSIFIER_ONLY_RANDOM_INITIAL_TOUR != 0);

    // --- Surrogate improvement (correlation with full DP cost)
    double surrogateReturnGamma = 0.5;  // 0 = off; >0 = add estimated return cost to the surrogate (gamma * estReturns * avgReturnCost). Recommended 0.3-0.5.

    // --- Research / statistical logging (intensifier binomial analysis, ablations)
    bool researchLog = false;  // 1 = write intensifier_binomial.csv (run, activation_index, cost_start, cost_end, improved, depth)
    bool researchSurrogateLog = false;  // 1 = write surrogate_correlation.csv (run, iter, ant_index, surrogate, dp_cost) to analyze surrogate vs full DP cost correlation

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
