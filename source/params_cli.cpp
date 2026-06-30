#include <params_cli.hpp>

#include <chrono>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <map>


namespace aco_cli
{
    
namespace
{

static inline long long defaultSeedFromTime()
{
    using namespace std::chrono;
    const milliseconds ms = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
    return static_cast<long long>(ms.count());
}

static inline bool startsWith(const std::string& s, const char* prefix)
{
    const std::string p(prefix);
    return s.size() >= p.size() && s.compare(0, p.size(), p) == 0;
}

static inline std::pair<std::string, std::string> splitKeyEqValue(const std::string& token)
{
    const std::size_t eq = token.find('=');
    if (eq == std::string::npos) return { token, {} };
    return { token.substr(0, eq), token.substr(eq + 1) };
}

static inline std::string stripLeadingDashes(std::string s)
{
    while (!s.empty() && s.front() == '-') s.erase(s.begin());
    return s;
}

static inline long long parseLL(const std::string& s, const std::string& what)
{
    std::size_t pos = 0;
    long long v = 0;
    try {
        v = std::stoll(s, &pos, 10);
    } catch (...) {
        throw ParseError("Invalid " + what + " (expected integer): '" + s + "'");
    }
    if (pos != s.size()) throw ParseError("Invalid " + what + " (trailing chars): '" + s + "'");
    return v;
}

static inline int parseInt(const std::string& s, const std::string& what)
{
    const long long v = parseLL(s, what);
    if (v < std::numeric_limits<int>::min() || v > std::numeric_limits<int>::max())
        throw ParseError("Invalid " + what + " (out of int range): '" + s + "'");
    return static_cast<int>(v);
}

static inline double parseDouble(const std::string& s, const std::string& what)
{
    std::size_t pos = 0;
    double v = 0.0;
    try {
        v = std::stod(s, &pos);
    } catch (...) {
        throw ParseError("Invalid " + what + " (expected floating point): '" + s + "'");
    }
    if (pos != s.size()) throw ParseError("Invalid " + what + " (trailing chars): '" + s + "'");
    return v;
}

static inline bool parseBool01(const std::string& v, const char* name)
{
    std::string s;
    s.reserve(v.size());
    for (char c : v) s.push_back((char)std::tolower((unsigned char)c));

    if (s == "1" || s == "true" || s == "yes" || s == "on")  return true;
    if (s == "0" || s == "false" || s == "no"  || s == "off") return false;

    throw ParseError(std::string(name) + ": expected bool (0/1, true/false, yes/no, on/off), got '" + v + "'");
}

static inline aco_cli::RestartTarget parseRestartTarget(const std::string& v, const char* name)
{
    std::string s;
    s.reserve(v.size());
    for (char c : v) s.push_back((char)std::tolower((unsigned char)c));

    if (s == "mid") return aco_cli::RestartTarget::Mid;
    if (s == "taumax" || s == "max") return aco_cli::RestartTarget::TauMax;

    throw ParseError(std::string(name) + ": expected 'mid' or 'tauMax', got '" + v + "'");
}

struct OptionSpec
{
    const char* longName;     // e.g. "carAlpha"
    const char* shortName;    // e.g. "ca" (used as -ca). nullptr if none.
    const char* typeHint;     // e.g. "double"
    const char* defaultHint;  // e.g. "1.1"
    const char* help;         // short text
    const char* category;     // e.g. "Quality", "Pheromone", "Local Search", etc. nullptr for no category
    std::function<void(ParamsData&, const std::string&)> set;
};

static inline std::vector<OptionSpec> makeSpecs()
{
    // CENTRAL PLACE TO ADD NEW PARAMS:
    // Add one OptionSpec line here.

    std::vector<OptionSpec> s;

    // === Basics ===
    s.push_back({ "filename",  "f",  "string",   DEFAULT_FILENAME, "Input instance file", "Basics",
        [](ParamsData& d, const std::string& v){ d.filename = v; } });

    s.push_back({ "seed",      "s",  "longlong", "time", "Random seed (default: from clock)", "Basics",
        [](ParamsData& d, const std::string& v){ d.seed = parseLL(v, "seed"); } });

    // === Algorithm ===
    s.push_back({ "nRuns",     "n",  "int",      "1",      "Number of runs", "Algorithm",
        [](ParamsData& d, const std::string& v){ d.nRuns = parseInt(v, "nRuns"); } });

    s.push_back({ "iterations",    "it",  "int",      "1000",  "Max iterations", "Algorithm",
        [](ParamsData& d, const std::string& v){ d.iterations = parseInt(v, "iterations"); } });

    s.push_back({ "antsN",     "k",  "int",      "100",    "Number of ants", "Algorithm",
        [](ParamsData& d, const std::string& v){ d.antsN = parseInt(v, "antsN"); } });

    s.push_back({ "stagnation","st", "int",      "30",     "Stagnation limit: iterations without global improvement before polish/kick. 0 = disabled. Recommended: 20-50.", "Algorithm",
        [](ParamsData& d, const std::string& v){ d.stagnation = parseInt(v, "stagnation"); } });

    s.push_back({ "eliteKAnts", "eka", "int", "1",
        "Keep top-K ants/solutions per iteration (1..antsN). Used for ranked pheromone update. Recommended: 5-10 for better exploration.", "Algorithm",
    [](ParamsData& d, const std::string& v){ d.eliteKAnts = parseInt(v, "eliteKAnts"); } });

    s.push_back({ "favorites", "fv", "int", "20",
    "Candidate list size: limits next-node choices in construction and speeds up local search. Controls exploration vs exploitation balance. Larger => better quality (more exploration), slower; smaller => faster, may degrade quality. Recommended: 0.5N to 0.67N (e.g., 25-33 for N=50, 50-67 for N=100, 150-200 for N=300). Default: 20.", "Algorithm",
    [](ParamsData& d, const std::string& v){ d.favorites = parseInt(v, "favorites"); } });

    s.push_back({ "policyQ0","pQ0","double", "0.1",
    "Cars policy: exploitation probability in [0,1]. With prob q0 choose argmax(weight), else roulette. Higher => more greedy, faster convergence but less exploration. Recommended: 0.05-0.2.", "Algorithm",
    [](ParamsData& d, const std::string& v){ d.policyQ0 = parseDouble(v, "q0Cars"); } });

    // Adaptive q0 parameters
    s.push_back({ "adaptiveQ0", "aq0", "bool", "0",
    "Enable adaptive q0: q0 increases linearly from q0Start to q0End during iterations. 0=disabled (use fixed policyQ0), 1=enabled. When enabled, q0 starts low (exploration) and increases to higher value (exploitation) as iterations progress.", "Algorithm",
    [](ParamsData& d, const std::string& v){ d.adaptiveQ0 = parseBool01(v, "adaptiveQ0"); } });
    
    s.push_back({ "q0Start", "q0s", "double", "0.0",
    "Adaptive q0: starting value for early iterations (more exploration). Range: [0,1]. Lower values = more exploration, higher = more exploitation. Recommended: 0.0-0.2.", "Algorithm",
    [](ParamsData& d, const std::string& v){ d.q0Start = parseDouble(v, "q0Start"); } });
    
    s.push_back({ "q0End", "q0e", "double", "0.9",
    "Adaptive q0: ending value for late iterations (more exploitation). Range: [0,1]. Higher values = more exploitation, lower = more exploration. Recommended: 0.7-0.95.", "Algorithm",
    [](ParamsData& d, const std::string& v){ d.q0End = parseDouble(v, "q0End"); } });


    // === Pheromone ===
    s.push_back({ "carAlpha",  "ca", "double",   "1.5",    "Alpha (cars): pheromone weight exponent. Higher => more influence of pheromone trails. Typical: 1.0-2.0.", "Pheromone",
        [](ParamsData& d, const std::string& v){ d.alphaCars = parseDouble(v, "carAlpha"); } });

    s.push_back({ "carBeta",   "cb", "double",   "1.1",    "Beta (cars): heuristic weight exponent. Higher => more influence of distance/cost heuristic. Typical: 1.0-5.0.", "Pheromone",
        [](ParamsData& d, const std::string& v){ d.betaCars = parseDouble(v, "carBeta"); } });

    s.push_back({ "carRho",    "cr", "double",   "0.1",    "Rho (cars): evaporation rate in (0,1). Higher => faster forgetting, more exploration. Lower => slower convergence, more exploitation. Recommended: 0.05-0.2.", "Pheromone",
        [](ParamsData& d, const std::string& v){ d.rhoCars = parseDouble(v, "carRho"); } });

    s.push_back({ "carMaxMin", "cm", "double",   "0",      "MaxMin (cars): tauMax/tauMin ratio. 0 disables. If >0, sets tauMin = tauMax/maxMin. Typical: 20-100. Overridden by pBestCars if pBestCars > 0.", "Pheromone",
        [](ParamsData& d, const std::string& v){ d.maxMinCars = parseDouble(v, "carMaxMin"); } });

    s.push_back({ "nodeAlpha", "na", "double",   "1.5",    "Alpha (nodes)", "Pheromone",
        [](ParamsData& d, const std::string& v){ d.alphaNodes = parseDouble(v, "nodeAlpha"); } });

    s.push_back({ "nodeBeta",  "nb", "double",   "1.1",    "Beta (nodes)", "Pheromone",
        [](ParamsData& d, const std::string& v){ d.betaNodes = parseDouble(v, "nodeBeta"); } });

    s.push_back({ "nodeRho",   "nr", "double",   "0.1",    "Rho (nodes)", "Pheromone",
        [](ParamsData& d, const std::string& v){ d.rhoNodes = parseDouble(v, "nodeRho"); } });

    s.push_back({ "nodeMaxMin","nm", "double",   "0",      "MaxMin (nodes), 0 disables", "Pheromone",
        [](ParamsData& d, const std::string& v){ d.maxMinNodes = parseDouble(v, "nodeMaxMin"); } });

    // === Pheromone (continued) ===
    s.push_back({ "globalBestPeriod","gb","int", "0", "Global best period: 0 auto, >0 fixed, -1 off", "Pheromone",
        [](ParamsData& d, const std::string& v){ d.gbPeriod = parseInt(v, "globalBestPeriod"); } });

    s.push_back({ "smGammaCars","sgc","double", "0.5", "Smoothing gamma (cars) in [0,1]", "Pheromone",
        [](ParamsData& d, const std::string& v){ d.smoothingGammaCars = parseDouble(v, "smGammaCars"); } });

    s.push_back({ "restartTargetCars","rtc","string", "mid", "Restart target for smoothing: mid | tauMax", "Pheromone",
        [](ParamsData& d, const std::string& v){ d.restartTargetCars = parseRestartTarget(v, "restartTargetCars"); } });

    s.push_back({ "pBestCars","pbc","double", "0",
        "MMAS p_best in (0,1). If > 0 enables pBest-based tauMin (overrides maxMinCars). Use 0 to disable.", "Pheromone",
        [](ParamsData& d, const std::string& v){ d.pBestCars = parseDouble(v, "pBestCars"); } });

    s.push_back({ "tbas","tb","bool", "0",
        "Pheromone model: 0 = MMAS (default, with evaporation), 1 = TBAS (Three Bounds Ant System, no evaporation, uses three bounds: tau_LB, tau_UB, tau_CB with clipping). TBAS uses Q variable updated as Q_i = Q_{i-1}/(1-rho) each iteration.", "Pheromone",
        [](ParamsData& d, const std::string& v){ d.useTBAS = parseBool01(v, "tbas"); } });

    // === Exploration ===
    s.push_back({ "exploreCars","exc","bool", "0", "Enable exploration after stagnation (0/1)", "Exploration",
        [](ParamsData& d, const std::string& v){ d.exploreEnabledCars = parseBool01(v, "exploreCars"); } });

    s.push_back({ "exploreItersCars","exi","int", "20", "Exploration iterations after stagnation", "Exploration",
        [](ParamsData& d, const std::string& v){ d.exploreItersCars = parseInt(v, "exploreItersCars"); } });

    s.push_back({ "exploreRhoMultCars","exrm","double", "2.0", "Exploration rho multiplier (>=1)", "Exploration",
        [](ParamsData& d, const std::string& v){ d.exploreRhoMultiplierCars = parseDouble(v, "exploreRhoMultCars"); } });

    // === Advanced ===
    s.push_back({ "epsilonEtaCars","eec","double", "1e-6", "Epsilon in eta=1/(eps+cost) for cars", "Advanced",
        [](ParamsData& d, const std::string& v){ d.epsilonEtaCars = parseDouble(v, "epsilonEtaCars"); } });

    s.push_back({ "returnFactorCars","rfc","double", "0.0",
    "Cars policy only: return influence factor. 0 => no returnCost/tauReturn in construction (recommended). Deposit still uses return edges via DP.", "Advanced",
    [](ParamsData& d, const std::string& v){ d.returnFactorCars = parseDouble(v, "returnFactorCars"); } });

    s.push_back({ "sampleCars","sc","bool", "1",
    "Cars policy: if 1, ant samples (nextNode, car) per edge; if 0, builds node tour only (cars assigned later by DP).", "Advanced",
    [](ParamsData& d, const std::string& v){ d.sampleCars = parseBool01(v, "sampleCars"); } });

    s.push_back({ "enforceNoRerent","enr","bool", "0",
    "Cars policy: if 1, enforce 'no re-rent after return' already during construction (extra overhead). DP always enforces it anyway.", "Advanced",
    [](ParamsData& d, const std::string& v){ d.enforceNoRerent = parseBool01(v, "enforceNoRerent"); } });

    // === Local Search ===
    s.push_back({ "lsTopW","lsw","int", "2",
    "Local search width: number of elite solutions to improve with LS per iteration. 0 = dynamic (starts at 2, grows only during stagnation), >0 = fixed. Default: 2.", "Local Search",
    [](ParamsData& d, const std::string& v){ d.lsTopW = parseInt(v, "lsTopW"); } });

    s.push_back({ "useFullPolish", "ufp", "bool", "0",
    "Use full polish algorithms (LK, 3-opt) instead of lite versions. 0 = lite (faster, default), 1 = full (slower but may find better solutions).", "Local Search",
    [](ParamsData& d, const std::string& v){ d.useFullPolish = parseBool01(v, "useFullPolish"); } });

    // === DP ===
    s.push_back({ "dpCacheCapacity","dpc","int", "5000",
    "DP cache capacity: number of cached tour evaluations. 0 = disabled. Larger => more memory, fewer recomputations. Recommended: 5000-20000 for N~300.", "DP",
    [](ParamsData& d, const std::string& v){ d.dpCacheCapacity = parseInt(v, "dpCacheCapacity"); } });

    s.push_back({ "dpMaxSegmentLenOffset","dpso","int", "5",
    "DP maxSegmentLen offset: maxSegmentLen = ceil(N/C) + offset. Larger => more accurate but slower. Use -1 for unlimited (exact DP).", "DP",
    [](ParamsData& d, const std::string& v){ d.dpMaxSegmentLenOffset = parseInt(v, "dpMaxSegmentLenOffset"); } });

    // === Intensifier ===
    s.push_back({ "intensifier","intf","bool", "1",
    "Intensifier on/off: 1 = enabled (extra LS on global best + related pheromone/archive reset logic), "
    "0 = disabled: no such component or extra LS on global best; LS on elite ants remains. Default: 1. "
    "CSV: RunRecorder header (param;value) includes intensifierEnabled=0|1 in all outputs that write a params header (e.g. iteration_best.csv, global_best.csv, weibulls.csv, weibull_hypotheses.csv, runs_summary.csv, etc.).",
    "Intensifier",
    [](ParamsData& d, const std::string& v){ d.intensifierEnabled = parseBool01(v, "intensifier"); } });
    s.push_back({ "intensifierMaxStagnationIterations","imsi","int", "3",
    "Intensifier max stagnation iterations: number of iterations after stagnation when intensifier is still active. Intensifier applies additional LS cycle on best solution when best is improved or within this many iterations after stagnation. Default: 3.", "Intensifier",
    [](ParamsData& d, const std::string& v){ d.intensifierMaxStagnationIterations = parseInt(v, "intensifierMaxStagnationIterations"); } });
    
    s.push_back({ "intensifierResetPheromonesCount","irph","int", "1",
    "Intensifier reset pheromones count: how many deactivations trigger pheromone reset. 0 or negative = never reset. 1 = reset after 1st deactivation, 2 = after 1st and 2nd, etc. Default: 1.", "Intensifier",
    [](ParamsData& d, const std::string& v){ d.intensifierResetPheromonesCount = parseInt(v, "intensifierResetPheromonesCount"); } });
    
    s.push_back({ "intensifierPheromoneReductionGamma","iprg","double", "0.3",
    "Intensifier pheromone reduction gamma: factor for selective pheromone reduction on best tour (0.0=full reset to tauBase, 1.0=no change). Used instead of full reset when resetPheromones is enabled. Lower values = stronger reduction. Default: 0.3.", "Intensifier",
    [](ParamsData& d, const std::string& v){ d.intensifierPheromoneReductionGamma = parseDouble(v, "intensifierPheromoneReductionGamma"); } });
    
    s.push_back({ "intensifierArchiveAndResetBest","iarab","bool", "1",
    "Intensifier archive and reset best: archive global best and reset it to infinity after first deactivation of intensifier. This allows complete restart while preserving best solution in archive. Default: 1 (enabled).", "Intensifier",
    [](ParamsData& d, const std::string& v){ d.intensifierArchiveAndResetBest = parseBool01(v, "intensifierArchiveAndResetBest"); } });
    
    s.push_back({ "ablationDisableAntDepositsDuringIntensification","adai","bool", "0",
    "Ablation (pheromones vs intensifier). 0=off: normal MMAS/TBAS end-of-iteration update (evaporation + deposits from ranked ants and global-best rules). 1=on: while intensifier is active in an iteration, skip the entire pheromone update (no evaporation, no deposits, no archive-memory deposit that iteration). Default: 0.",
    "Ablation",
    [](ParamsData& d, const std::string& v){ d.ablationDisableAntDepositsDuringIntensification = parseBool01(v, "ablationDisableAntDepositsDuringIntensification"); } });
    
    s.push_back({ "ablationResetPheromonesOnIntensifierActivation","aria","bool", "0",
    "Ablation (pheromones vs intensifier). 0=off: no extra reset on intensifier activation. 1=on: each time intensifier goes inactive->active, call pheromone reset() and skip the entire pheromone update in that iteration (no deposits). Default: 0.",
    "Ablation",
    [](ParamsData& d, const std::string& v){ d.ablationResetPheromonesOnIntensifierActivation = parseBool01(v, "ablationResetPheromonesOnIntensifierActivation"); } });

    s.push_back({ "ablationIntensifierOnlyNoAnts","aion","bool", "0",
    "Ablation (no ant colony construction). 0=off: full ACO — ants build tours each iteration, ranked elite, normal pipeline. 1=on: no ant construction; iter 0 starts from one tour only (see -airit); later iters copy current global best as sole elite; pheromone updates are always skipped. Compare intensifier+LS vs full ACO. Default: 0.",
    "Ablation",
    [](ParamsData& d, const std::string& v){ d.ablationIntensifierOnlyNoAnts = parseBool01(v, "ablationIntensifierOnlyNoAnts"); } });
    s.push_back({ "ablationIntensifierOnlyRandomInitialTour","airit","bool", "1",
    "Only when -aion 1 (iter 0 seed tour). 0=identity Hamilton cycle 0,1,...,N-1. 1=random: shuffle nodes 1..N-1, depot 0 fixed. Default: 1.",
    "Ablation",
    [](ParamsData& d, const std::string& v){ d.ablationIntensifierOnlyRandomInitialTour = parseBool01(v, "ablationIntensifierOnlyRandomInitialTour"); } });

    s.push_back({ "surrogateReturnGamma","srg","double", "0.5",
    "Surrogate return term: add gamma * (est. returns) * (avg return cost) to surrogate for better correlation with full DP cost. 0 = off (old surrogate), 0.3-0.5 = recommended. Default: 0.5.", "Surrogate",
    [](ParamsData& d, const std::string& v){ d.surrogateReturnGamma = parseDouble(v, "surrogateReturnGamma"); } });

    s.push_back({ "researchLog","rlog","bool", "0",
    "Research logging: write intensifier_binomial.csv (run, activation_index, cost_start, cost_end, improved, depth) for binomial analysis and ablations; depth = count of global-improvement stagnation counter resets during that activation while intensifier is active. Default: 0 (off).", "Research",
    [](ParamsData& d, const std::string& v){ d.researchLog = parseBool01(v, "researchLog"); } });

    s.push_back({ "researchSurrogateLog","rsur","bool", "0",
    "Research surrogate vs DP cost: write surrogate_correlation.csv (run, iteration, ant_index, surrogate, dp_cost) to validate correlation (~0.8). Slows run: evaluates full DP for every ant. Suggested: instances up to 100 nodes, 20-30 it, 2-3 runs. Default: 0 (off).", "Research",
    [](ParamsData& d, const std::string& v){ d.researchSurrogateLog = parseBool01(v, "researchSurrogateLog"); } });

    // === LS Quality ===
    s.push_back({ "lsQuality", "lsq", "int", "3",
    "Local Search quality mode: 1=low (fast, less precise), 2=medium, 3=high (default, balanced), 4=extreme (slowest, most precise). Sets maxPasses for TwoOpt #1, Relocation, TwoOpt #2, and DP retries.",
    "Quality",
    [](ParamsData& d, const std::string& v){ d.lsQualityMode = parseInt(v, "lsQuality"); } });

    s.push_back({ "lsTwoOpt1Passes", "ls2o1p", "int", "-1",
    "Override TwoOpt #1 maxPasses (-1 = use lsQuality mode, >0 = override value).",
    "Quality",
    [](ParamsData& d, const std::string& v){ d.lsTwoOpt1Passes = parseInt(v, "lsTwoOpt1Passes"); } });

    s.push_back({ "lsRelocPasses", "lsrp", "int", "-1",
    "Override Relocation maxPasses (-1 = use lsQuality mode, >0 = override value).",
    "Quality",
    [](ParamsData& d, const std::string& v){ d.lsRelocPasses = parseInt(v, "lsRelocPasses"); } });

    s.push_back({ "lsTwoOpt2Passes", "ls2o2p", "int", "-1",
    "Override TwoOpt #2 maxPasses (-1 = use lsQuality mode, >0 = override value).",
    "Quality",
    [](ParamsData& d, const std::string& v){ d.lsTwoOpt2Passes = parseInt(v, "lsTwoOpt2Passes"); } });

    s.push_back({ "lsDpRetries", "lsdr", "int", "-1",
    "Override DP fallback maxRetries (-1 = use lsQuality mode, >0 = override value).",
    "Quality",
    [](ParamsData& d, const std::string& v){ d.lsDpRetries = parseInt(v, "lsDpRetries"); } });

    s.push_back({ "useReinsertion", "ri", "bool", "0",
    "Enable Reinsertion local search (single-node relocation with candidate list). LS chain becomes: 2-opt -> Relocation -> Reinsertion -> 2-opt. Default: 0 (disabled).",
    "Quality",
    [](ParamsData& d, const std::string& v){ d.useReinsertion = parseBool01(v, "useReinsertion"); } });

    s.push_back({ "reinsertionMaxPasses", "rimp", "int", "5",
    "Reinsertion: max passes when enabled. Default: 5.",
    "Quality",
    [](ParamsData& d, const std::string& v){ d.reinsertionMaxPasses = parseInt(v, "reinsertionMaxPasses"); } });

    return s;
}

static inline const OptionSpec* findByLong(const std::vector<OptionSpec>& specs, const std::string& longName)
{
    for (const auto& sp : specs)
        if (longName == sp.longName) return &sp;
    return nullptr;
}

static inline const OptionSpec* findByShort(const std::vector<OptionSpec>& specs, const std::string& shortName)
{
    for (const auto& sp : specs)
        if (sp.shortName && shortName == sp.shortName) return &sp;
    return nullptr;
}

static inline void validate(const ParamsData& d)
{
    if (d.filename.empty()) throw ParseError("filename must not be empty");

    if (d.nRuns < 1) throw ParseError("nRuns must be >= 1");
    if (d.nRuns > MAX_NRUNS) throw ParseError("nRuns must be <= " + std::to_string(MAX_NRUNS));

    if (d.iterations <= 0) throw ParseError("Number of iterations must be > 0");
    if (d.antsN <= 0) throw ParseError("antsN must be > 0");

    if (d.stagnation < 0) throw ParseError("stagnation must be >= 0");
    
    if (d.eliteKAnts < 1)
        throw ParseError("eliteKAnts must be >= 1");
    if (d.eliteKAnts > d.antsN)
        throw ParseError("eliteKAnts must be <= antsN");

    
    if (d.favorites < 0) throw ParseError("favorites must be >= 0");

    if (d.gbPeriod < -1) throw ParseError("globalBestPeriod must be >= -1");

    // Optional sanity (better for MMAS)
    if (!(d.rhoCars > 0.0 && d.rhoCars < 1.0)) throw ParseError("carRho must be in (0,1)");
    if (!(d.rhoNodes > 0.0 && d.rhoNodes < 1.0)) throw ParseError("nodeRho must be in (0,1)");

    // Smoothing gamma
    if (!(d.smoothingGammaCars >= 0.0 && d.smoothingGammaCars <= 1.0))
        throw ParseError("smGammaCars must be in [0,1]");

    // Exploration
    if (d.exploreItersCars < 0)
        throw ParseError("exploreItersCars must be >= 0");

    if (!std::isfinite(d.exploreRhoMultiplierCars) || d.exploreRhoMultiplierCars < 1.0)
        throw ParseError("exploreRhoMultCars must be finite and >= 1");

    if (!std::isfinite(d.policyQ0) || d.policyQ0 < 0.0 || d.policyQ0 > 1.0)
        throw ParseError("policyQ0 must be finite and in [0,1]");
    
    // Adaptive q0 validation
    if (d.adaptiveQ0) {
        if (!std::isfinite(d.q0Start) || d.q0Start < 0.0 || d.q0Start > 1.0)
            throw ParseError("q0Start must be finite and in [0,1]");
        if (!std::isfinite(d.q0End) || d.q0End < 0.0 || d.q0End > 1.0)
            throw ParseError("q0End must be finite and in [0,1]");
        if (d.q0Start > d.q0End)
            throw ParseError("q0Start must be <= q0End");
    }

    // epsilonEta
    if (!std::isfinite(d.epsilonEtaCars) || d.epsilonEtaCars < 0.0)
        throw ParseError("epsilonEtaCars must be finite and >= 0");

    // OPTIONAL: pBest (disabled if == 0)
    if (!std::isfinite(d.pBestCars))
        throw ParseError("pBestCars must be finite");

    if (d.pBestCars < 0.0)
        throw ParseError("pBestCars must be >= 0 (use 0 to disable)");

    if (d.pBestCars > 0.0 && !(d.pBestCars < 1.0))
        throw ParseError("pBestCars must be in (0,1) when enabled");

    // returnFactorCars
    if (!std::isfinite(d.returnFactorCars) || d.returnFactorCars < 0.0)
        throw ParseError("returnFactorCars must be finite and >= 0");

    // Local search and DP tuning
    if (d.lsTopW < 0)
        throw ParseError("lsTopW must be >= 0 (0 = auto)");

    if (d.dpCacheCapacity < 0)
        throw ParseError("dpCacheCapacity must be >= 0 (0 = disabled)");

    if (d.dpMaxSegmentLenOffset < -1)
        throw ParseError("dpMaxSegmentLenOffset must be >= -1 (-1 = unlimited)");

    if (d.reinsertionMaxPasses < 1)
        throw ParseError("reinsertionMaxPasses must be >= 1");

    if (d.intensifierMaxStagnationIterations < 0)
        throw ParseError("intensifierMaxStagnationIterations must be >= 0");

    if (d.intensifierPheromoneReductionGamma < 0.0 || d.intensifierPheromoneReductionGamma > 1.0)
        throw ParseError("intensifierPheromoneReductionGamma must be in [0.0, 1.0]");

    // Intensifier validation
    // (no additional validation needed for bool parameters)

    if (d.surrogateReturnGamma < 0.0)
        throw ParseError("surrogateReturnGamma must be >= 0.0");

}

} // namespace

Parser::Parser(ParserOptions opt) : opt_(opt) {}

ParamsData Parser::parse(int argc, const char* const* argv) const
{
    if (opt_.printHelpOnEmpty && argc <= 1)
        throw ParseError("HELP_REQUESTED");

    ParamsData data;

    // Defaults (explicit)
    data.seed = defaultSeedFromTime();
    data.filename = DEFAULT_FILENAME;

    data.nRuns = DEFAULT_NRUNS;
    data.iterations = DEFAULT_ITERATIONS;
    data.antsN = DEFAULT_ANT_NUMBER;

    data.alphaCars = DEFAULT_ALPHA;
    data.betaCars  = DEFAULT_BETA;
    data.rhoCars   = DEFAULT_RHO;

    data.alphaNodes = DEFAULT_ALPHA;
    data.betaNodes  = DEFAULT_BETA;
    data.rhoNodes   = DEFAULT_RHO;

    data.maxMinCars  = DEFAULT_MAXMIN;
    data.maxMinNodes = DEFAULT_MAXMIN;

    data.stagnation = DEFAULT_STAGNATION;

    data.eliteKAnts = DEFAULT_ELITE_K_ANTS;
    data.favorites  = DEFAULT_FAVORITES;
    data.gbPeriod = DEFAULT_GB_PERIOD;

    // Pass (if used elsewhere)
    data.alphaPass = DEFAULT_ALPHA;
    data.betaPass  = DEFAULT_BETA;
    data.rhoPass   = DEFAULT_RHO;

    // New defaults (explicit)
    data.smoothingGammaCars = DEFAULT_SMOOTHING_GAMMA;
    data.restartTargetCars = RestartTarget::Mid;

    data.exploreEnabledCars = (DEFAULT_EXPLORE_ENABLED != 0);
    data.exploreItersCars = DEFAULT_EXPLORE_ITERS;
    data.exploreRhoMultiplierCars = DEFAULT_EXPLORE_RHO_MULT;

    data.policyQ0 = DEFAULT_POLICY_Q0;
    data.adaptiveQ0 = (DEFAULT_ADAPTIVE_Q0 != 0);
    data.q0Start = DEFAULT_Q0_START;
    data.q0End = DEFAULT_Q0_END;

    data.epsilonEtaCars = DEFAULT_EPSILON_ETA;

    data.pBestCars = DEFAULT_PBEST;

    data.returnFactorCars = DEFAULT_RETURN_FACTOR_CARS;
    data.sampleCars = (DEFAULT_SAMPLE_CARS != 0);
    data.enforceNoRerent = (DEFAULT_ENFORCE_NO_RERENT != 0);

    data.lsTopW = DEFAULT_LS_TOP_W;
    data.dpCacheCapacity = DEFAULT_DP_CACHE_CAPACITY;
    data.dpMaxSegmentLenOffset = DEFAULT_DP_MAX_SEGMENT_LEN_OFFSET;
    data.intensifierEnabled = true;
    data.intensifierMaxStagnationIterations = DEFAULT_INTENSIFIER_MAX_STAGNATION_ITERATIONS;
    data.intensifierResetPheromonesCount = DEFAULT_INTENSIFIER_RESET_PHEROMONES_COUNT;
    data.intensifierPheromoneReductionGamma = DEFAULT_INTENSIFIER_PHEROMONE_REDUCTION_GAMMA;
    data.intensifierArchiveAndResetBest = (DEFAULT_INTENSIFIER_ARCHIVE_AND_RESET_BEST != 0);
    data.surrogateReturnGamma = 0.5;
    data.useTBAS = false;  // default: MMAS


    const auto specs = makeSpecs();

    auto needValue = [&](const std::string& shownKey, int i) -> std::string 
    {
        if (i + 1 >= argc)
            throw ParseError("Missing value for parameter '" + shownKey + "'");
        return std::string(argv[i + 1]);
    };

    for (int i = 1; i < argc; ++i)
    {
        std::string tok(argv[i]);

        if (tok == "--help" || tok == "-h")
            throw ParseError("HELP_REQUESTED");

        if (tok == "--") // stop parsing options; rest are positional (currently not supported)
        {
            if (opt_.strictUnknownParams && i + 1 < argc)
                throw ParseError("Unexpected positional argument: '" + std::string(argv[i + 1]) + "'");
            break;
        }

        if (!startsWith(tok, "-"))
        {
            if (opt_.strictUnknownParams)
                throw ParseError("Unexpected positional argument: '" + tok + "'");
            continue;
        }

        // Determine if long (--x) or short (-x / -ca / -st ...)
        const bool isLong = startsWith(tok, "--");
        auto [k0, v0] = splitKeyEqValue(tok);

        std::string key = stripLeadingDashes(k0);
        std::string value = v0.empty() ? needValue(tok, i) : v0;
        if (v0.empty()) ++i; // consumed next argv as value

        const OptionSpec* sp = isLong ? findByLong(specs, key) : findByShort(specs, key);

        if (!sp)
        {
            if (opt_.strictUnknownParams)
            {
                const std::string prefix = isLong ? "--" : "-";
                throw ParseError("Unknown parameter '" + prefix + key + "'. Use --help.");
            }
            continue;
        }

        sp->set(data, value);
    }

    validate(data);
    return data;
}

void Parser::printUsage(std::ostream& os, const std::string& exeName) const
{
    const auto specs = makeSpecs();

    os << "Usage:\n"
       << "  " << exeName << " [options]\n\n"
       << "Options:\n"
       << "  -h, --help                Show this help and exit\n\n";

    // Group by category
    std::map<std::string, std::vector<const OptionSpec*>> byCategory;
    std::vector<const OptionSpec*> noCategory;

    for (const auto& sp : specs)
    {
        if (sp.category && std::string(sp.category).size() > 0)
        {
            byCategory[sp.category].push_back(&sp);
        }
        else
        {
            noCategory.push_back(&sp);
        }
    }

    // Print categorized options
    const char* categories[] = {
        "Basics",
        "Algorithm",
        "Pheromone",
        "Quality",
        "Local Search",
        "DP",
        "Intensifier",
        "Ablation",
        "Exploration",
        "Surrogate",
        "Research",
        "Advanced"
    };

    for (const char* cat : categories)
    {
        auto it = byCategory.find(cat);
        if (it != byCategory.end() && !it->second.empty())
        {
            os << "\n" << cat << ":\n";
            for (const auto* sp : it->second)
            {
                std::ostringstream names;
                if (sp->shortName) names << "-" << sp->shortName << ", ";
                else names << "    ";

                names << "--" << sp->longName;

                os << "  " << std::left << std::setw(24) << names.str()
                   << " <" << std::left << std::setw(8) << sp->typeHint << ">  "
                   << sp->help;

                if (sp->defaultHint && std::string(sp->defaultHint).size() > 0)
                    os << " (default: " << sp->defaultHint << ")";

                os << "\n";
            }
        }
    }

    // Print uncategorized options
    if (!noCategory.empty())
    {
        os << "\nOther:\n";
        for (const auto* sp : noCategory)
        {
            std::ostringstream names;
            if (sp->shortName) names << "-" << sp->shortName << ", ";
            else names << "    ";

            names << "--" << sp->longName;

            os << "  " << std::left << std::setw(24) << names.str()
               << " <" << std::left << std::setw(8) << sp->typeHint << ">  "
               << sp->help;

            if (sp->defaultHint && std::string(sp->defaultHint).size() > 0)
                os << " (default: " << sp->defaultHint << ")";

            os << "\n";
        }
    }

    os << "\nNotes:\n"
       << "  * Short options here are '-x' or multi-letter like '-ca' (not bundled like '-abc').\n"
       << "  * You can pass values as '--key value' or '--key=value' (same for short forms).\n"
       << "  * CSV files under outputData/ that start with a param;value header list resolved settings using internal field names as keys (e.g. -intf/--intensifier appears as intensifierEnabled).\n";
}

} // namespace aco_cli
