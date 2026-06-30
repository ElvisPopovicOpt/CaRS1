#pragma once

/**
 * Istraživački modul: logiranje intenzifikatora za binomnu analizu i statistička mjerenja.
 *
 * Snima po aktivaciji intenzifikatora: cost na početku aktivacije, cost na deaktivaciji,
 * je li došlo do poboljšanja (improved), te depth (broj reset-a stagnacijskog brojača zbog
 * globalnog poboljšanja tijekom te aktivacije, isto što i stupanj intensifier_events
 * stagnation_reset_global_improvement dok je intenzifikator aktivan u toj iteraciji).
 * Iz toga se može izvesti (n, X) po runu za
 * model Bin(n, p) — broj pokusa n = broj aktivacija, X = broj uspjeha (poboljšanja).
 *
 * Namjena: ablacije, statistička mjerenja, validacija binomne raspodjele rješenja
 * intenzifikatora. CSV se piše u outputData/ (ili drugi dir) za kasniju analizu.
 */

#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace aco
{

struct IntensifierActivationRecord
{
    int runIndex = -1;
    int activationIndex = 0;   // 1, 2, 3, ... (redni broj aktivacije u tom runu)
    double costAtStart = 0.0;  // cost besta na početku ove aktivacije
    double costAtEnd = 0.0;    // cost besta na deaktivaciji (kraj ove aktivacije)
    bool improved = false;     // je li u ovoj aktivaciji došlo do poboljšanja besta
    int depth = 0;             // broj stagnation_reset_global_improvement događaja tijekom aktivacije
};

class IntensifierBinomialLogger
{
public:
    IntensifierBinomialLogger() = default;

    /** Poziva se kada intenzifikator započne novu aktivaciju (opcionalno, za potpuniji trag). */
    void onActivationStart(int runIdx, int activationIndex, double costAtStart);

    /** Poništi brojač depth za run pri početku nove aktivacije (colony: kad activationCount poraste). */
    void resetActivationDepth(int runIdx);

    /**
     * Jedan događaj kao intensifier_events stagnation_reset_global_improvement (globalno poboljšanje)
     * u iteraciji u kojoj je intenzifikator aktivan; povećava depth za tekuću aktivaciju tog runa.
     */
    void noteGlobalStagnationResetDuringIntensifier(int runIdx);

    /**
     * Poziva se pri svakoj deaktivaciji intenzifikatora.
     * costAtStart = cost besta na početku te aktivacije, costAtEnd = na kraju, improved = je li best poboljšan.
     */
    void onDeactivation(int runIdx, int activationIndex,
                        double costAtStart, double costAtEnd, bool improved);

    /** Piše CSV: run;activation_index;cost_start;cost_end;improved;depth (samo redovi deaktivacije). */
    void writeCsv(const std::string& filename, char sep = ';') const;

    /** Vraća broj runova za koje ima barem jedan zapis (korisno za provjeru). */
    size_t numRunsWithData() const;

    /** Vraća sve zapise (thread-safe kopija). */
    std::vector<IntensifierActivationRecord> getRecords() const;

private:
    mutable std::mutex mtx_;
    std::vector<IntensifierActivationRecord> records_;
    std::unordered_map<int, int> activationDepthByRun_;
};

} // namespace aco
