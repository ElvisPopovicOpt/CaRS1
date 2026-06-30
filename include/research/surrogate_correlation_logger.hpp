#pragma once

/**
 * Istraživački modul: logiranje surogat vs. puna DP cijena za analizu korelacije.
 *
 * Tvrdnja: korelacija surogata (min travel cost po bridu) s pravom cijenom (CarAssignerDP)
 * je oko 0,8. Ovaj logger snima (surrogate, dp_cost) po konstruiranom rješenju (po mravu
 * po iteraciji) kako bi se mogla izračunati Pearson korelacija i validirati tvrdnja.
 *
 * Namjena: dokaz/validacija korelacije surogata s punim DP costom; instance do ~100 čvorova,
 * manji broj iteracija (npr. 20–30), 2–3 runa. Brzina nije prioritet — za svakog mrava
 * poziva se DP evaluateCost.
 */

#include <mutex>
#include <string>
#include <vector>

namespace aco
{

struct SurrogateDpRecord
{
    int runIndex = -1;
    int iteration = -1;
    int antIndex = -1;   // redni broj mrava u toj iteraciji (0..antsN-1)
    double surrogate = 0.0;  // suma min travel cost po bridu (surogat za LS)
    double dpCost = 0.0;     // puna cijena iz CarAssignerDP (evaluateCost)
};

class SurrogateCorrelationLogger
{
public:
    SurrogateCorrelationLogger() = default;

    /** Zapiši jedan par (surogat, dp_cost) za danog mrava u danoj iteraciji. */
    void log(int runIdx, int iter, int antIdx, double surrogate, double dpCost);

    /** Piše CSV: run;iteration;ant_index;surrogate;dp_cost. */
    void writeCsv(const std::string& filename, char sep = ';') const;

    /** Broj zapisa (parova). */
    size_t numRecords() const;

    /** Thread-safe kopija zapisa (za kasniju analizu korelacije). */
    std::vector<SurrogateDpRecord> getRecords() const;

private:
    mutable std::mutex mtx_;
    std::vector<SurrogateDpRecord> records_;
};

} // namespace aco
