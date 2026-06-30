#pragma once
#include <interfaces.hpp>
#include <memory>

namespace cars_tsplib { struct Instance; }

namespace aco 
{

/**
 * Intensifier: Komponenta za intenzifikaciju (dublju eksploataciju) najboljeg rješenja.
 * 
 * OVAJ KOMPONENTA:
 * - Primjenjuje dodatni Local Search (LS) ciklus na global best rješenje
 * - Aktivira se kada je best upravo promijenjen ili u prvih nekoliko iteracija nakon stagnacije
 * - Omogućava dublju eksploataciju bez prevelikog utroška vremena (radi samo na best rješenju)
 * 
 * KADA SE AKTIVIRA:
 * 1. Kada je best upravo promijenjen u ovoj iteraciji (improvedByIter == true)
 *    - Ovo znači da smo našli novo bolje rješenje i želimo ga dodatno poboljšati
 * 2. Ili kada smo u prvih maxStagnationIterations_ iteracija nakon stagnacije
 *    - Default: 3 iteracije nakon stagnacije (može se override-ati preko konstruktora)
 *    - Ovo omogućava dodatno poboljšanje čak i kada nema novog globalnog poboljšanja
 * 
 * KAKO SE PRIKLJUČUJE:
 * - Kreira se u main.cpp factory funkciji (lambda) prije kreiranja Colony objekta
 * - Prosljeđuje se u Colony konstruktor kao zadnji parametar (opcionalan, default nullptr)
 * - Colony automatski kreira intenzifikator ako nije proslijeđen ali postoje potrebne komponente
 *   (localSearch_ i dp_)
 * - Koristi se u Colony::run() metodi nakon što se ažurira rr.best s iterBestCost
 * 
 * PARAMETRI:
 * - localSearch_: Local Search komponenta koja se koristi za poboljšanje (npr. ChainedLocalSearch)
 * - dp_: Dynamic Programming komponenta za optimalno dodjeljivanje automobila (opcionalno)
 * - inst_: Instance problema (potrebno za pristup N, cars, itd.)
 * - maxStagnationIterations_: Broj iteracija nakon stagnacije kada se još aktivira (default: 3)
 * 
 * LOGIKA RADA:
 * 1. Provjerava je li potrebna intenzifikacija (improvedByIter ili stagnacija)
 * 2. Kopira best rješenje i primjenjuje LS ciklus
 * 3. Ako je dostupan DP, konsolidira car assignment (optimalno dodjeljivanje automobila)
 * 4. Ako je poboljšanje (cost < best.cost - improveEps), ažurira best i vraća novi cost
 * 5. Inače vraća infinity (nije bilo poboljšanja)
 */
class Intensifier final : public IIntensifier
{
public:
    /**
     * Konstruktor s DP komponentom.
     * 
     * @param localSearch Local Search komponenta za poboljšanje rješenja
     * @param dp Dynamic Programming komponenta za optimalno dodjeljivanje automobila
     * @param inst Instance problema
     * @param maxStagnationIterations Broj iteracija nakon stagnacije kada se još aktivira (default: 3)
     * @param resetPheromonesOnDeactivationCount Broj deaktivacija nakon kojih se resetiraju feromoni (0 ili negativno = nikad, 1 = nakon 1. gašenja, 2 = nakon 1. i 2., itd., default: 1)
     */
    Intensifier(std::shared_ptr<const ILocalSearch> localSearch,
                std::shared_ptr<const IFixedTourCarAssigner> dp,
                std::shared_ptr<const cars_tsplib::Instance> inst,
                int maxStagnationIterations = 3,
                int resetPheromonesOnDeactivationCount = 1);
    
    /**
     * Konstruktor bez DP komponente (samo LS).
     * 
     * @param localSearch Local Search komponenta za poboljšanje rješenja
     * @param inst Instance problema
     * @param maxStagnationIterations Broj iteracija nakon stagnacije kada se još aktivira (default: 3)
     * @param resetPheromonesOnDeactivationCount Broj deaktivacija nakon kojih se resetiraju feromoni (0 ili negativno = nikad, 1 = nakon 1. gašenja, itd., default: 1)
     */
    Intensifier(std::shared_ptr<const ILocalSearch> localSearch,
                std::shared_ptr<const cars_tsplib::Instance> inst,
                int maxStagnationIterations = 3,
                int resetPheromonesOnDeactivationCount = 1);
    
    /**
     * Primjenjuje intenzifikaciju na best rješenje.
     * 
     * @param best Global best rješenje (može biti ažurirano ako postoji poboljšanje)
     * @param improvedByIter Je li best upravo promijenjen u ovoj iteraciji
     * @param itersSinceGlobalImprovement Broj iteracija od zadnjeg globalnog poboljšanja
     * @param improveEps Epsilon za provjeru poboljšanja (numerička preciznost)
     * @return Cost poboljšanja (ili infinity ako nije bilo poboljšanja)
     */
    double intensify(EvaluatedSolution& best,
                    bool improvedByIter,
                    int itersSinceGlobalImprovement,
                    double improveEps) const override;

    /**
     * Vraća broj puta koliko je intenzifikator bio uključen.
     * Brojač se povećava svaki put kada se intenzifikator prebaci iz neaktivnog u aktivno stanje.
     * 
     * @return Broj aktivacija intenzifikatora
     */
    int getActivationCount() const override { return activationCount_; }

    /**
     * Provjerava je li intenzifikator upravo deaktiviran (bio aktivan, sada nije).
     * 
     * @return True ako je intenzifikator upravo deaktiviran, false inače
     */
    bool wasJustDeactivated() const override;
    
    /**
     * Provjerava treba li resetirati feromone.
     * Vraća true ako je intenzifikator upravo deaktiviran (bio aktivan, sada nije),
     * parametar resetPheromonesOnDeactivationCount > 0 i brojač aktivacija <= parametru
     * (npr. parametar 1 = reset samo nakon 1. gašenja, parametar 2 = nakon 1. i 2. gašenja).
     * 
     * @return True ako treba resetirati feromone, false inače
     */
    bool shouldResetPheromones() const override;
    
    /**
     * Ažurira interno stanje za sljedeći poziv.
     * Treba se pozvati nakon što se provjeri deaktivacija i resetiranje feromona.
     */
    void updateState() const override;

private:
    std::shared_ptr<const ILocalSearch> localSearch_;  // LS komponenta za poboljšanje
    std::shared_ptr<const IFixedTourCarAssigner> dp_;  // DP komponenta za car assignment (opcionalno)
    std::shared_ptr<const cars_tsplib::Instance> inst_; // Instance problema
    int maxStagnationIterations_;  // Broj iteracija nakon stagnacije kada se još aktivira (default: 3)
    int resetPheromonesOnDeactivationCount_;   // Broj deaktivacija nakon kojih se resetiraju feromoni (0 ili negativno = nikad, 1 = nakon 1. gašenja, 2 = nakon 1. i 2., itd.)
    
    // Brojač aktivacija (mutable jer se mijenja u const metodi)
    mutable int activationCount_ = 0;  // Broj puta koliko je intenzifikator bio uključen
    mutable bool wasActiveLastTime_ = false;  // Je li bio aktivan u prethodnom pozivu (za detekciju prelaska iz neaktivnog u aktivno)
    mutable bool wasActiveThisTime_ = false;  // Je li aktivan u ovom pozivu (za detekciju deaktivacije)
};

} // namespace aco
