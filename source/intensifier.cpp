#include <intensifier.hpp>
#include <parser.hpp>
#include <limits>
#include <algorithm>

namespace aco 
{

// Konstruktor s DP komponentom
Intensifier::Intensifier(std::shared_ptr<const ILocalSearch> localSearch,
                        std::shared_ptr<const IFixedTourCarAssigner> dp,
                        std::shared_ptr<const cars_tsplib::Instance> inst,
                        int maxStagnationIterations,
                        int resetPheromonesOnDeactivationCount)
    : localSearch_(std::move(localSearch))
    , dp_(std::move(dp))
    , inst_(std::move(inst))
    , maxStagnationIterations_(maxStagnationIterations)
    , resetPheromonesOnDeactivationCount_(resetPheromonesOnDeactivationCount)
{
}

// Konstruktor bez DP komponente (samo LS)
Intensifier::Intensifier(std::shared_ptr<const ILocalSearch> localSearch,
                        std::shared_ptr<const cars_tsplib::Instance> inst,
                        int maxStagnationIterations,
                        int resetPheromonesOnDeactivationCount)
    : localSearch_(std::move(localSearch))
    , dp_(nullptr)
    , inst_(std::move(inst))
    , maxStagnationIterations_(maxStagnationIterations)
    , resetPheromonesOnDeactivationCount_(resetPheromonesOnDeactivationCount)
{
}

/**
 * Implementacija intenzifikacije.
 * 
 * LOGIKA:
 * 1. Provjerava je li intenzifikacija potrebna:
 *    - Mora postojati localSearch_ komponenta
 *    - best.cost mora biti finite
 *    - Treba zadovoljiti uvjete aktivacije (improvedByIter ili stagnacija)
 * 
 * 2. Uvjeti aktivacije:
 *    - improvedByIter == true: best je upravo promijenjen u ovoj iteraciji
 *    - ILI itersSinceGlobalImprovement > 0 && <= maxStagnationIterations_:
 *      smo u prvih maxStagnationIterations_ iteracija nakon stagnacije
 * 
 * 3. Ako se aktivira:
 *    - Kopira best rješenje (ne mijenja original dok ne provjeri poboljšanje)
 *    - Primjenjuje LS ciklus na kopiju (localSearch_->improve())
 *    - Ako je dostupan DP, konsolidira car assignment (optimalno dodjeljivanje automobila)
 *    - Ako je poboljšanje (cost < best.cost - improveEps), ažurira best i vraća novi cost
 *    - Inače vraća infinity (nije bilo poboljšanja)
 * 
 * 4. Greške se ignoriraju (ne blokiraju algoritam) - intenzifikacija je opcionalna optimizacija
 */
double Intensifier::intensify(EvaluatedSolution& best,
                              bool improvedByIter,
                              int itersSinceGlobalImprovement,
                              double improveEps) const
{
    // Provjeri je li intenzifikacija potrebna (osnovni uvjeti)
    if (!localSearch_ || !std::isfinite(best.cost))
    {
        return std::numeric_limits<double>::infinity();
    }
    
    // Provjeri uvjete aktivacije:
    // 1. Best je upravo promijenjen (improvedByIter == true)
    // 2. ILI smo u prvih maxStagnationIterations_ iteracija nakon stagnacije
    //    (mala stagnacija - možda još možemo poboljšati)
    const bool shouldIntensify = improvedByIter || 
                                (itersSinceGlobalImprovement > 0 && 
                                 itersSinceGlobalImprovement <= maxStagnationIterations_);
    
    // BROJAČ AKTIVACIJA:
    // Povećaj brojač ako se prebacujemo iz neaktivnog u aktivno stanje
    // (prvi puta kada se uključi, ili kada je bio isključen pa se ponovno uključi)
    if (shouldIntensify && !wasActiveLastTime_)
    {
        ++activationCount_;
    }
    
    // Zapamti trenutno stanje za provjeru deaktivacije
    // VAŽNO: wasActiveThisTime_ se postavlja prije provjere deaktivacije, ali wasActiveLastTime_ se ažurira
    // nakon što se provjeri deaktivacija (u wasJustDeactivated() ili shouldResetPheromones())
    wasActiveThisTime_ = shouldIntensify;
    
    if (!shouldIntensify)
    {
        // Ako nije aktivan, ažuriraj wasActiveLastTime_ za sljedeći poziv
        // (ali samo ako nije već ažuriran u wasJustDeactivated() ili shouldResetPheromones())
        // Zapravo, ne trebamo ga ažurirati ovdje jer će se ažurirati nakon provjere deaktivacije
        return std::numeric_limits<double>::infinity();
    }
    
    // Napravi kopiju best rješenja (ne mijenjamo original dok ne provjerimo poboljšanje)
    Solution bestCopy = best.sol;
    double bestCostCopy = best.cost;
    
    try 
    {
        // Primijeni dodatni LS ciklus na best kopiju
        // LS može mijenjati tour (node[]) i cost, ali ne mijenja car[] direktno
        localSearch_->improve(bestCopy, bestCostCopy);
        
        // DP konsolidacija (ako je dostupan)
        // DP optimalno dodjeljuje automobile za zadani tour (node[])
        // Ovo je važno jer LS može promijeniti tour, pa car assignment možda više nije optimalan
        if (dp_ && std::isfinite(bestCostCopy))
        {
            const int N = inst_->n();
            // Osiguraj da car[] ima ispravnu veličinu
            if ((int)bestCopy.car.size() != N)
                bestCopy.car.assign((size_t)N, 0);
            
            // DP optimalno dodjeljuje automobile i vraća optimalan cost
            // Uvijek koristimo DP cost jer je on optimalan za zadani tour
            const double newCost = dp_->reassignCars(bestCopy.node, bestCopy.car);
            if (std::isfinite(newCost))
                bestCostCopy = newCost;
        }
        
        // Ažuriraj best ako je poboljšanje (numerička preciznost: improveEps)
        // VAŽNO: best se ažurira direktno ovdje (pass-by-reference), Colony samo provjerava je li se promijenio
        if (std::isfinite(bestCostCopy) && bestCostCopy < best.cost - improveEps)
        {
            best.sol = std::move(bestCopy);
            best.cost = bestCostCopy;
            return bestCostCopy;  // Vrati novi cost (Colony koristi ovo za tracking)
        }
    } 
    catch (const std::exception& ex) 
    {
        // Ignoriraj greške u intensification fazi (ne blokiraj algoritam)
        // Intenzifikacija je opcionalna optimizacija - ako ne uspije, algoritam nastavlja normalno
        // Možemo dodati logging ako je potrebno za debugging
    }
    
    // Nije bilo poboljšanja - vrati infinity
    return std::numeric_limits<double>::infinity();
}

bool Intensifier::wasJustDeactivated() const
{
    // Provjeri je li intenzifikator upravo deaktiviran (bio aktivan u prethodnom pozivu, sada nije aktivan)
    return wasActiveLastTime_ && !wasActiveThisTime_;
}

bool Intensifier::shouldResetPheromones() const
{
    // Resetiraj feromone ako:
    // 1. Parametar > 0 (0 ili negativno = nikad resetirati)
    // 2. Intenzifikator je upravo deaktiviran (bio aktivan u prethodnom pozivu, sada nije aktivan)
    // 3. Brojač aktivacija <= parametru (npr. parametar 1 = reset nakon 1. gašenja, parametar 2 = nakon 1. i 2. gašenja)
    return resetPheromonesOnDeactivationCount_ > 0 &&
           wasJustDeactivated() &&
           activationCount_ <= resetPheromonesOnDeactivationCount_;
}

void Intensifier::updateState() const
{
    // Ažuriraj wasActiveLastTime_ za sljedeći poziv (nakon što se provjeri deaktivacija i resetiranje feromona)
    wasActiveLastTime_ = wasActiveThisTime_;
}

} // namespace aco
