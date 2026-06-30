# Research moduli (istraživačko logiranje)

## 1. Intensifier binomial (`intensifier_binomial_logger`)

- **Parametar:** `--researchLog 1` (ili `-rlog 1`)
- **Izlaz:** `outputData/intensifier_binomial.csv` — stupci: `run`, `activation_index`, `cost_start`, `cost_end`, `improved`
- **Namjena:** binomna analiza broja poboljšanja po aktivacijama intenzifikatora

## 2. Surrogate vs DP cost (`surrogate_correlation_logger`)

- **Parametar:** `--researchSurrogateLog 1` (ili `-rsur 1`)
- **Izlaz:** `outputData/surrogate_correlation.csv` — stupci: `run`, `iteration`, `ant_index`, `surrogate`, `dp_cost`
- **Namjena:** validacija tvrdnje da je korelacija surogata (min travel cost po bridu) s punom DP cijenom oko 0,8. Za svakog mrava u svakoj iteraciji snima se surogat i puna cijena (evaluateCost), pa se iz CSV-a može izračunati Pearson korelacija.

### Preporučeni parametri za surrogate analizu

- **Instance:** do ~100 čvorova (npr. kroB150n je 150 — može, ali sporije; Bra50ns, BrasilCO40n, d198n ako je manji subset, itd.)
- **Iteracije:** **20–30** (dovoljno parova za korelaciju, ne predugo)
- **Runovi:** **2–3** (npr. `-n 2` ili `-n 3`)
- **Mravi:** npr. 50–100 (više mrava = više redova u CSV-u po iteraciji)

Primjer:

```bash
./build/main --filename ./inputData/Bra50ns.car -n 2 -it 25 -k 80 --researchSurrogateLog 1
```

Brzina nije prioritet: s uključenim `researchSurrogateLog` za svakog mrava poziva se puni DP (`evaluateCost`), pa je run sporiji.
