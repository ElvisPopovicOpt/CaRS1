# Experimental Data — CaRS1 Hybrid ACO Algorithm

This folder contains all run-level measurement data used in the paper:

> **A Hybrid ACO Algorithm for the Car Rental Salesman Problem**
> Elvis Popović — Springer Nature Soft Computing (submitted 2026)

## Folder Structure

### `measurements/Intensifier_phero_ablation/`
Ablation study results (Table: tab:ablation).
5 instances × 6 algorithmic variants (Normal, aria=1, adai=1, adai+aria, No ACO, Weak pheromone).
Configuration: nRuns=51, base_seed=1, useReinsertion=0, alpha=1.5.

### `measurements/Varijante mjerenje 2/`
Convergence analysis results (Table: tab:convergence, Figure: convergence_gb).
5 instances × 4 variants (Full method, Without intensifier, Without forced restart, Without pheromone guidance).
Configuration: nRuns=51, base_seed=1, useReinsertion=0, alpha=1.5.

### `measurements/Varijante mjerenje 3/`
Intensifier activation analysis (Tables: tab:intensifier_contribution, tab:int_activation).
5 instances × 2 variants (Full model, Without pheromone guidance).
Configuration: nRuns=51, base_seed=1, useReinsertion=1, alpha=1.5, -rlog 1.

### `measurements/Usporedba random/`
Attribution analysis — ACO (Normal) vs Random construction (Table: tab:attribution, Figure: attribution_comparison).
13 instances, each with Normal hi_phero/ and Random/ subfolders.
Configuration: nRuns=51, base_seed=1, useReinsertion=0.

## File Types per Variant

Each variant subfolder typically contains:
- `global_best.csv` — global-best trajectory over iterations, all 51 runs
- `iteration_best.csv` — iteration-best trajectory, all 51 runs
- `runs_summary.csv` — per-run summary statistics
- `solutions.csv` — best solution found per run (node sequence)
- `dynamic_params.csv` — algorithm parameter log
- `weibulls.csv` — Weibull shape parameter per iteration
- `intensifier_binomial.csv` — per-activation intensifier log (Measurement 3 only)

XLSX files (per instance) are compiled summaries used for direct table verification.

## Instances

Instances originate from real Brazilian logistics data and TSPLIB:
BrasilPR25n (n=25), BrasilMG30n (30), BrasilCO40n (40), BrasilNO45n (45),
BrasilNE50n (50), rd100nB (100), Londrina100n (100), Osasco100n (100),
kroB150n (150), d198n (198), Aracaju200n (200), Teresina200n (200), Curitiba300n (300).

## Note on Folder Name Typo
The folder `Usporedba random/BrasilPR25m/` has a typo (m instead of n).
The instance is correctly named **BrasilPR25n** throughout the paper and XLSX files.

## Date
Data collected: 2025–2026. Uploaded to repository: June 2026.
