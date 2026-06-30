# ACO for the CaRS and TSP Problems

Ant Colony Optimization (ACO) algorithm for solving the **CaRS (Car Rental Salesman)** and **TSP (Traveling Salesman Problem)** problems.

## 📋 Overview

The project implements the **MMAS (Max-Min Ant System)** and **TBAS (Three Bounds Ant System)** variants of the ACO algorithm with advanced optimizations. The model is selected via the CLI parameter `--tbas` (0 = MMAS, default; 1 = TBAS).

- **CaRS problem**: a Hamiltonian tour with cars that can be swapped along the tour (the old car returns to its rental node, the new one continues)
- **TSP problem**: classic TSP (treated as CaRS with a single car)

## 🏗️ Architecture

### Main components

1. **Colony** (`colony.cpp`) - main algorithm loop
   - Two-phase ant ranking (surrogate + DP evaluation)
   - Elite management
   - Stagnation handling (LK-lite, 3-opt-lite, kick)

2. **Pheromone Models** (`pheromone_model.cpp`):
   - **PheromoneModelMMASMove** - MMAS pheromone model
     - Move pheromones: `τ(car, i→j)`
     - Return pheromones (for CaRS): `τ(car, from→to)`
     - Ranked update (top-W solutions)
     - Archive memory
     - Trail smoothing on stagnation
     - Pheromone evaporation
   - **PheromoneModelTBASMove** - TBAS pheromone model
     - Three bounds: `τ_LB`, `τ_UB`, `τ_CB`
     - Q variable: `Q_i = Q_{i-1} / (1 - ρ)`
     - Clipping procedure (no evaporation)
     - Ranked update (top-W solutions)
     - Archive memory
     - Trail smoothing (no evaporation)

3. **FixedTourCarAssignerDP** (`FixedTourCarAssignerDP.cpp`) - DP for optimal car assignment
   - Optimal tour segmentation per car
   - View methods for local search evaluations
   - Fallback mechanism (maxSegmentLen)

4. **CachedFixedTourCarAssigner** (`CachedFixedTourCarAssigner.cpp`) - DP cache
   - Canonical key (rotation + reverse)
   - FIFO eviction
   - Thread-safe option

5. **Local Search Chain**:
   - **TwoOptLocalSearch** - 2-opt with DP view evaluations
   - **RelocationLocalSearch** - relocation (or-opt) with DP view evaluations
   - Chain: TwoOpt → Relocation → TwoOpt

6. **Polish algorithms**:
   - **LkLiteLocalSearch** - Lin-Kernighan lite variant
   - **ThreeOptLiteLocalSearch** - 3-opt lite variant
   - Auto-selection based on instance symmetry

7. **AntPolicyCandidateListRoulette** (`antPolicy.cpp`) - tour construction
   - Candidate list support
   - q0 exploitation/exploration
   - Return cost aware (for CaRS)

## 🚀 Build

```bash
# Build
make

# Or directly:
g++ -std=c++17 -O3 -o build/main source/*.cpp -I include
```

## 💻 Usage

### Basic example

```bash
./build/main --filename inputData/BrasilMG30n.car \
             --iterations 1000 \
             --antsN 100 \
             --favorites 20 \
             --stagnation 30 \
             --carRho 0.1 \
             --carMaxMin 50
```

### Main parameters

| Parameter | Short | Default | Description |
|-----------|-------|---------|------|
| `filename` | `f` | - | Input instance (.car, .tsp, .atsp) |
| `iterations` | `it` | 1000 | Number of iterations per run |
| `antsN` | `k` | 100 | Number of ants per iteration |
| `favorites` | `fv` | 20 | Candidate list size (recommended: 0.5N to 0.67N) |
| `stagnation` | `st` | 30 | Stagnation threshold for polish/kick |
| `carRho` | `cr` | 0.1 | Pheromone evaporation rate (0,1) |
| `carMaxMin` | `cm` | 0 | Max-min ratio (0 = disabled) |
| `eliteKAnts` | - | 1 | Number of elite ants for ranked update |
| `policyQ0` | `pQ0` | 0.1 | q0 parameter (exploitation probability) |
| `adaptiveQ0` | `aq0` | 0 | Enable dynamic q0 (0=disabled, 1=enabled) |
| `q0Start` | `q0s` | 0.0 | Initial q0 value (exploration) |
| `q0End` | `q0e` | 0.9 | Final q0 value (exploitation) |

### Advanced parameters

- `--pBestCars <p>` - pBest override for tauMin (0<p<1)
- `--smGammaCars <γ>` - Trail smoothing strength [0,1]
- `--restartTargetCars <mid|tauMax>` - Restart target
- `--exploreCars <0|1>` - Exploration phase after stagnation
- `--globalBestPeriod <n>` - Global-best deposit period (-1=off, 0=auto, >0=fixed)
- `--dpCacheCapacity <n>` - DP cache capacity (default: 5000)
- `--dpMaxSegmentLenOffset <n>` - Offset for max segment length (default: 5)
- `--lsTopW <n>` / `--lsw <n>` - Local search width: number of elite solutions to improve with LS per iteration (0=dynamic: starts at 2, grows only during stagnation, >0=fixed, default: 2)

Full parameter list: `./build/main --help`

## ✨ Key features

### 1. Two-phase ant ranking

- **Phase 1**: surrogate evaluation for all ants (min travel cost per edge)
- **Phase 2**: DP `evaluateCost()` only for shortlisted candidates (max 30% of ants)
- Wildcards for diversity (10% of ants, randomly chosen from the rest)
- **Optimization**: min travel cost is precomputed once at the start (10-20% faster)

**Benefit**: speeds up the algorithm for large instances (N≈300, antsN=200) without losing quality.

### 2. DP Cache

- Canonical key (rotation + reverse canonicalization)
- FIFO eviction (default: 5000 entries)
- Thread-safe option

**Benefit**: avoids repeated DP evaluations after LS/polish.

### 3. Local Search optimizations

- **DP view methods**: LS evaluates moves without a full `reassignCars` call
- **Surrogate filtering**: Relocation uses a surrogate for pre-filtering
- **Chain optimization**: TwoOpt → Relocation → TwoOpt

### 4. Stagnation handling

- **LK-lite** or **3-opt-lite** (auto-selected based on symmetry)
- **Double-bridge kick** (if polish fails)
- **Trail smoothing** in the pheromone model (after stagnation)

### 5. Dynamic parameters

The algorithm uses several dynamic mechanisms that adapt during execution:

#### 5.1. Dynamic q0

- **Adaptive q0**: linear increase over iterations
- **Early iterations**: lower q0 → more exploration
- **Late iterations**: higher q0 → more exploitation
- **Effect**: 5-10% faster convergence

**Usage**: `--adaptiveQ0 1` (or `-aq0 1`) with optional `--q0Start` and `--q0End`

#### 5.2. Dynamic LSW (Local Search Width)

- **Normally**: LSW = 2 (faster execution)
- **During stagnation**: LSW grows up to min(K, 5) (more ants improved with LS)
- **After stagnation**: LSW returns to 2

**Usage**: `--lsTopW 0` (dynamic, default) or `--lsTopW N` (fixed)

#### 5.3. Trail Smoothing (Pheromone Smoothing)

- Activates during stagnation (2× stagnation threshold)
- Reduces pheromones toward a base value ($\tau_{\text{base}}$)
- Enables exploration of new regions

**Parameters**: `--smGammaCars <γ>`, `--restartTargetCars <mid|tauMax>`

#### 5.4. Exploration Phase (Increased Evaporation)

- Activates after trail smoothing
- Increases effective evaporation ($\rho_{\text{eff}} = \rho \cdot \text{multiplier}$)
- Lasts `exploreItersCars` iterations

**Parameters**: `--exploreCars <0|1>`, `--exploreItersCars <n>`, `--exploreRhoMultiplier <m>`

### 6. Ranked update

- Top-W solutions (default: `eliteKAnts`) deposit pheromones
- Weighted deposit (better ants = more pheromone)
- Archive memory for long-term retention

### 7. Thread-Safety

- Each run has its own `Colony` instance with its own RNG
- Pheromones are read-only during construction
- DP cache supports a thread-safe option
- Safe for parallel runs

## 📁 Project structure

```
ACO_CARS/
├── include/          # Header files
├── source/           # Implementation files
├── inputData/        # Input instances (.car, .tsp, .atsp)
├── outputData/       # CSV outputs (global_best.csv, iteration_best.csv, runs_summary.csv)
├── data/             # Measurement data accompanying the paper (see data/README.md)
└── run_*.sh          # Shell scripts for running experiments
```

## 📊 Output data

The algorithm generates 3 CSV files in `outputData/`:

1. **global_best.csv** - global best per iteration
2. **iteration_best.csv** - iteration best per iteration
3. **runs_summary.csv** - summary of all runs

## 🔧 Advanced options

### Self-test (DP consistency)

```bash
ACO_SELFTEST=1 ./build/main --filename inputData/BrasilMG30n.car
```

### Multi-threaded

The algorithm automatically uses `hardware_concurrency - 1` threads for parallel runs.

## 🐛 Known limitations

- The DP supports a maximum of 20 cars (bitmask limit)
- Cache capacity is fixed (default: 5000) - adjustable in `main.cpp`

## 👤 Author

Elvis Popović — Faculty of Organization and Informatics, University of Zagreb.

## 📄 License

Source code is released under the [MIT License](LICENSE). Measurement data in `data/` is released under [CC BY 4.0](data/README.md#license).
