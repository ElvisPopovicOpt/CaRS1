# ACO za CaRS i TSP Probleme

Ant Colony Optimization (ACO) algoritam za rješavanje **CaRS (Car Rental Salesman)** i **TSP (Traveling Salesman Problem)** problema.

## 📋 Opis

Projekt implementira **MMAS (Max-Min Ant System)** i **TBAS (Three Bounds Ant System)** varijante ACO algoritma s naprednim optimizacijama. Model se bira pomoću CLI parametra `--tbas` (0 = MMAS, default; 1 = TBAS).

- **CaRS problem**: Hamiltonova tura s automobilima koji se mogu mijenjati na turi (stari se vraća u čvor iznajmljivanja, novi ide dalje)
- **TSP problem**: Klasični TSP (tretira se kao CaRS s jednim automobilom)

## 🏗️ Arhitektura

### Glavne komponente

1. **Colony** (`colony.cpp`) - glavna petlja algoritma
   - Dvofazno rangiranje mrava (surrogate + DP evaluacija)
   - Elite management
   - Stagnation handling (LK-lite, 3-opt-lite, kick)

2. **Feromonski Modeli** (`pheromone_model.cpp`):
   - **PheromoneModelMMASMove** - MMAS feromon model
     - Move feromoni: `τ(car, i→j)`
     - Return feromoni (za CaRS): `τ(car, from→to)`
     - Ranked update (top-W rješenja)
     - Archive memory
     - Trail smoothing na stagnaciju
     - Evaporacija feromona
   - **PheromoneModelTBASMove** - TBAS feromon model
     - Tri granice: `τ_LB`, `τ_UB`, `τ_CB`
     - Q varijabla: `Q_i = Q_{i-1} / (1 - ρ)`
     - Clipping procedura (bez evaporacije)
     - Ranked update (top-W rješenja)
     - Archive memory
     - Trail smoothing (bez evaporacije)

3. **FixedTourCarAssignerDP** (`FixedTourCarAssignerDP.cpp`) - DP za optimalan car assignment
   - Optimalna segmentacija ture po automobilima
   - View metode za local search evaluacije
   - Fallback mehanizam (maxSegmentLen)

4. **CachedFixedTourCarAssigner** (`CachedFixedTourCarAssigner.cpp`) - DP cache
   - Canonical key (rotation + reverse)
   - FIFO eviction
   - Thread-safe opcija

5. **Local Search Chain**:
   - **TwoOptLocalSearch** - 2-opt s DP view evaluacijama
   - **RelocationLocalSearch** - relocation (or-opt) s DP view evaluacijama
   - Chain: TwoOpt → Relocation → TwoOpt

6. **Polish algoritmi**:
   - **LkLiteLocalSearch** - Lin-Kernighan lite varijanta
   - **ThreeOptLiteLocalSearch** - 3-opt lite varijanta
   - Auto-odabir ovisno o simetriji instance

7. **AntPolicyCandidateListRoulette** (`antPolicy.cpp`) - konstrukcija tura
   - Candidate list support
   - q0 exploitation/exploration
   - Return cost aware (za CaRS)

## 🚀 Kompilacija

```bash
# Kompajliranje
make

# Ili direktno:
g++ -std=c++17 -O3 -o build/main source/*.cpp -I include
```

## 💻 Korištenje

### Osnovni primjer

```bash
./build/main --filename inputData/BrasilMG30n.car \
             --iterations 1000 \
             --antsN 100 \
             --favorites 20 \
             --stagnation 30 \
             --carRho 0.1 \
             --carMaxMin 50
```

### Glavni parametri

| Parametar | Short | Default | Opis |
|-----------|-------|---------|------|
| `filename` | `f` | - | Ulazna instanca (.car, .tsp, .atsp) |
| `iterations` | `it` | 1000 | Broj iteracija po runu |
| `antsN` | `k` | 100 | Broj mrava po iteraciji |
| `favorites` | `fv` | 20 | Candidate list veličina (preporučeno: 0.5N do 0.67N) |
| `stagnation` | `st` | 30 | Prag stagnacije za polish/kick |
| `carRho` | `cr` | 0.1 | Evaporacija feromona (0,1) |
| `carMaxMin` | `cm` | 0 | Max-min omjer (0 = disabled) |
| `eliteKAnts` | - | 1 | Broj elite mrava za ranked update |
| `policyQ0` | `pQ0` | 0.1 | q0 parametar (exploitation vjerojatnost) |
| `adaptiveQ0` | `aq0` | 0 | Omogući dinamički q0 (0=disabled, 1=enabled) |
| `q0Start` | `q0s` | 0.0 | Početna vrijednost q0 (eksploracija) |
| `q0End` | `q0e` | 0.9 | Konačna vrijednost q0 (eksploatacija) |

### Napredni parametri

- `--pBestCars <p>` - pBest override za tauMin (0<p<1)
- `--smGammaCars <γ>` - Trail smoothing jačina [0,1]
- `--restartTargetCars <mid|tauMax>` - Restart cilj
- `--exploreCars <0|1>` - Exploration faza nakon stagnacije
- `--globalBestPeriod <n>` - Global-best deposit period (-1=off, 0=auto, >0=fixed)
- `--dpCacheCapacity <n>` - DP cache kapacitet (default: 5000)
- `--dpMaxSegmentLenOffset <n>` - Offset za max segment length (default: 5)
- `--lsTopW <n>` / `--lsw <n>` - Local search width: number of elite solutions to improve with LS per iteration (0=dynamic: starts at 2, grows only during stagnation, >0=fixed, default: 2)

Puni popis parametara: `./build/main --help`

## ✨ Ključne značajke

### 1. Dvofazno rangiranje mrava

- **Faza 1**: Surrogate evaluacija za sve mrave (min travel cost po bridu)
- **Faza 2**: DP `evaluateCost()` samo za shortlist kandidate (max 30% mrava)
- Wildcards za diverzitet (10% mrava, slučajno iz ostatka)
- **Optimizacija**: Min travel cost se precompute-ira jednom na početku (10-20% brže)

**Prednost**: Ubrzava algoritam za velike instance (N≈300, antsN=200) bez gubitka kvalitete.

### 2. DP Cache

- Canonical key (rotation + reverse canonicalizacija)
- FIFO eviction (default: 5000 zapisa)
- Thread-safe opcija

**Prednost**: Izbjegava ponovne DP evaluacije nakon LS/polisha.

### 3. Local Search optimizacije

- **DP view metode**: LS evaluira poteze bez punog `reassignCars` poziva
- **Surrogate filtering**: Relocation koristi surrogate za pre-filtering
- **Chain optimizacija**: TwoOpt → Relocation → TwoOpt

### 4. Stagnation handling

- **LK-lite** ili **3-opt-lite** (auto-odabir ovisno o simetriji)
- **Double-bridge kick** (ako polish ne uspije)
- **Trail smoothing** u feromon modelu (nakon stagnacije)

### 5. Dinamičke Veličine

Algoritam koristi nekoliko dinamičkih mehanizama koji se prilagođavaju tokom izvršavanja:

#### 5.1. Dinamički q0

- **Adaptivni q0**: Linearno povećanje tijekom iteracija
- **Rane iteracije**: Niži q0 → više eksploracije
- **Kasne iteracije**: Viši q0 → više eksploatacije
- **Učinak**: 5-10% brža konvergencija

**Korištenje**: `--adaptiveQ0 1` (ili `-aq0 1`) s opcionalnim `--q0Start` i `--q0End`

#### 5.2. Dinamički LSW (Local Search Width)

- **Normalno**: LSW = 2 (brže izvršavanje)
- **Tokom stagnacije**: LSW raste do min(K, 5) (više mrava se poboljšava LS-om)
- **Nakon stagnacije**: LSW se vraća na 2

**Korištenje**: `--lsTopW 0` (dinamički, default) ili `--lsTopW N` (fiksni)

#### 5.3. Trail Smoothing (Feromon Smoothing)

- Aktivira se tokom stagnacije (2× stagnation prag)
- Smanjuje feromone prema baznoj vrijednosti ($\tau_{\text{base}}$)
- Omogućava eksploraciju novih područja

**Parametri**: `--smGammaCars <γ>`, `--restartTargetCars <mid|tauMax>`

#### 5.4. Exploration Faza (Pojačana Evaporacija)

- Aktivira se nakon trail smoothing-a
- Povećava efektivnu evaporaciju ($\rho_{\text{eff}} = \rho \cdot \text{multiplier}$)
- Traje `exploreItersCars` iteracija

**Parametri**: `--exploreCars <0|1>`, `--exploreItersCars <n>`, `--exploreRhoMultiplier <m>`

Detaljnije: vidi `docs/DOKUMENTACIJA.md` - poglavlje "Dinamičke Veličine"

### 6. Ranked update

- Top-W rješenja (default: `eliteKAnts`) deponiraju feromone
- Weighted deposit (bolji mravi = više feromona)
- Archive memory za dugoročno pamćenje

### 7. Thread-Safety

- Svaki run ima svoju `Colony` instancu s vlastitim RNG-om
- Feromoni su read-only tijekom konstrukcije
- DP cache podržava thread-safe opciju
- Sigurno za paralelne runove

## 📁 Struktura projekta

```
ACO_2026_AI/
├── include/          # Header fajlovi
├── source/           # Implementacije
├── inputData/        # Ulazne instance (.car, .tsp, .atsp)
├── outputData/       # CSV izlazi (global_best.csv, iteration_best.csv, runs_summary.csv)
├── docs/             # Detaljna dokumentacija
│   └── ACO_CaRS_TSP.md  # MMAS dokumentacija
├── analiza/          # Excel analize
└── run_*.sh          # Shell skripte za pokretanje
```

## 📊 Izlazni podaci

Algoritam generira 3 CSV fajla u `outputData/`:

1. **global_best.csv** - Global best po iteraciji
2. **iteration_best.csv** - Iteration best po iteraciji
3. **runs_summary.csv** - Sažetak svih runova

## 🔧 Napredne opcije

### Self-test (DP konzistentnost)

```bash
ACO_SELFTEST=1 ./build/main --filename inputData/BrasilMG30n.car
```

### Multi-threaded

Algoritam automatski koristi `hardware_concurrency - 1` dretvi za paralelne runove.

## 📚 Dokumentacija

- **Detaljna MMAS dokumentacija**: `docs/ACO_CaRS_TSP.md`
- **TBAS dokumentacija**: `docs/DOKUMENTACIJA.tex` i `docs/DOKUMENTACIJA_EN.tex` (sekcija "Feromonski Modeli")
- **Analiza jezgre**: `ANALIZA_JEZGRE.md`
- **Analiza dvofaznog rangiranja**: `ANALIZA_DVOFAZNOG_RANGIRANJA.md`

## 🐛 Poznati problemi / Ograničenja

- DP podržava maksimalno 20 automobila (bitmask ograničenje)
- Cache capacity je fiksan (default: 5000) - može se podesiti u `main.cpp`

## 📝 Changelog

### Najnovije promjene

- ✅ **Dinamički q0**: Linearno povećanje q0 tijekom iteracija (5-10% brža konvergencija)
- ✅ **Precompute MinTravel**: Cache za min travel cost (10-20% brže surrogate evaluacija)
- ✅ **Thread-safety**: DP cache thread-safe opcija za paralelne runove
- ✅ **Dvofazno rangiranje**: Surrogate evaluacija (min travel cost) umjesto `costModel_->evaluate()`
- ✅ **ShortlistM**: Smanjeno na max 30% mrava (umjesto 60%)
- ✅ **DP pozivi**: Dodani `reassignCars` nakon polisha (LK-lite, 3-opt-lite)
- ✅ **DP pozivi**: Osigurana konzistentnost nakon LS lanca
- ✅ **TwoOptLocalSearch**: Koristi `evaluateCostViewScratch` umjesto `evaluateCostView`

## 👤 Autor

Projekt za ACO 2026 - Ant Colony Optimization za CaRS i TSP probleme.

## 📄 Licenca

[Ovdje dodajte licencu ako je potrebno]
