#pragma once

#include <stdexcept>
#include <string>
#include <istream>
#include <optional>
#include <vector>

namespace cars_tsplib 
{

enum class ProblemType { TSP, ATSP, CaRS, CaRSP };

struct PassengerRequest 
{
    int origin = 0;       // 0..N-1
    int destination = 0;  // 0..N-1
    double budget = 0.0;
};

struct DenseMatrix 
{
    int n = 0;
    std::vector<double> a; // row-major, size n*n

    DenseMatrix() = default;
    explicit DenseMatrix(int n_) : n(n_), a(static_cast<std::size_t>(n_) * static_cast<std::size_t>(n_), 0.0) {}

    double& at(int i, int j) 
    {
        return a[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) + static_cast<std::size_t>(j)];
    }
    const double& at(int i, int j) const 
    {
        return a[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) + static_cast<std::size_t>(j)];
    }
};

struct Instance 
{
    // Meta
    std::string name;
    ProblemType type = ProblemType::TSP;
    std::string comment;
    

    // Core
    int dimension = 0;     // N
    int carsNumber = 1;    // TSP: 1, CaRS/CaRSP: CARS_NUMBER

    // Always explicit dense matrices
    std::vector<DenseMatrix> travelCostPerCar; // size = carsNumber
    std::vector<DenseMatrix> returnCostPerCar; // TSP: empty, CaRS/CaRSP: size = carsNumber

    // CaRSP extras
    int passengersNumber = 1; // interface-unified: TSP/CaRS => 1, CaRSP => PASSENGERS_NUMBER
    std::vector<int> passengersLimitPerCar;    // CaRSP: size = carsNumber
    std::vector<PassengerRequest> passengers;  // CaRSP: size = PASSENGERS_NUMBER, else empty

    // Optional metadata (set only if tags exist)
    std::optional<bool> edgeWeightIsSymmetric;  // EDGE_WEIGHT : SYMMETRIC
    std::optional<bool> returnRateIsAsymmetric; // RETURN_RATE : ASYMMETRIC

    // ---- Accessors for algorithms ----
    int n() const noexcept { return dimension; }
    int cars() const noexcept { return carsNumber; }
    bool hasReturnCosts() const noexcept { return !returnCostPerCar.empty(); }
    bool hasExplicitPassengers() const noexcept { return !passengers.empty(); }
    bool isCaRSP() const noexcept { return type == ProblemType::CaRSP; }

    const DenseMatrix& travel(int car) const { return travelCostPerCar.at(static_cast<std::size_t>(car)); }
    const DenseMatrix& ret(int car) const { return returnCostPerCar.at(static_cast<std::size_t>(car)); }

    double travelCost(int car, int i, int j) const { return travel(car).at(i, j); }
    // TSP has no return cost matrices even though it formally has 1 car
    double returnCost(int car, int i, int j) const { return hasReturnCosts()?ret(car).at(i, j):0.0; }

    int carCapacity(int car) const { return passengersLimitPerCar.at(static_cast<std::size_t>(car)); } // CaRSP only
};

struct ParseError : std::runtime_error 
{
    using std::runtime_error::runtime_error;
};

struct ParserOptions 
{
    // TSPLIB
    bool tsplibComputeEuc2d = true;        // EUC_2D -> explicit matrix (TSPLIB rounding)
    bool tsplibSupportTriangular = true;   // UPPER_ROW, LOWER_DIAG_ROW, ... -> FULL_MATRIX

    // Validations
    bool validateSymmetryWhenClaimed = true; // only if tag says SYMMETRIC
};

class Parser 
{
public:
    explicit Parser(ParserOptions opt = {});

    // spec can be "dir/file.tsp" or "dir?file.tsp" (we normalize '?' -> '/')
    Instance parseFileSpec(const std::string& spec) const;

    // for tests / already-open streams
    Instance parseStream(std::istream& in) const;

private:
    ParserOptions opt_;
};

void printMatrix(std::ostream& os, const DenseMatrix& m, const std::string& title = {});
void printInstance(std::ostream& os, const Instance& inst, bool printMatrices = true);

// By default, operator<< prints only a summary (matrices off).
std::ostream& operator<<(std::ostream& os, const Instance& inst);

} // namespace cars_tsplib
