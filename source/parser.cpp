#include <parser.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <utility>

namespace cars_tsplib 
{

namespace 
{

// ---------- string helpers ----------
static std::string trim(std::string s) 
{
    auto notSpace = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());
    return s;
}

static std::string upper(std::string s) 
{
    for (auto& ch : s) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    return s;
}

static std::vector<std::string> splitWS(const std::string& s) 
{
    std::istringstream iss(s);
    std::vector<std::string> out;
    std::string tok;
    while (iss >> tok) out.push_back(tok);
    return out;
}

// ---------- numeric parsing ----------
static bool parseDoubleStrict(const std::string& tok, double& out) 
{
    try {
        std::size_t pos = 0;
        out = std::stod(tok, &pos);
        return pos == tok.size();
    } catch (...) { return false; }
}

static bool parseIntStrict(const std::string& tok, int& out) 
{
    try {
        std::size_t pos = 0;
        long long v = std::stoll(tok, &pos);
        if (pos != tok.size()) return false;
        if (v < std::numeric_limits<int>::min() || v > std::numeric_limits<int>::max()) return false;
        out = static_cast<int>(v);
        return true;
    } catch (...) { return false; }
}

// ---------- sections ----------
enum class Section 
{
    EdgeWeight,
    NodeCoord,
    ReturnRate,
    PassengersLimit,
    OriginsDestBudgets,
    DisplayData,
    Eof
};

static std::optional<Section> classifySectionLine(const std::string& lineTrimmedUpper) 
{
    if (lineTrimmedUpper == "EDGE_WEIGHT_SECTION") return Section::EdgeWeight;
    if (lineTrimmedUpper == "NODE_COORD_SECTION") return Section::NodeCoord;
    if (lineTrimmedUpper == "RETURN_RATE_SECTION") return Section::ReturnRate;
    if (lineTrimmedUpper == "PASSENGERS_LIMIT") return Section::PassengersLimit;
    if (lineTrimmedUpper == "ORIGINS_DESTINATIONS_AND_FINANCIAL_LIMITS") return Section::OriginsDestBudgets;
    if (lineTrimmedUpper == "DISPLAY_DATA_SECTION") return Section::DisplayData;
    if (lineTrimmedUpper == "EOF") return Section::Eof;
    return std::nullopt;
}

struct SectionSpan 
{
    Section kind{};
    int headerLine = -1; // line index of section header
    int endLine = -1;    // exclusive (next section header, or EOF line)
};

static std::vector<SectionSpan> computeSectionSpansOrThrow(const std::vector<std::string>& lines) {
    std::vector<std::pair<Section, int>> headers;
    headers.reserve(16);

    for (int i = 0; i < static_cast<int>(lines.size()); ++i) {
        const auto t = trim(lines[i]);
        if (t.empty()) continue;
        auto sec = classifySectionLine(upper(t));
        if (sec) headers.push_back({*sec, i});
    }
    if (headers.empty() || headers.back().first != Section::Eof)
        throw ParseError("Missing EOF marker");

    std::vector<SectionSpan> spans;
    spans.reserve(headers.size());
    for (std::size_t k = 0; k < headers.size(); ++k) {
        SectionSpan s;
        s.kind = headers[k].first;
        s.headerLine = headers[k].second;
        s.endLine = (k + 1 < headers.size()) ? headers[k + 1].second : static_cast<int>(lines.size());
        spans.push_back(s);
    }
    return spans;
}

static const SectionSpan* findSpan(const std::vector<SectionSpan>& spans, Section kind) {
    for (const auto& s : spans) if (s.kind == kind) return &s;
    return nullptr;
}

// ---------- header parsing ----------
static std::unordered_map<std::string, std::string> parseHeadersBeforeFirstSection(
    const std::vector<std::string>& lines,
    int firstSectionLine
) {
    std::unordered_map<std::string, std::string> h;
    for (int i = 0; i < firstSectionLine; ++i) {
        auto t = trim(lines[i]);
        if (t.empty()) continue;

        // KEY : VALUE (TSPLIB style). Also accept "KEY: VALUE".
        auto pos = t.find(':');
        if (pos == std::string::npos) continue;

        auto key = upper(trim(t.substr(0, pos)));
        auto val = trim(t.substr(pos + 1));
        if (!key.empty()) h[key] = val;
    }
    return h;
}

static ProblemType parseProblemTypeOrThrow(const std::unordered_map<std::string, std::string>& h) {
    auto it = h.find("TYPE");
    if (it == h.end()) throw ParseError("Missing TYPE header");
    const auto t = upper(trim(it->second));
    if (t == "TSP") return ProblemType::TSP;
    if (t == "ATSP") return ProblemType::ATSP;
    // CaRS/CaRSP are matched case-insensitively (already uppercased above)
    if (t == "CARS") return ProblemType::CaRS;
    if (t == "CARSP") return ProblemType::CaRSP;

    throw ParseError("Unsupported TYPE: " + it->second);
}

static int getIntHeaderOrThrow(const std::unordered_map<std::string, std::string>& h, const std::string& key) {
    auto it = h.find(key);
    if (it == h.end()) throw ParseError("Missing header: " + key);
    int v = 0;
    if (!parseIntStrict(trim(it->second), v))
        throw ParseError("Invalid integer for header " + key + ": " + it->second);
    return v;
}

static std::optional<bool> parseEdgeWeightSymmetricClaim(const std::unordered_map<std::string, std::string>& h) {
    auto it = h.find("EDGE_WEIGHT");
    if (it == h.end()) return std::nullopt;
    auto v = upper(trim(it->second));
    if (v == "SYMMETRIC") return true;
    if (v == "ASYMMETRIC") return false;
    return std::nullopt;
}

static std::optional<bool> parseReturnRateAsymmetricClaim(const std::unordered_map<std::string, std::string>& h) {
    auto it = h.find("RETURN_RATE");
    if (it == h.end()) return std::nullopt;
    auto v = upper(trim(it->second));
    if (v == "ASYMMETRIC") return true;
    if (v == "SYMMETRIC") return false;
    return std::nullopt;
}

static void validateSymmetricOrThrow(const DenseMatrix& m, const std::string& what) {
    for (int i = 0; i < m.n; ++i) {
        for (int j = i + 1; j < m.n; ++j) {
            if (m.at(i, j) != m.at(j, i)) {
                throw ParseError(what + " claimed SYMMETRIC but differs at (" +
                                 std::to_string(i) + "," + std::to_string(j) + ")");
            }
        }
    }
}

// ---------- TSPLIB EXPLICIT formats ----------
enum class ExplicitFormat {
    FullMatrix,
    UpperRow,
    LowerRow,
    UpperDiagRow,
    LowerDiagRow,
    UpperCol,
    LowerCol,
    UpperDiagCol,
    LowerDiagCol
};

static ExplicitFormat parseExplicitFormatOrThrow(const std::unordered_map<std::string, std::string>& h) {
    auto it = h.find("EDGE_WEIGHT_FORMAT");
    if (it == h.end()) throw ParseError("Missing EDGE_WEIGHT_FORMAT (EXPLICIT)");
    const auto f = upper(trim(it->second));
    if (f == "FULL_MATRIX") return ExplicitFormat::FullMatrix;
    if (f == "UPPER_ROW") return ExplicitFormat::UpperRow;
    if (f == "LOWER_ROW") return ExplicitFormat::LowerRow;
    if (f == "UPPER_DIAG_ROW") return ExplicitFormat::UpperDiagRow;
    if (f == "LOWER_DIAG_ROW") return ExplicitFormat::LowerDiagRow;
    if (f == "UPPER_COL") return ExplicitFormat::UpperCol;
    if (f == "LOWER_COL") return ExplicitFormat::LowerCol;
    if (f == "UPPER_DIAG_COL") return ExplicitFormat::UpperDiagCol;
    if (f == "LOWER_DIAG_COL") return ExplicitFormat::LowerDiagCol;
    throw ParseError("Unsupported EDGE_WEIGHT_FORMAT: " + it->second);
}

static std::vector<double> readNumericTokensInSpanOrThrow(
    const std::vector<std::string>& lines,
    int fromLine,
    int toLineExclusive
) {
    std::vector<double> vals;
    vals.reserve(1024);

    for (int i = fromLine; i < toLineExclusive; ++i) {
        const auto t = trim(lines[i]);
        if (t.empty()) continue;

        for (const auto& tok : splitWS(t)) {
            double x = 0.0;
            if (!parseDoubleStrict(tok, x)) {
                throw ParseError("Expected numeric token, got '" + tok +
                                 "' at line " + std::to_string(i + 1));
            }
            vals.push_back(x);
        }
    }
    return vals;
}

static DenseMatrix buildMatrixFromExplicitTokensOrThrow(
    int n,
    ExplicitFormat fmt,
    const std::vector<double>& v
) {
    DenseMatrix m(n);

    auto needCount = [&](ExplicitFormat f) -> std::size_t {
        const std::size_t N = static_cast<std::size_t>(n);
        switch (f) {
            case ExplicitFormat::FullMatrix:    return N * N;
            case ExplicitFormat::UpperRow:
            case ExplicitFormat::LowerRow:
            case ExplicitFormat::UpperCol:
            case ExplicitFormat::LowerCol:      return N * (N - 1) / 2;
            case ExplicitFormat::UpperDiagRow:
            case ExplicitFormat::LowerDiagRow:
            case ExplicitFormat::UpperDiagCol:
            case ExplicitFormat::LowerDiagCol:  return N * (N + 1) / 2;
        }
        return 0;
    };

    const auto expected = needCount(fmt);
    if (v.size() != expected) {
        throw ParseError("EDGE_WEIGHT_SECTION token count mismatch: expected " +
                         std::to_string(expected) + ", got " + std::to_string(v.size()));
    }

    // Initialize diagonal (common expectation for non-diag formats)
    for (int i = 0; i < n; ++i) m.at(i, i) = 0.0;

    std::size_t k = 0;
    const auto next = [&]() -> double { return v.at(k++); };

    switch (fmt) {
        case ExplicitFormat::FullMatrix: {
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    m.at(i, j) = next();
            break;
        }

        case ExplicitFormat::UpperRow: { // i<j by rows
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    m.at(i, j) = next();
                    m.at(j, i) = m.at(i, j);
                }
            }
            break;
        }
        case ExplicitFormat::LowerRow: { // i>j by rows
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < i; ++j) {
                    m.at(i, j) = next();
                    m.at(j, i) = m.at(i, j);
                }
            }
            break;
        }
        case ExplicitFormat::UpperDiagRow: { // i<=j by rows
            for (int i = 0; i < n; ++i) {
                for (int j = i; j < n; ++j) {
                    m.at(i, j) = next();
                    m.at(j, i) = m.at(i, j);
                }
            }
            break;
        }
        case ExplicitFormat::LowerDiagRow: { // i>=j by rows
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j <= i; ++j) {
                    m.at(i, j) = next();
                    m.at(j, i) = m.at(i, j);
                }
            }
            break;
        }

        case ExplicitFormat::UpperCol: { // i<j by cols
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < j; ++i) {
                    m.at(i, j) = next();
                    m.at(j, i) = m.at(i, j);
                }
            }
            break;
        }
        case ExplicitFormat::LowerCol: { // i>j by cols
            for (int j = 0; j < n; ++j) {
                for (int i = j + 1; i < n; ++i) {
                    m.at(i, j) = next();
                    m.at(j, i) = m.at(i, j);
                }
            }
            break;
        }
        case ExplicitFormat::UpperDiagCol: { // i<=j by cols
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i <= j; ++i) {
                    m.at(i, j) = next();
                    m.at(j, i) = m.at(i, j);
                }
            }
            break;
        }
        case ExplicitFormat::LowerDiagCol: { // i>=j by cols
            for (int j = 0; j < n; ++j) {
                for (int i = j; i < n; ++i) {
                    m.at(i, j) = next();
                    m.at(j, i) = m.at(i, j);
                }
            }
            break;
        }
    }

    return m;
}

// ---------- TSPLIB EUC_2D ----------
struct Coord { double x=0.0, y=0.0; };

static std::vector<Coord> readNodeCoordsOrThrow(
    int n,
    const std::vector<std::string>& lines,
    const SectionSpan& span
) {
    // Accept either "id x y" or "x y" (id implicit by order).
    std::vector<std::optional<Coord>> coords(static_cast<std::size_t>(n));
    int assigned = 0;
    int implicitId = 0;

    int minId = std::numeric_limits<int>::max();
    int maxId = std::numeric_limits<int>::min();
    bool sawExplicitId = false;

    // First pass: gather raw entries (id,x,y)
    struct Entry { int id; double x; double y; };
    std::vector<Entry> entries;

    for (int i = span.headerLine + 1; i < span.endLine; ++i) {
        const auto t = trim(lines[i]);
        if (t.empty()) continue;

        auto toks = splitWS(t);
        if (toks.size() == 3) {
            int id = 0;
            double x=0.0, y=0.0;
            if (!parseIntStrict(toks[0], id) ||
                !parseDoubleStrict(toks[1], x) ||
                !parseDoubleStrict(toks[2], y)) {
                throw ParseError("Invalid NODE_COORD_SECTION line at " + std::to_string(i + 1));
            }
            sawExplicitId = true;
            minId = std::min(minId, id);
            maxId = std::max(maxId, id);
            entries.push_back({id, x, y});
        } else if (toks.size() == 2) {
            double x=0.0, y=0.0;
            if (!parseDoubleStrict(toks[0], x) || !parseDoubleStrict(toks[1], y)) {
                throw ParseError("Invalid NODE_COORD_SECTION line at " + std::to_string(i + 1));
            }
            entries.push_back({implicitId++, x, y});
        } else {
            throw ParseError("Invalid NODE_COORD_SECTION token count at line " + std::to_string(i + 1));
        }
    }

    if (static_cast<int>(entries.size()) != n) {
        throw ParseError("NODE_COORD_SECTION count mismatch: expected " +
                         std::to_string(n) + ", got " + std::to_string(entries.size()));
    }

    // Determine mapping to 0..N-1:
    // - if explicit IDs and look like 1..N => idx = id-1
    // - if explicit IDs and look like 0..N-1 => idx = id
    // - otherwise: sort unique IDs and map by ascending order (strictly N unique)
    std::vector<int> ids;
    ids.reserve(entries.size());
    for (auto& e : entries) ids.push_back(e.id);
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    if (static_cast<int>(ids.size()) != n) {
        throw ParseError("NODE_COORD_SECTION IDs are not unique / not N entries");
    }

    auto idToIndex = [&](int id) -> int {
        if (!sawExplicitId) return id; // implicit already 0..n-1
        if (minId == 1 && maxId == n) return id - 1;
        if (minId == 0 && maxId == n - 1) return id;

        // fallback: ordered mapping
        auto it = std::lower_bound(ids.begin(), ids.end(), id);
        if (it == ids.end() || *it != id) return -1;
        return static_cast<int>(std::distance(ids.begin(), it));
    };

    for (auto& e : entries) {
        int idx = idToIndex(e.id);
        if (idx < 0 || idx >= n) {
            throw ParseError("NODE_COORD_SECTION ID out of range / cannot map to 0..N-1: " + std::to_string(e.id));
        }
        if (coords[static_cast<std::size_t>(idx)].has_value()) {
            throw ParseError("Duplicate node ID after mapping: " + std::to_string(e.id));
        }
        coords[static_cast<std::size_t>(idx)] = Coord{e.x, e.y};
        assigned++;
    }

    if (assigned != n) throw ParseError("NODE_COORD_SECTION internal assignment failed");

    std::vector<Coord> out(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        if (!coords[static_cast<std::size_t>(i)].has_value())
            throw ParseError("Missing coordinate for node index " + std::to_string(i));
        out[static_cast<std::size_t>(i)] = *coords[static_cast<std::size_t>(i)];
    }
    return out;
}

static DenseMatrix buildEuc2dRoundedMatrix(int n, const std::vector<Coord>& c) {
    DenseMatrix m(n);
    for (int i = 0; i < n; ++i) {
        m.at(i, i) = 0.0;
        for (int j = i + 1; j < n; ++j) {
            const double dx = c[static_cast<std::size_t>(i)].x - c[static_cast<std::size_t>(j)].x;
            const double dy = c[static_cast<std::size_t>(i)].y - c[static_cast<std::size_t>(j)].y;
            const double dist = std::sqrt(dx*dx + dy*dy);
            const double rounded = std::round(dist); // TSPLIB EUC_2D rule
            m.at(i, j) = rounded;
            m.at(j, i) = rounded;
        }
    }
    return m;
}

// ---------- CaRS/CaRSP car-block matrices ----------
static std::vector<DenseMatrix> readCarBlockMatricesOrThrow(
    int cars,
    int n,
    const std::vector<std::string>& lines,
    const SectionSpan& span,
    const std::string& sectionNameForErrors
) {
    std::vector<DenseMatrix> out;
    out.reserve(static_cast<std::size_t>(cars));

    int line = span.headerLine + 1;

    auto skipEmpty = [&]() {
        while (line < span.endLine && trim(lines[line]).empty()) line++;
    };

    for (int car = 0; car < cars; ++car) {
        skipEmpty();
        if (line >= span.endLine) {
            throw ParseError(sectionNameForErrors + ": unexpected end while expecting car id " + std::to_string(car));
        }

        // Car id line: can have whitespace; assume first token is integer id.
        auto toks = splitWS(trim(lines[line]));
        if (toks.empty()) {
            throw ParseError(sectionNameForErrors + ": missing car id at line " + std::to_string(line + 1));
        }
        int carId = -1;
        if (!parseIntStrict(toks[0], carId)) {
            throw ParseError(sectionNameForErrors + ": invalid car id token '" + toks[0] +
                             "' at line " + std::to_string(line + 1));
        }
        if (carId != car) {
            throw ParseError(sectionNameForErrors + ": expected car id " + std::to_string(car) +
                             ", got " + std::to_string(carId) + " at line " + std::to_string(line + 1));
        }
        line++;

        DenseMatrix m(n);
        const std::size_t need = static_cast<std::size_t>(n) * static_cast<std::size_t>(n);
        std::size_t k = 0;

        while (k < need) {
            if (line >= span.endLine) {
                throw ParseError(sectionNameForErrors + ": not enough matrix values for car " +
                                 std::to_string(car) + " (need " + std::to_string(need) +
                                 ", got " + std::to_string(k) + ")");
            }
            auto t = trim(lines[line]);
            line++;
            if (t.empty()) continue;

            for (const auto& tok : splitWS(t)) {
                double x = 0.0;
                if (!parseDoubleStrict(tok, x)) {
                    throw ParseError(sectionNameForErrors + ": expected numeric token, got '" + tok +
                                     "' at line " + std::to_string(line));
                }
                if (k >= need) {
                    throw ParseError(sectionNameForErrors + ": too many matrix values for car " + std::to_string(car));
                }
                m.a[k++] = x;
            }
        }

        out.push_back(std::move(m));
    }

    // Strictness: only allow whitespace until next section header.
    while (line < span.endLine) {
        if (!trim(lines[line]).empty()) {
            throw ParseError(sectionNameForErrors + ": unexpected extra content at line " + std::to_string(line + 1));
        }
        line++;
    }

    return out;
}

static std::vector<int> readKIntsOrThrow(
    int kNeed,
    const std::vector<std::string>& lines,
    const SectionSpan& span,
    const std::string& sectionNameForErrors
) {
    std::vector<int> out;
    out.reserve(static_cast<std::size_t>(kNeed));

    for (int i = span.headerLine + 1; i < span.endLine; ++i) {
        const auto t = trim(lines[i]);
        if (t.empty()) continue;

        for (const auto& tok : splitWS(t)) {
            int v = 0;
            if (!parseIntStrict(tok, v)) {
                throw ParseError(sectionNameForErrors + ": expected integer token, got '" + tok +
                                 "' at line " + std::to_string(i + 1));
            }
            out.push_back(v);
            if (static_cast<int>(out.size()) > kNeed) {
                throw ParseError(sectionNameForErrors + ": too many integers (need " +
                                 std::to_string(kNeed) + ")");
            }
        }
    }

    if (static_cast<int>(out.size()) != kNeed) {
        throw ParseError(sectionNameForErrors + ": count mismatch (need " +
                         std::to_string(kNeed) + ", got " + std::to_string(out.size()) + ")");
    }
    return out;
}

static std::vector<PassengerRequest> readPassengerTriplesOrThrow(
    int pNeed,
    int n,
    const std::vector<std::string>& lines,
    const SectionSpan& span
) {
    // Parse numeric tokens and group them as (origin, destination, budget).
    // We require exactly pNeed * 3 tokens (so EOF early => error, extra data => error).
    std::vector<double> tokens = readNumericTokensInSpanOrThrow(lines, span.headerLine + 1, span.endLine);
    const std::size_t needTok = static_cast<std::size_t>(pNeed) * 3;

    if (tokens.size() != needTok) {
        throw ParseError("ORIGINS_DESTINATIONS_AND_FINANCIAL_LIMITS token count mismatch: expected " +
                         std::to_string(needTok) + ", got " + std::to_string(tokens.size()));
    }

    std::vector<PassengerRequest> out;
    out.reserve(static_cast<std::size_t>(pNeed));

    for (int p = 0; p < pNeed; ++p) {
        const double oD = tokens[static_cast<std::size_t>(p) * 3 + 0];
        const double dD = tokens[static_cast<std::size_t>(p) * 3 + 1];
        const double b  = tokens[static_cast<std::size_t>(p) * 3 + 2];

        // Origins/destinations should be integers in the file. Enforce "integer-like".
        const double oRound = std::round(oD);
        const double dRound = std::round(dD);
        if (oD != oRound || dD != dRound) {
            throw ParseError("Passenger origin/destination must be integer-like, got (" +
                             std::to_string(oD) + ", " + std::to_string(dD) + ")");
        }

        const int o = static_cast<int>(oRound);
        const int d = static_cast<int>(dRound);

        if (o < 0 || o >= n || d < 0 || d >= n) {
            throw ParseError("Passenger origin/destination out of range: (" +
                             std::to_string(o) + ", " + std::to_string(d) + "), N=" + std::to_string(n));
        }
        if (b < 0.0) {
            throw ParseError("Passenger budget must be >= 0, got " + std::to_string(b));
        }

        out.push_back(PassengerRequest{o, d, b});
    }

    return out;
}

// ---------- high-level parsing per TYPE ----------
static std::string getHeaderString(const std::unordered_map<std::string, std::string>& h,
                                   const std::string& key) {
    auto it = h.find(key);
    return (it == h.end()) ? std::string{} : it->second;
}

static std::string requireHeaderString(const std::unordered_map<std::string, std::string>& h,
                                       const std::string& key) {
    auto it = h.find(key);
    if (it == h.end()) throw ParseError("Missing header: " + key);
    return it->second;
}

static Instance parseTSPLIBOrThrow(
    const std::unordered_map<std::string, std::string>& h,
    const std::vector<std::string>& lines,
    const std::vector<SectionSpan>& spans,
    const ParserOptions& opt
) {
    Instance inst;
    inst.name = getHeaderString(h, "NAME");
    inst.comment = getHeaderString(h, "COMMENT");
    inst.type = parseProblemTypeOrThrow(h);
    inst.dimension = getIntHeaderOrThrow(h, "DIMENSION");

    inst.carsNumber = 1;
    inst.passengersNumber = 1;

    // TSPLIB must not have return rate section (per your requirement).
    if (findSpan(spans, Section::ReturnRate) != nullptr) {
        throw ParseError("TSPLIB instance must not contain RETURN_RATE_SECTION");
    }

    // EDGE_WEIGHT_TYPE
    const auto ewt = upper(trim(requireHeaderString(h, "EDGE_WEIGHT_TYPE")));

    DenseMatrix travel(inst.dimension);

    if (ewt == "EUC_2D") {
        if (!opt.tsplibComputeEuc2d) {
            throw ParseError("TSPLIB EUC_2D disabled by options");
        }
        const auto* nodeSpan = findSpan(spans, Section::NodeCoord);
        if (!nodeSpan) throw ParseError("TSPLIB EUC_2D requires NODE_COORD_SECTION");

        // If an explicit EDGE_WEIGHT_SECTION also exists, treat as error (avoid ambiguity).
        if (findSpan(spans, Section::EdgeWeight) != nullptr) {
            throw ParseError("TSPLIB EUC_2D instance contains EDGE_WEIGHT_SECTION (ambiguous)");
        }

        auto coords = readNodeCoordsOrThrow(inst.dimension, lines, *nodeSpan);
        travel = buildEuc2dRoundedMatrix(inst.dimension, coords);
    }
    else if (ewt == "EXPLICIT") {
        const auto* ewSpan = findSpan(spans, Section::EdgeWeight);
        if (!ewSpan) throw ParseError("TSPLIB EXPLICIT requires EDGE_WEIGHT_SECTION");

        auto fmt = parseExplicitFormatOrThrow(h);

        // For ATSP, triangular formats are generally not valid; require FULL_MATRIX.
        if (inst.type == ProblemType::ATSP && fmt != ExplicitFormat::FullMatrix) {
            throw ParseError("ATSP requires EDGE_WEIGHT_FORMAT FULL_MATRIX (triangular formats not allowed)");
        }

        if (!opt.tsplibSupportTriangular && fmt != ExplicitFormat::FullMatrix) {
            throw ParseError("Triangular TSPLIB formats disabled by options");
        }

        auto tokens = readNumericTokensInSpanOrThrow(lines, ewSpan->headerLine + 1, ewSpan->endLine);
        travel = buildMatrixFromExplicitTokensOrThrow(inst.dimension, fmt, tokens);
    }
    else {
        throw ParseError("Unsupported TSPLIB EDGE_WEIGHT_TYPE: " + ewt);
    }

    inst.travelCostPerCar.clear();
    inst.travelCostPerCar.push_back(std::move(travel));

    // Optional symmetry claim validation (only if tag exists).
    inst.edgeWeightIsSymmetric = parseEdgeWeightSymmetricClaim(h);
    if (opt.validateSymmetryWhenClaimed && inst.edgeWeightIsSymmetric.has_value() && *inst.edgeWeightIsSymmetric) {
        validateSymmetricOrThrow(inst.travelCostPerCar[0], "EDGE_WEIGHT");
    }

    return inst;
}

static Instance parseCaRSOrThrow(
    const std::unordered_map<std::string, std::string>& h,
    const std::vector<std::string>& lines,
    const std::vector<SectionSpan>& spans,
    const ParserOptions& opt,
    bool isCaRSP
) {
    Instance inst;
    inst.name = getHeaderString(h, "NAME");
    inst.comment = getHeaderString(h, "COMMENT");
    inst.type = parseProblemTypeOrThrow(h);
    inst.dimension = getIntHeaderOrThrow(h, "DIMENSION");
    inst.carsNumber = getIntHeaderOrThrow(h, "CARS_NUMBER");

    // CaRS/CaRSP only supports EXPLICIT + FULL_MATRIX per spec.
    const auto ewt = upper(trim(requireHeaderString(h, "EDGE_WEIGHT_TYPE")));
    if (ewt != "EXPLICIT") {
        throw ParseError("CaRS/CaRSP only supports EDGE_WEIGHT_TYPE EXPLICIT (got " + ewt + ")");
    }

    const auto fmt = upper(trim(requireHeaderString(h, "EDGE_WEIGHT_FORMAT")));
    if (fmt != "FULL_MATRIX") {
        throw ParseError("CaRS/CaRSP only supports EDGE_WEIGHT_FORMAT FULL_MATRIX (got " + fmt + ")");
    }

    const auto* ewSpan = findSpan(spans, Section::EdgeWeight);
    if (!ewSpan) throw ParseError("Missing EDGE_WEIGHT_SECTION");

    const auto* rrSpan = findSpan(spans, Section::ReturnRate);
    if (!rrSpan) throw ParseError("Missing RETURN_RATE_SECTION");

    inst.travelCostPerCar = readCarBlockMatricesOrThrow(inst.carsNumber, inst.dimension, lines, *ewSpan, "EDGE_WEIGHT_SECTION");
    inst.returnCostPerCar = readCarBlockMatricesOrThrow(inst.carsNumber, inst.dimension, lines, *rrSpan, "RETURN_RATE_SECTION");

    // Optional claims (validate only if present).
    inst.edgeWeightIsSymmetric = parseEdgeWeightSymmetricClaim(h);
    inst.returnRateIsAsymmetric = parseReturnRateAsymmetricClaim(h);

    if (opt.validateSymmetryWhenClaimed && inst.edgeWeightIsSymmetric.has_value() && *inst.edgeWeightIsSymmetric) {
        for (int c = 0; c < inst.carsNumber; ++c) {
            validateSymmetricOrThrow(inst.travelCostPerCar[static_cast<std::size_t>(c)], "EDGE_WEIGHT");
        }
    }

    if (!isCaRSP) {
        inst.passengersNumber = 1; // interface-unified
        return inst;
    }

    // --- CaRSP extras ---
    inst.passengersNumber = getIntHeaderOrThrow(h, "PASSENGERS_NUMBER");

    const auto* plSpan = findSpan(spans, Section::PassengersLimit);
    if (!plSpan) throw ParseError("Missing PASSENGERS_LIMIT section for CaRSP");

    const auto* odSpan = findSpan(spans, Section::OriginsDestBudgets);
    if (!odSpan) throw ParseError("Missing ORIGINS_DESTINATIONS_AND_FINANCIAL_LIMITS section for CaRSP");

    inst.passengersLimitPerCar = readKIntsOrThrow(inst.carsNumber, lines, *plSpan, "PASSENGERS_LIMIT");
    inst.passengers = readPassengerTriplesOrThrow(inst.passengersNumber, inst.dimension, lines, *odSpan);

    return inst;
}

} // namespace

// ---------- Parser public API ----------
Parser::Parser(ParserOptions opt) : opt_(opt) {}

Instance Parser::parseFileSpec(const std::string& spec) const {
    std::string path = spec;
    std::replace(path.begin(), path.end(), '?', '/');

    std::ifstream in(path);
    if (!in) throw ParseError("Cannot open file: " + path);

    return parseStream(in);
}

Instance Parser::parseStream(std::istream& in) const {
    std::vector<std::string> lines;
    lines.reserve(2048);

    std::string line;
    while (std::getline(in, line)) {
        // Keep raw lines; we trim when needed.
        lines.push_back(line);
    }
    if (lines.empty()) throw ParseError("Empty input");

    const auto spans = computeSectionSpansOrThrow(lines);

    // Determine where headers stop: first section header line (excluding EOF).
    int firstSectionLine = std::numeric_limits<int>::max();
    for (const auto& s : spans) {
        if (s.kind == Section::Eof) continue;
        firstSectionLine = std::min(firstSectionLine, s.headerLine);
    }
    if (firstSectionLine == std::numeric_limits<int>::max()) {
        throw ParseError("No sections found before EOF");
    }

    const auto headers = parseHeadersBeforeFirstSection(lines, firstSectionLine);
    const auto type = parseProblemTypeOrThrow(headers);

    // Route to correct parser.
    if (type == ProblemType::TSP || type == ProblemType::ATSP) {
        // Ensure CaRS/CaRSP-only sections are not present:
        if (findSpan(spans, Section::PassengersLimit) != nullptr ||
            findSpan(spans, Section::OriginsDestBudgets) != nullptr) {
            throw ParseError("TSPLIB instance contains CaRSP-only sections");
        }
        return parseTSPLIBOrThrow(headers, lines, spans, opt_);
    }

    if (type == ProblemType::CaRS) {
        // CaRS must not contain CaRSP-only passenger sections.
        if (findSpan(spans, Section::PassengersLimit) != nullptr ||
            findSpan(spans, Section::OriginsDestBudgets) != nullptr) {
            throw ParseError("CaRS instance contains CaRSP-only passenger sections");
        }
        return parseCaRSOrThrow(headers, lines, spans, opt_, /*isCaRSP=*/false);
    }

    if (type == ProblemType::CaRSP) {
        return parseCaRSOrThrow(headers, lines, spans, opt_, /*isCaRSP=*/true);
    }

    throw ParseError("Unsupported problem TYPE");
}

void printMatrix(std::ostream& os, const DenseMatrix& m, const std::string& title) 
{
    if (!title.empty()) os << title << "\n";
    os << "n=" << m.n << "\n";
    os << std::fixed << std::setprecision(0);

    for (int i = 0; i < m.n; ++i) 
    {
        for (int j = 0; j < m.n; ++j) 
        {
            os << m.at(i, j);
            if (j + 1 < m.n) os << ' ';
        }
        os << "\n";
    }
}

static const char* toString(ProblemType t) 
{
    switch (t) 
    {
        case ProblemType::TSP:  return "TSP";
        case ProblemType::ATSP: return "ATSP";
        case ProblemType::CaRS: return "CaRS";
        case ProblemType::CaRSP:return "CaRSP";
    }
    return "UNKNOWN";
}

void printInstance(std::ostream& os, const Instance& inst, bool printMatrices) 
{
    os << "NAME: " << inst.name << "\n";
    os << "TYPE: " << toString(inst.type) << "\n";
    if (!inst.comment.empty()) os << "COMMENT: " << inst.comment << "\n";
    os << "DIMENSION: " << inst.dimension << "\n";
    os << "CARS_NUMBER: " << inst.carsNumber << "\n";
    os << "PASSENGERS_NUMBER (interface): " << inst.passengersNumber << "\n";

    if (inst.edgeWeightIsSymmetric.has_value()) 
    {
        os << "EDGE_WEIGHT claim: " << (*inst.edgeWeightIsSymmetric ? "SYMMETRIC" : "ASYMMETRIC") << "\n";
    }
    if (inst.returnRateIsAsymmetric.has_value()) 
    {
        os << "RETURN_RATE claim: " << (*inst.returnRateIsAsymmetric ? "ASYMMETRIC" : "SYMMETRIC") << "\n";
    }

    os << "Travel matrices: " << inst.travelCostPerCar.size() << "\n";
    os << "Return matrices: " << inst.returnCostPerCar.size() << "\n";
    os << "Passengers (explicit): " << inst.passengers.size() << "\n";
    os << "Capacities: " << inst.passengersLimitPerCar.size() << "\n";

    if (inst.isCaRSP()) 
    {
        if (!inst.passengersLimitPerCar.empty()) 
        {
            os << "PASSENGERS_LIMIT per car: ";
            for (std::size_t i = 0; i < inst.passengersLimitPerCar.size(); ++i) {
                os << inst.passengersLimitPerCar[i];
                if (i + 1 < inst.passengersLimitPerCar.size()) os << ' ';
            }
            os << "\n";
        }

        if (!inst.passengers.empty()) 
        {
            os << "ORIGINS_DESTINATIONS_AND_FINANCIAL_LIMITS:\n";
            os << std::fixed << std::setprecision(6);
            for (const auto& p : inst.passengers) {
                os << p.origin << " " << p.destination << " " << p.budget << "\n";
            }
        }
    }

    if (!printMatrices) return;

    for (std::size_t c = 0; c < inst.travelCostPerCar.size(); ++c) 
    {
        printMatrix(os, inst.travelCostPerCar[c], "EDGE_WEIGHT car " + std::to_string(c));
    }
    for (std::size_t c = 0; c < inst.returnCostPerCar.size(); ++c) 
    {
        printMatrix(os, inst.returnCostPerCar[c], "RETURN_RATE car " + std::to_string(c));
    }
}

std::ostream& operator<<(std::ostream& os, const Instance& inst) 
{
    // Default: summary only
    printInstance(os, inst, /*printMatrices=*/true);
    return os;
}

} // namespace cars_tsplib
