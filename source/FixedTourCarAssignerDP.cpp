#include <FixedTourCarAssignerDP.hpp>
#include <parser.hpp>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>
#include <cmath>

namespace localSearch 
{

FixedTourCarAssignerDP::FixedTourCarAssignerDP(std::shared_ptr<const cars_tsplib::Instance> inst,
                                               CarAssignmentDPOptions opt)
    : inst_(std::move(inst)), opt_(opt)
{
    if (!inst_) throw std::runtime_error("FixedTourCarAssignerDP: inst is null");
}

void FixedTourCarAssignerDP::validateTourOrThrow_(const std::vector<int>& nodes) const
{
    const int N = inst_->n();
    if (static_cast<int>(nodes.size()) != N)
        throw std::runtime_error("FixedTourCarAssignerDP: nodes.size() != inst->n()");

    if (N <= 0) return;

    if (nodes[0] != 0)
        throw std::runtime_error("FixedTourCarAssignerDP: node[0] must be 0 (fixed start).");

    // permutation check (optional)
    std::vector<char> seen(static_cast<std::size_t>(N), 0);
    for (int k = 0; k < N; ++k)
    {
        const int vv = nodes[static_cast<std::size_t>(k)];
        if (vv < 0 || vv >= N)
            throw std::runtime_error("FixedTourCarAssignerDP: node value out of range.");
        if (seen[static_cast<std::size_t>(vv)])
            throw std::runtime_error("FixedTourCarAssignerDP: duplicate node in tour.");
        seen[static_cast<std::size_t>(vv)] = 1;
    }
}

int FixedTourCarAssignerDP::normalizeMaxSegLen_(int L) const
{
    const int N = inst_->n();
    const int C = inst_->cars();

    if (L <= 0) return -1;   // <=0 means "no limit"
    if (N <= 0) return -1;

    // Necessary feasibility condition with a limit: N edges must be covered by at most
    // C segments, each of length at most L => C*L >= N.
    if (C > 0)
    {
        const int minL = (N + C - 1) / C; // ceil(N/C)
        if (L < minL) L = minL;
    }

    if (L >= N) return -1;   // effectively unlimited
    return L;
}


void Scratch::ensure(int N_, int C_)
{
    if (N == N_ && C == C_) return;

    N = N_;
    C = C_;
    if (C <= 0) throw std::runtime_error("Scratch::ensure: C<=0");
    if (C > 20) throw std::runtime_error("Scratch::ensure: C>20 not supported");

    M = 1 << C;

    pref.assign((std::size_t)C * (std::size_t)(N + 1), 0.0);
    dp.assign((std::size_t)(N + 1) * (std::size_t)M, 0.0);

    tour.resize((std::size_t)N + 1);
    seenStamp.assign((std::size_t)N, 0);
    stamp = 1;
}


CarAssignmentDPResult FixedTourCarAssignerDP::optimizeWithMaxSegLen_(const std::vector<int>& nodes,
                                                                    int maxSegmentLen) const
{
    const int N = inst_->n();
    const int C = inst_->cars();

    CarAssignmentDPResult res;

    if (N <= 0)
    {
        res.cost = 0.0;
        res.carPerEdge.clear();
        return res;
    }
    if (C <= 0) throw std::runtime_error("FixedTourCarAssignerDP: inst->cars() <= 0");

    if (C > 20)
        throw std::runtime_error("FixedTourCarAssignerDP: too many cars for bitmask DP (C > 20).");

    // Fast path: pure TSP (C==1) => no switching and no returnCost.
    if (C == 1)
    {
        double sum = 0.0;
        for (int k = 0; k < N; ++k)
        {
            const int a = nodes[(std::size_t)k];
            const int b = (k + 1 < N) ? nodes[(std::size_t)(k + 1)] : nodes[0];
            sum += inst_->travelCost(0, a, b);
        }
        res.cost = sum;
        res.carPerEdge.assign((std::size_t)N, 0);
        return res;
    }

    maxSegmentLen = normalizeMaxSegLen_(maxSegmentLen);

    const int M = 1 << C;
    const double INF = std::numeric_limits<double>::infinity();

    // v[0..N], closure v[N]=v[0]
    std::vector<int> v(static_cast<std::size_t>(N + 1));
    for (int i = 0; i < N; ++i) v[static_cast<std::size_t>(i)] = nodes[static_cast<std::size_t>(i)];
    v[static_cast<std::size_t>(N)] = nodes[0];

    // pref[c][t] = sum_{k=0..t-1} travelCost(c, v[k], v[k+1])
    std::vector<std::vector<double>> pref(static_cast<std::size_t>(C),
                                          std::vector<double>(static_cast<std::size_t>(N + 1), 0.0));
    for (int c = 0; c < C; ++c)
    {
        for (int k = 0; k < N; ++k)
        {
            const int a = v[static_cast<std::size_t>(k)];
            const int b = v[static_cast<std::size_t>(k + 1)];
            pref[static_cast<std::size_t>(c)][static_cast<std::size_t>(k + 1)] =
                pref[static_cast<std::size_t>(c)][static_cast<std::size_t>(k)] + inst_->travelCost(c, a, b);
        }
    }

    auto segCost = [&](int c, int s, int t) -> double
    {
        const double travel =
            pref[static_cast<std::size_t>(c)][static_cast<std::size_t>(t)] -
            pref[static_cast<std::size_t>(c)][static_cast<std::size_t>(s)];

        const int endNode   = v[static_cast<std::size_t>(t)];
        const int startNode = v[static_cast<std::size_t>(s)];

        const double ret = (C == 1) ? 0.0 : inst_->returnCost(c, endNode, startNode);
        return travel + ret;

    };

    // dp[t][mask] = minimal cost to cover edges 0..t-1,
    // using cars from mask (each car at most once).
    std::vector<double> dp(static_cast<std::size_t>((N + 1) * M), INF);

    struct Parent
    {
        int prevT = -1;
        int prevMask = -1;
        int car = -1;
    };
    std::vector<Parent> parent(static_cast<std::size_t>((N + 1) * M));

    auto idx = [&](int t, int mask) -> std::size_t
    {
        return static_cast<std::size_t>(t * M + mask);
    };

    dp[idx(0, 0)] = 0.0;

    for (int s = 0; s < N; ++s)
    {
        int tMax = N;
        if (maxSegmentLen > 0)
            tMax = std::min(N, s + maxSegmentLen);

        for (int mask = 0; mask < M; ++mask)
        {
            const double base = dp[idx(s, mask)];
            if (base == INF) continue;

            // pick a "new" car not yet used
            for (int c = 0; c < C; ++c)
            {
                const int bit = 1 << c;
                if (mask & bit) continue;

                const int nmask = mask | bit;

                // pick the end of the segment
                for (int t = s + 1; t <= tMax; ++t)
                {
                    const double cand = base + segCost(c, s, t);
                    const std::size_t it = idx(t, nmask);
                    if (cand < dp[it])
                    {
                        dp[it] = cand;
                        Parent p;
                        p.prevT = s;
                        p.prevMask = mask;
                        p.car = c;
                        parent[it] = p;
                    }
                }
            }
        }
    }

    // best ending: dp[N][mask] for any mask
    double bestCost = INF;
    int bestMask = -1;
    for (int mask = 0; mask < M; ++mask)
    {
        const double val = dp[idx(N, mask)];
        if (val < bestCost)
        {
            bestCost = val;
            bestMask = mask;
        }
    }

    res.cost = bestCost;
    res.carPerEdge.assign(static_cast<std::size_t>(N), 0);

    if (bestCost == INF)
    {
        // infeasible (could be due to maxSegmentLen or a genuine lack of cars)
        return res;
    }

    // Reconstruction: walk backward from (N, bestMask), filling carPerEdge on [s,t)
    int t = N;
    int mask = bestMask;
    while (t > 0)
    {
        const Parent p = parent[idx(t, mask)];
        if (p.prevT < 0 || p.prevMask < 0 || p.car < 0)
            throw std::runtime_error("FixedTourCarAssignerDP: DP reconstruction failed.");

        const int s = p.prevT;
        const int c = p.car;

        for (int k = s; k < t; ++k)
            res.carPerEdge[static_cast<std::size_t>(k)] = c;

        t = s;
        mask = p.prevMask;
    }

    return res;
}

CarAssignmentDPResult FixedTourCarAssignerDP::optimizeWithFallback_(const std::vector<int>& nodes) const
{
    const int N = inst_->n();
    const double INF = std::numeric_limits<double>::infinity();

    // If maxSegmentLen <=0 or fallback is disabled, run only once.
    if (!opt_.fallbackToUnlimitedOnInf)
        return optimizeWithMaxSegLen_(nodes, opt_.maxSegmentLen);

    int L = opt_.maxSegmentLen;

    // First attempt with the user-supplied L (may be -1)
    {
        auto res = optimizeWithMaxSegLen_(nodes, L);
        if (std::isfinite(res.cost) || L <= 0) return res;
    }

    // Retry with a larger L
    L = std::max(1, L);
    for (int r = 0; r < std::max(0, opt_.fallbackMaxRetries); ++r)
    {
        if (L >= N) break;
        const int nextL = std::min(N, L * std::max(2, opt_.fallbackGrowthFactor));
        if (nextL == L) break;
        L = nextL;

        auto res = optimizeWithMaxSegLen_(nodes, L);
        if (std::isfinite(res.cost)) return res;
    }

    // Final: safe -1 (exact optimum, slowest)
    {
        auto res = optimizeWithMaxSegLen_(nodes, -1);
        return res; // may be INF if truly infeasible (too few cars)
    }
}

CarAssignmentDPResult FixedTourCarAssignerDP::optimize(const std::vector<int>& nodes) const
{
    if (opt_.validateTour)
    {
        validateTourOrThrow_(nodes);
    }
    else
    {
        // minimal checks
        const int N = inst_->n();
        if (static_cast<int>(nodes.size()) != N)
            throw std::runtime_error("FixedTourCarAssignerDP: nodes.size() != inst->n()");
        if (N > 0 && nodes[0] != 0)
            throw std::runtime_error("FixedTourCarAssignerDP: node[0] must be 0 (fixed start).");
    }

    return optimizeWithFallback_(nodes);
}

double FixedTourCarAssignerDP::reassignCars(const std::vector<int>& nodes,
                                            std::vector<int>& carPerEdgeOut) const
{
    auto res = optimize(nodes);
    carPerEdgeOut = std::move(res.carPerEdge);
    return res.cost;
}

// -------
// cost-only DP, view variant

double FixedTourCarAssignerDP::evaluateCostViewWithMaxSegLen_(const std::function<int(int)>& getNode,
                                                             int maxSegmentLen) const
{
    const int N = inst_->n();
    const int C = inst_->cars();

    if (N <= 0) return 0.0;
    if (C <= 0) throw std::runtime_error("FixedTourCarAssignerDP: inst->cars() <= 0");
    if (C > 20) throw std::runtime_error("FixedTourCarAssignerDP: too many cars for bitmask DP (C > 20).");

    // Fast path: pure TSP (C==1) => no switching and no returnCost.
    if (C == 1)
    {
        double sum = 0.0;
        for (int k = 0; k < N; ++k)
        {
            const int a = getNode(k);
            const int b = (k + 1 < N) ? getNode(k + 1) : getNode(0);
            sum += inst_->travelCost(0, a, b);
        }
        return sum;
    }


    maxSegmentLen = normalizeMaxSegLen_(maxSegmentLen);

    const int M = 1 << C;
    const double INF = std::numeric_limits<double>::infinity();

    auto nodeAt = [&](int idx) -> int {
        // closure: v[N] = v[0]
        if (idx == N) return getNode(0);
        return getNode(idx);
    };

    // pref[c][t] = sum_{k=0..t-1} travelCost(c, v[k], v[k+1])
    std::vector<std::vector<double>> pref((std::size_t)C,
                                          std::vector<double>((std::size_t)N + 1, 0.0));
    for (int c = 0; c < C; ++c)
    {
        for (int k = 0; k < N; ++k)
        {
            const int a = nodeAt(k);
            const int b = (k + 1 < N) ? getNode(k + 1) : getNode(0);
            pref[(std::size_t)c][(std::size_t)k + 1] =
                pref[(std::size_t)c][(std::size_t)k] + inst_->travelCost(c, a, b);
        }
    }

    auto segCost = [&](int c, int s, int t) -> double
    {
        const double travel =
            pref[(std::size_t)c][(std::size_t)t] - pref[(std::size_t)c][(std::size_t)s];

        const int endNode   = nodeAt(t);
        const int startNode = nodeAt(s);

        const double ret = (C == 1) ? 0.0 : inst_->returnCost(c, endNode, startNode);
        return travel + ret;

    };

    std::vector<double> dp((std::size_t)(N + 1) * M, INF);

    auto idx = [&](int t, int mask) -> std::size_t {
        return (std::size_t)(t * M + mask);
    };

    dp[idx(0, 0)] = 0.0;

    for (int s = 0; s < N; ++s)
    {
        int tMax = N;
        if (maxSegmentLen > 0)
            tMax = std::min(N, s + maxSegmentLen);

        for (int mask = 0; mask < M; ++mask)
        {
            const double base = dp[idx(s, mask)];
            if (base == INF) continue;

            for (int c = 0; c < C; ++c)
            {
                const int bit = 1 << c;
                if (mask & bit) continue;

                const int nmask = mask | bit;

                for (int t = s + 1; t <= tMax; ++t)
                {
                    const double cand = base + segCost(c, s, t);
                    const std::size_t it = idx(t, nmask);
                    if (cand < dp[it]) dp[it] = cand;
                }
            }
        }
    }

    double bestCost = INF;
    for (int mask = 0; mask < M; ++mask)
        bestCost = std::min(bestCost, dp[idx(N, mask)]);

    return bestCost;
}

double FixedTourCarAssignerDP::evaluateCostViewWithMaxSegLenScratch_(
    const std::function<int(int)>& getNode,
    int maxSegmentLen,
    Scratch& scratch) const
{
    const int N = inst_->n();
    const int C = inst_->cars();
    if (N <= 0) return 0.0;
    if (C <= 0) throw std::runtime_error("FixedTourCarAssignerDP: inst->cars() <= 0");
    if (C > 20) throw std::runtime_error("FixedTourCarAssignerDP: too many cars for bitmask DP (C > 20).");

    // Fast path: pure TSP (C==1) => no switching and no returnCost.
    if (C == 1)
    {
        double sum = 0.0;
        for (int k = 0; k < N; ++k)
        {
            const int a = getNode(k);
            const int b = getNode(k + 1); // expects closure: getNode(N) == getNode(0)
            sum += inst_->travelCost(0, a, b);
        }
        return sum;
    }


    maxSegmentLen = normalizeMaxSegLen_(maxSegmentLen);

    scratch.ensure(N, C);

    const int M = scratch.M;
    const double INF = std::numeric_limits<double>::infinity();

    auto nodeAt = [&](int idx) -> int {
        return getNode(idx); // expects idx to range over 0..N
    };

    auto prefAt = [&](int c, int t) -> double& {
        return scratch.pref[(std::size_t)c * (std::size_t)(N + 1) + (std::size_t)t];
    };

    // Fill pref (no need to fill it entirely first, since it is fully overwritten)
    for (int c = 0; c < C; ++c)
    {
        prefAt(c, 0) = 0.0;
        for (int k = 0; k < N; ++k)
        {
            const int a = nodeAt(k);
            const int b = nodeAt(k + 1);
            prefAt(c, k + 1) = prefAt(c, k) + inst_->travelCost(c, a, b);
        }
    }

    auto segCost = [&](int c, int s, int t) -> double
    {
        const double travel = prefAt(c, t) - prefAt(c, s);
        const int endNode   = nodeAt(t);
        const int startNode = nodeAt(s);
        return travel + inst_->returnCost(c, endNode, startNode);
    };

    // Fill dp with INF
    std::fill(scratch.dp.begin(), scratch.dp.end(), INF);

    auto idx = [&](int t, int mask) -> std::size_t {
        return (std::size_t)t * (std::size_t)M + (std::size_t)mask;
    };

    scratch.dp[idx(0, 0)] = 0.0;

    for (int s = 0; s < N; ++s)
    {
        int tMax = N;
        if (maxSegmentLen > 0)
            tMax = std::min(N, s + maxSegmentLen);

        for (int mask = 0; mask < M; ++mask)
        {
            const double base = scratch.dp[idx(s, mask)];
            if (base == INF) continue;


            for (int c = 0; c < C; ++c)
            {
                const int bit = 1 << c;
                if (mask & bit) continue;

                const int nmask = mask | bit;

                for (int t = s + 1; t <= tMax; ++t)
                {
                    const double cand = base + segCost(c, s, t);
                    const std::size_t it = idx(t, nmask);
                    if (cand < scratch.dp[it]) scratch.dp[it] = cand;
                }
            }
        }
    }

    double bestCost = INF;
    for (int mask = 0; mask < M; ++mask)
        bestCost = std::min(bestCost, scratch.dp[idx(N, mask)]);

    return bestCost;
}


double FixedTourCarAssignerDP::evaluateCostViewWithFallback_(const std::function<int(int)>& getNode) const
{
    const int N = inst_->n();
    const double INF = std::numeric_limits<double>::infinity();

    if (!opt_.fallbackToUnlimitedOnInf)
        return evaluateCostViewWithMaxSegLen_(getNode, opt_.maxSegmentLen);

    int L = opt_.maxSegmentLen;

    // First attempt
    {
        const double c0 = evaluateCostViewWithMaxSegLen_(getNode, L);
        if (std::isfinite(c0) || L <= 0) return c0;
    }

    // Retry with a larger L
    L = std::max(1, L);
    for (int r = 0; r < std::max(0, opt_.fallbackMaxRetries); ++r)
    {
        if (L >= N) break;
        const int nextL = std::min(N, L * std::max(2, opt_.fallbackGrowthFactor));
        if (nextL == L) break;
        L = nextL;

        const double cc = evaluateCostViewWithMaxSegLen_(getNode, L);
        if (std::isfinite(cc)) return cc;
    }

    // Final: safe -1
    return evaluateCostViewWithMaxSegLen_(getNode, -1);
}

double FixedTourCarAssignerDP::evaluateCostViewScratch(const std::function<int(int)>& getNode,
                                                       Scratch& scratch) const
{
    const int N = inst_->n();
    const int C = inst_->cars();
    if (N <= 0) return 0.0;
    if (C <= 0) throw std::runtime_error("FixedTourCarAssignerDP: inst->cars() <= 0");
    if (C > 20) throw std::runtime_error("FixedTourCarAssignerDP: too many cars for bitmask DP (C > 20).");

    scratch.ensure(N, C);

    // 1) Materialize the tour once (including closure)
    for (int i = 0; i < N; ++i) scratch.tour[(std::size_t)i] = getNode(i);
    scratch.tour[(std::size_t)N] = scratch.tour[0];

    // 2) Validation without allocations (only if enabled)
    if (opt_.validateTour)
    {
        if (scratch.tour[0] != 0)
            throw std::runtime_error("FixedTourCarAssignerDP: node[0] must be 0 (fixed start).");

        int st = scratch.stamp++;
        if (scratch.stamp == std::numeric_limits<int>::max())
        {
            std::fill(scratch.seenStamp.begin(), scratch.seenStamp.end(), 0);
            scratch.stamp = 1;
        }

        for (int k = 0; k < N; ++k)
        {
            const int vv = scratch.tour[(std::size_t)k];
            if (vv < 0 || vv >= N)
                throw std::runtime_error("FixedTourCarAssignerDP: node value out of range.");
            if (scratch.seenStamp[(std::size_t)vv] == st)
                throw std::runtime_error("FixedTourCarAssignerDP: duplicate node in tour.");
            scratch.seenStamp[(std::size_t)vv] = st;
        }
    }
    else
    {
        // minimal check (cheap now since we have scratch.tour)
        if (N > 0 && scratch.tour[0] != 0)
            throw std::runtime_error("FixedTourCarAssignerDP: node[0] must be 0 (fixed start).");
    }

    // 3) DP uses a "fast view" (array access), no getNode calls inside the loops
    auto getNodeFast = [&](int i) -> int { return scratch.tour[(std::size_t)i]; };

    // fallback logic, but uses scratch on every attempt
    if (!opt_.fallbackToUnlimitedOnInf)
        return evaluateCostViewWithMaxSegLenScratch_(getNodeFast, opt_.maxSegmentLen, scratch);

    int L = opt_.maxSegmentLen;

    {
        const double c0 = evaluateCostViewWithMaxSegLenScratch_(getNodeFast, L, scratch);
        if (std::isfinite(c0) || L <= 0) return c0;
    }

    L = std::max(1, L);
    for (int r = 0; r < std::max(0, opt_.fallbackMaxRetries); ++r)
    {
        if (L >= N) break;
        const int nextL = std::min(N, L * std::max(2, opt_.fallbackGrowthFactor));
        if (nextL == L) break;
        L = nextL;

        const double cc = evaluateCostViewWithMaxSegLenScratch_(getNodeFast, L, scratch);
        if (std::isfinite(cc)) return cc;
    }

    return evaluateCostViewWithMaxSegLenScratch_(getNodeFast, -1, scratch);
}


double FixedTourCarAssignerDP::evaluateCostView(const std::function<int(int)>& getNode) const
{
    // uses the internal per-instance scratch_
    return evaluateCostViewScratch(getNode, scratch_);
}

double FixedTourCarAssignerDP::evaluateCost(const std::vector<int>& nodes) const
{
    const int N = inst_->n();
    if (opt_.validateTour)
    {
        validateTourOrThrow_(nodes);
    }
    else
    {
        if ((int)nodes.size() != N)
            throw std::runtime_error("FixedTourCarAssignerDP: nodes.size() != inst->n()");
    }

    if ((int)nodes.size() != N)
        throw std::runtime_error("FixedTourCarAssignerDP: nodes.size() != inst->n()");

    return evaluateCostView([&](int i) { return nodes[(std::size_t)i]; });
}


} // namespace localSearch

