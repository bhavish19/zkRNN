#include <atomic>
#include <iostream>
#include <cmath>
#include <omp.h>
#include "poly_commit.h"
#include "utils.hpp"
#include "GKR.h"
#include "logup.hpp"
#include "bench.hpp"
#include "expanders.h"


double logup_commit_time = 0.0;
double logup_prove_time  = 0.0;
double logup_total_time  = 0.0;

extern bool __encode_initialized; // poly_commit 内部用

namespace
{
    struct TVecRootCache {
        std::mutex mu;
        std::unordered_map<uint64_t, std::vector<F>> roots; // key -> root_fields

        static inline uint64_t key(int logN, int level) {
            return (uint64_t(uint32_t(level)) << 32) | uint32_t(logN);
        }

        bool get(int logN, int level, std::vector<F>& out) {
            std::lock_guard<std::mutex> lk(mu);
            auto it = roots.find(key(logN, level));
            if (it == roots.end()) return false;
            out = it->second; // 复制很小（只有 1 或 8 个 F）
            return true;
        }

        void set(int logN, int level, const std::vector<F>& val) {
            std::lock_guard<std::mutex> lk(mu);
            roots[key(logN, level)] = val;
        }

        void clear() {
            std::lock_guard<std::mutex> lk(mu);
            roots.clear();
        }
    };

    static TVecRootCache g_t_root_cache;

static std::atomic<bool> g_orion_inited{false};
static std::atomic<bool> g_expand_inited{false};

#ifdef DEBUG_LOGUP
#define LOGUP_DBG(msg) do { std::cerr << "[logup] " << msg << std::endl; } while(0)
#else
#define LOGUP_DBG(msg) do {} while(0)
#endif

    static inline int next_pow2(int x)
    {
        int p = 1;
        while (p < x)
            p <<= 1;
        return p;
    }
    static inline int lg2i(int x)
    {
        int r = 0;
        while (x >>= 1)
            ++r;
        return r;
    }

    static inline long long to_index_from_F(const F &x, int logN)
    {
        const long long N = 1LL << logN;
        long long v = (long long)(x.toint128() % N);
        if (v < 0)
            v += N;
        return v;
    }

    static inline F random_F()
    {
        return F::random();
    }

    static inline F poly_eval(const F &x0, const F &x1, const F &x2, const F &x3, const F &u)
    {
        // (1/6) * [ (-x0)(u-1)(u-2)(u-3) + 3x1*u(u-2)(u-3) - 3x2*u(u-1)(u-3) + x3*u(u-1)(u-2) ]
        F c6_inv = F(6).inv();
        F term0 = (F_ZERO - x0) * (u - F(1)) * (u - F(2)) * (u - F(3));
        F term1 = F(3) * x1 * u * (u - F(2)) * (u - F(3));
        F term2 = F(3) * x2 * u * (u - F(1)) * (u - F(3));
        F term3 = x3 * u * (u - F(1)) * (u - F(2));
        return c6_inv * (term0 + term1 - term2 + term3);
    }

    // —— Sumcheck（度 1）:   sum f_i = S
    struct SC_Return
    {
        vector<F> random;
        F claim_f = F_ZERO;
        F claim_g = F_ZERO;
    };

    // sumcheck degree 1
    static SC_Return sumcheck_deg1(int l, std::vector<F> f, F S)
    {
        SC_Return s;
        s.random.resize(l);
        for (int i = l; i >= 1; --i)
        {
            F sum0 = F_ZERO, sum1 = F_ZERO;
            for (int j = 0; j < (1 << i); ++j)
                ((j & 1) ? sum1 : sum0) += f[j];
            F u = random_F();
            s.random[l - i] = u;
            S = u * (sum1 - sum0) + sum0;
            std::vector<F> nf(1 << (i - 1));
            #pragma omp parallel for
            for (int j = 0; j < (1 << (i - 1)); ++j)
                nf[j] = (F_ONE - u) * f[2 * j] + u * f[2 * j + 1];
            f.swap(nf);
        }
        s.claim_f = f[0];
        return s;
    }

    // sumcheck degree 3
    static SC_Return sumcheck_deg3(int l, const std::vector<F> &rbits, std::vector<F> f, std::vector<F> g, F S)
    {
        std::vector<F> lag(1 << l, F_ONE);
        for (int k = 0; k < l; ++k)
        {
            int B = 1 << k;
            for (int j = 0; j < (1 << l); j += (B << 1))
            {
                for (int t = 0; t < B; ++t)
                {
                    F a = lag[j + t];
                    lag[j + t] = a * (F_ONE - rbits[k]); // bit=0
                    lag[j + t + B] = a * rbits[k];       // bit=1
                }
            }
        }

        SC_Return s;
        s.random.resize(l);

        for (int i = l; i >= 1; --i)
        {
            int half = 1 << (i - 1);
            std::vector<F> S0(half, F_ZERO), S1(half, F_ZERO), S2(half, F_ZERO), S3(half, F_ZERO);

            #pragma omp parallel for
            for (int j = 0; j < (1 << i); j += 2)
            {
                F l0 = lag[j], l1 = lag[j + 1];
                F f0 = f[j], f1 = f[j + 1];
                F g0 = g[j], g1 = g[j + 1];

                S0[j >> 1] = l0 * f0 * g0;
                S1[j >> 1] = l1 * f1 * g1;
                S2[j >> 1] = (l1 + l1 - l0) * (f1 + f1 - f0) * (g1 + g1 - g0);
                S3[j >> 1] = (l1 + l1 + l1 - l0 - l0) * (f1 + f1 + f1 - f0 - f0) * (g1 + g1 + g1 - g0 - g0);
            }

            F sum0 = F_ZERO, sum1 = F_ZERO, sum2 = F_ZERO, sum3 = F_ZERO;
            for (int j = 0; j < half; ++j)
            {
                sum0 += S0[j];
                sum1 += S1[j];
                sum2 += S2[j];
                sum3 += S3[j];
            }

            F u = random_F();
            s.random[l - i] = u;
            S = poly_eval(sum0, sum1, sum2, sum3, u);

            std::vector<F> nlag(half), nf(half), ng(half);
#pragma omp parallel for
            for (int j = 0; j < half; ++j)
            {
                nlag[j] = lag[2 * j] + u * (lag[2 * j + 1] - lag[2 * j]);
                nf[j] = f[2 * j] + u * (f[2 * j + 1] - f[2 * j]);
                ng[j] = g[2 * j] + u * (g[2 * j + 1] - g[2 * j]);
            }
            lag.swap(nlag);
            f.swap(nf);
            g.swap(ng);
        }
        s.claim_f = f[0];
        s.claim_g = g[0];
        return s;
    }

    // 把 commitment “变成若干 F”放进 proof.sig（最小可用版：只导出根）
    // - 若是 MIMC：直接用 com.hashes_f 最后一层的根（F）
    // - 若是 SHA：导出 8 个 32bit 整数到 F（需要 merkle_tree::get_int32）
    static std::vector<F> commitment_root_fields(const commitment &com)
    {
        std::vector<F> out;
        if (!com.hashes_f.empty())
        {
            out.push_back(com.hashes_f.back()[0]); // MiMC
            return out;
        }
        if (!com.hashes_sha.empty())
        {
            std::vector<unsigned int> ints(8, 0);
            __hhash_digest root = com.hashes_sha.back()[0];
            merkle_tree::get_int32(ints, root);
            out.reserve(8);
            for (int i = 0; i < 8; ++i)
                out.emplace_back(F((long long)ints[i]));
        }
        return out;
    }
}

    void logup_reset_round_timers()
    {
        logup_commit_time = 0.0;
        logup_prove_time = 0.0;
        logup_total_time = 0.0;
    }

    void logup_accumulate_round()
    {
        logup_acc_commit_time += logup_commit_time;
        logup_acc_prove_time += logup_prove_time;
        logup_acc_total_time += logup_total_time;
        
        total_prove += logup_prove_time;
        total_logup += logup_total_time;
        total_logup_prove += logup_prove_time;
    }

    void logup_reset_all_accumulators()
    {
        logup_acc_commit_time = 0.0;
        logup_acc_prove_time = 0.0;
        logup_acc_total_time = 0.0;
    }


void init_hyrax_logup_once(int threads)
{
    if (g_orion_inited.exchange(true))
        return;
    omp_set_num_threads(threads);
    __encode_initialized = false;
}

static void commit_and_track(const vector<F> &vec,
                             vector<vector<F>> &enc,
                             commitment &com,
                             int level)
{
    double before = g_commit_logup;
    vector<F> commit_vec(vec.begin(), vec.end());
    size_t n = commit_vec.size();
    if (n <= 1) {
        F pad = (n == 0) ? F_ZERO : commit_vec[0];
        commit_vec.resize(2, pad);
        n = commit_vec.size();
    }
    int max_level = (n <= 1) ? 0 : static_cast<int>(std::floor(std::log2(static_cast<double>(n))));
    int effective_level = std::min(level, max_level);
    LOGUP_DBG("commit_and_track: requested level=" << level
              << " vec_size=" << vec.size()
              << " padded_size=" << n
              << " effective_level=" << effective_level);
    poly_commit(commit_vec, enc, com, effective_level, CommitSrc::Logup);
    double after = g_commit_logup;
    logup_commit_time += (after - before);
}

struct proof prove_logup(const ExpBatch &batch, int table_logN, int threads)
{

    printf("Prove logup.\n");
    init_hyrax_logup_once(threads);
    const int N = 1 << table_logN;
    const int M = (int)batch.Xq.size();

    if (!g_expand_inited.load(std::memory_order_acquire)) {
        LOGUP_DBG("expander_init(" << N << ")");
        expander_init(N);
        g_expand_inited.store(true, std::memory_order_release);
    }

    struct proof Pr;
    Pr.type = LOGUP_PROOF;

    if (M == 0)
    {
        Pr.output = {F(table_logN), F(0), F(0)};
        Pr.final_eval = F_ZERO;
        Pr.eval_point.clear();
        return Pr;
    }

    logup_reset_round_timers();
    
    clock_t start = clock();
    vector<F> f_vec;
    f_vec.reserve(M);
    for (int i = 0; i < M; ++i)
    {
        f_vec.push_back(F(to_index_from_F(batch.Xq[i], table_logN)));
    }
    int m_pad = next_pow2(M);
    f_vec.resize(m_pad, F_ZERO);
    LOGUP_DBG("vectors prepared M=" << M << " m_pad=" << m_pad << " N=" << N);

    vector<F> t_vec(N);
    #pragma omp parallel for
    for (int j = 0; j < N; ++j)
        t_vec[j] = F(j);

    vector<F> c_vec(N, F_ZERO);
    #pragma omp parallel for
    for (int i = 0; i < M; ++i)
    {
        long long idx = (long long)f_vec[i].toint128();
        c_vec[(size_t)idx] += F_ONE;
    }
    logup_prove_time = double(clock() - start) / CLOCKS_PER_SEC;

    vector<F> t_root_fields;
    if (!g_t_root_cache.get(table_logN, /*level=*/2, t_root_fields)) {
        vector<vector<F>> enc_t;
        commitment com_t;
        LOGUP_DBG("commit_and_track(t_vec)");
        commit_and_track(t_vec, enc_t, com_t, /*level=*/2);
        t_root_fields = commitment_root_fields(com_t);
        g_t_root_cache.set(table_logN, /*level=*/2, t_root_fields);
        // enc_t/com_t 只为第一次承诺服务，后续不再需要
    }
    vector<vector<F>> enc_f, enc_c;
    commitment com_f, com_c;
    LOGUP_DBG("commit_and_track(f_vec)");
    commit_and_track(f_vec, enc_f, com_f, 2);
    LOGUP_DBG("commit_and_track(c_vec)");
    commit_and_track(c_vec, enc_c, com_c, 2);

    F r = random_F();
    Pr.eval_point = {r};
    Pr.final_r = {r};

    start = clock();
    vector<F> K(m_pad), G(m_pad);
    #pragma omp parallel for
    for (int i = 0; i < m_pad; ++i)
    {
        K[i] = r + f_vec[i];
        G[i] = F_ONE / K[i];
    }

    vector<F> Hp(N), H(N);
    for (int j = 0; j < N; ++j)
    {
        Hp[j] = r + t_vec[j];
        H[j] = c_vec[j] / Hp[j];
    }
    logup_prove_time += double(clock() - start) / CLOCKS_PER_SEC;

    vector<vector<F>> enc_G, enc_H;
    commitment com_G, com_H;
    
    LOGUP_DBG("commit_and_track(G)");
    commit_and_track(G, enc_G, com_G, 2);
    LOGUP_DBG("commit_and_track(H)");
    commit_and_track(H, enc_H, com_H, 2);

    start = clock();
    F S = F_ZERO;
    for (int i = 0; i < m_pad; ++i)
        S += G[i];
    vector<F> rp1(lg2i(N)), rp2(lg2i(m_pad));
    for (auto &x : rp1)
        x = random_F();
    for (auto &x : rp2)
        x = random_F();
    Pr.sc_challenges = {rp1, rp2};
    F c_eval = evaluate_vector(c_vec, rp1);

    SC_Return ret1 = sumcheck_deg3(lg2i(N), rp1, H, Hp, c_eval);   // Σ eq(rp1,j) H_j Hp_j = c(rp1)
    SC_Return ret2 = sumcheck_deg3(lg2i(m_pad), rp2, G, K, F_ONE); // Σ eq(rp2,i) G_i K_i = 1
    SC_Return ret3 = sumcheck_deg1(lg2i(m_pad), G, S);             // Σ_i G_i = S
    SC_Return ret4 = sumcheck_deg1(lg2i(N), H, S);                 // Σ_j H_j = S
    logup_prove_time += double(clock() - start) / CLOCKS_PER_SEC;

    Pr.output = {F(table_logN), F(M), F(m_pad)};
    Pr.final_eval = S;

    Pr.sig.clear();
    Pr.sig.push_back(commitment_root_fields(com_f));
    Pr.sig.push_back(t_root_fields);
    Pr.sig.push_back(commitment_root_fields(com_c));
    Pr.sig.push_back(commitment_root_fields(com_G));
    Pr.sig.push_back(commitment_root_fields(com_H));

    Pr.final_claims_v = {
        {ret1.claim_f, ret1.claim_g, c_eval},
        {ret2.claim_f, ret2.claim_g},
        {ret3.claim_f},
        {ret4.claim_f}};

    logup_total_time = logup_commit_time + logup_prove_time;
    logup_accumulate_round();

    return Pr;
}