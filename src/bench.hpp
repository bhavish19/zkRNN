#pragma once
#include <cstdio>
#include <string>
#include <ctime>
#include <unistd.h>
#include <sys/stat.h>

enum class CommitSrc { General, Logup, Aggregation };

extern double g_commit_general;
extern double g_commit_logup;
extern double g_commit_aggr;


inline void commit_accumulate(CommitSrc src, double sec) {
    switch (src) {
      case CommitSrc::Logup:       g_commit_logup   += sec; break;
      case CommitSrc::Aggregation: g_commit_aggr    += sec; break;
      default:                     g_commit_general += sec; break;
    }
}

inline void reset_commit_timers() {
    g_commit_general = g_commit_logup = g_commit_aggr = 0.0;
}
struct BenchmarkRow {
    // 运行环境
    string timestamp;     // ISO8601
    string hostname;
    string git_rev;       // 可留空

    // 模型/系统配置
    int T, n, m, k;
    int PC_scheme, Commitment_hash, levels, arity;
    long long param_count;     // n*m + m*m + k*m + m + k

    // 核心指标（秒 / MiB）
    double data_prove_time;               // PoGD dataset proving time
    double avg_update, avg_backward, avg_forward, avg_logup, avg_logup_prove, avg_pogd;
    double agg_recursion_total, agg_recursion_amort;    // 递归（Sumcheck+Merklize）总/摊销
    double aggr_time_total, aggr_prove_total;           // 聚合阶段（aggregate + prove_aggregation）的粗分
    double aggr_time_avg, aggr_prove_avg;
    double commit_time_avg;               // (g_commit_general+g_commit_aggr+g_commit_logup)/arity^2
    double recursion_verifier_time;
    double verifier_time;                 // 最终校验时间
    double pogd_proof_size;           // 第一轮 PoGD 证明体积
    double final_proof_size;          // 最终递归后的整体证明体积

    double peak_memory;
};

inline bool file_exists(const std::string& path) {
    struct stat st{}; return ::stat(path.c_str(), &st) == 0;
}

inline std::string now_iso8601() {
    char buf[64];
    std::time_t t = std::time(nullptr);
    std::tm tm{};
    gmtime_r(&t, &tm);
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm);
    return buf;
}

inline std::string get_hostname() {
    char buf[256]; buf[0] = '\0';
    if (::gethostname(buf, sizeof(buf)) == 0) return std::string(buf);
    return "unknown-host";
}

// 追加一行到 CSV；立即 fflush + fsync 防止崩溃丢日志
inline void append_benchmark_csv(const std::string& path, const BenchmarkRow& r) {
    const bool need_header = !file_exists(path);
    FILE* fp = std::fopen(path.c_str(), "a");
    if (!fp) return;

    if (need_header) {
        std::fprintf(fp,
            "timestamp,hostname,git_rev,"
            "T,n,m,k,PC_scheme,Commitment_hash,levels,arity,param_count,"
            "data_prove_time,avg_update,avg_backward,avg_forward,avg_logup,avg_logup_prove,avg_pogd,"
            "agg_recursion_total,agg_recursion_amort,"
            "aggr_time_total,aggr_time_avg,"
            "aggr_prove_total,aggr_prove_avg,"
            "commit_time_avg,verifier_time,recursion_verifier_time,"
            "pogd_proof_size,final_proof_size, peak_memory\n");
    }

    std::fprintf(fp,
        "%s,%s,%s,"
        "%d,%d,%d,%d,%d,%d,%d,%d,%lld,"
        "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,"
        "%.6f,%.6f,"
        "%.6f,%.6f,"
        "%.6f,%.6f,"
        "%.6f,%.6f,%.6f,"
        "%.6f,%.6f,%.6f\n",
        r.timestamp.c_str(), r.hostname.c_str(), r.git_rev.c_str(),
        r.T, r.n, r.m, r.k, r.PC_scheme, r.Commitment_hash, r.levels, r.arity, r.param_count,
        r.data_prove_time, r.avg_update, r.avg_backward, r.avg_forward, r.avg_logup, r.avg_logup_prove,r.avg_pogd,
        r.agg_recursion_total, r.agg_recursion_amort, r.aggr_time_total, r.aggr_time_avg,
        r.aggr_prove_total, r.aggr_prove_avg,
        r.commit_time_avg, r.verifier_time, r.recursion_verifier_time,
        r.pogd_proof_size, r.final_proof_size, r.peak_memory
    );

    std::fflush(fp);
    ::fsync(::fileno(fp));   // ★关键：立刻落盘
    std::fclose(fp);
}