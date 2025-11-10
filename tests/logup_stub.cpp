#include "GKR.h"
#include "logup.hpp"

#include <omp.h>

// Provide stubbed accumulators for tests

double logup_acc_commit_time = 0.0;
double logup_acc_prove_time  = 0.0;
double logup_acc_total_time  = 0.0;

double total_logup        = 0.0;
double total_logup_prove  = 0.0;
double total_prove        = 0.0;

double g_commit_general   = 0.0;
double g_commit_logup     = 0.0;
double g_commit_aggr      = 0.0;

bool __encode_initialized = false;

void logup_accumulate_round() {}

void logup_reset_all_accumulators() {
  logup_acc_commit_time = 0.0;
  logup_acc_prove_time  = 0.0;
  logup_acc_total_time  = 0.0;
}

void init_hyrax_logup_once(int threads) {
  (void)threads;
}

namespace merkle_tree {
void get_int32(std::vector<unsigned int> &, __hhash_digest) {}
}

proof prove_logup(const ExpBatch &batch, int table_logN, int threads) {
  (void)threads;
  proof Pr;
  Pr.type = LOGUP_PROOF;
  Pr.output = {F(table_logN), F(batch.Xq.size()), F(batch.Yq.size())};
  Pr.sig.resize(1, std::vector<F>{F_ZERO});
  Pr.final_claims_v.resize(1, std::vector<F>{F_ZERO});
  Pr.randomness.resize(1, std::vector<F>{F_ZERO});
  Pr.vr = {F_ZERO, F_ZERO};
  Pr.gr = {F_ZERO, F_ZERO};

  layer_proof layer;
  layer.type = LOGUP_PROOF;
  layer.randomness.push_back(std::vector<F>{F_ZERO});
  layer.c_poly.push_back(cubic_poly(F_ZERO, F_ZERO, F_ZERO, F_ZERO));
  layer.q_poly.push_back(quadratic_poly(F_ZERO, F_ZERO, F_ZERO));
  layer.vr = {F_ZERO, F_ZERO};
  Pr.proofs.push_back(layer);

  return Pr;
}

proof prove_logup_batch(const ExpBatch &batch,
                        const std::vector<F> &table_vals,
                        int table_logN,
                        int threads) {
  (void)table_vals;
  return prove_logup(batch, table_logN, threads);
}
