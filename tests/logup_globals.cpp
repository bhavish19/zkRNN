#include "logup.hpp"
#include "bench.hpp"

double logup_acc_commit_time = 0.0;
double logup_acc_prove_time  = 0.0;
double logup_acc_total_time  = 0.0;

double total_logup       = 0.0;
double total_logup_prove = 0.0;
double total_prove       = 0.0;

double g_commit_general = 0.0;
double g_commit_logup   = 0.0;
double g_commit_aggr    = 0.0;

double aggregation_time = 0.0;

int PC_scheme        = 1;
int Commitment_hash  = 1;
int arity            = 2;

