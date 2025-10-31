#pragma once
#include "RNN.h"
#include <time.h>

extern double logup_acc_commit_time;
extern double logup_acc_prove_time;
extern double logup_acc_total_time;

extern double total_logup;
extern double total_logup_prove;
extern double total_prove;

struct CachedT { 
    vector<vector<F>> enc; 
    commitment com; 
};

static unordered_map<int, CachedT> g_t_cache;



void logup_accumulate_round();       // 本轮结束后累计到全局
void logup_reset_all_accumulators(); // 训练/epoch 开始时清零累计

void init_hyrax_logup_once(int threads);

struct proof prove_logup(const ExpBatch& batch, int table_logN, int threads);

struct proof prove_logup_batch(const ExpBatch& batch, const vector<F>& table_vals, int table_logN, int threads);
