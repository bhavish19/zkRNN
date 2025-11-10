#pragma once

#include <string>
#include <vector>
#include "types.h"    // your FieldVector typedef, etc.
#include "GKR.h"      // struct proof
#include "RNN.h"      // ExpBatch definition

struct LeafResult {
  struct proof step_proof;  // legacy placeholder (compressed proof if available)
  std::vector<struct proof> transcript;        // matmul/range proofs for the timestep
  std::vector<struct proof> logup_proofs;      // activation LOGUP proofs retained separately
};

struct TimeStepInput {
  FieldVector x_t;    // input vector in field representation
  FieldVector h_prev; // previous hidden state
};

struct RNNWeights {
  std::vector<std::vector<F>> W_x;  // [m][n] input-to-hidden weights
  std::vector<std::vector<F>> W_h;  // [m][m] hidden-to-hidden weights
  std::vector<std::vector<F>> W_y;  // [k][m] hidden-to-output weights
  FieldVector b1;                    // [m] hidden bias
  FieldVector b2;                    // [k] output bias
};

struct TimeStepOutput {
  FieldVector h_t;    // hidden state output
  FieldVector a_t;    // pre-activation (before tanh)
  FieldVector z_t;    // pre-softmax (before softmax)
  FieldVector yHat_t; // output (after softmax)
  ExpBatch exp_tanh;      // lookup witness for tanh (two_x, e2x)
  ExpBatch exp_softmax;   // lookup witness for softmax (z_shift, exp_shift)
};

struct ProveParams {
  int quantization_bits = 32;    // quantization bits (default) - renamed from Q to avoid macro conflict
  int pc_type = 1; // which PCS, default Orion
};

// Prove a single RNN timestep
// Combines linear operations (matrix-vector products) using ProveSumMatMul
// and non-linear operations (activations) using ProveActivationGKR
LeafResult ProveLeaf(const ProveParams &p,
                     const TimeStepInput &in,
                     const RNNWeights &weights,
                     const TimeStepOutput &out);