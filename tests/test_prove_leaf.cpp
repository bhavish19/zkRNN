#include "prove_leaf.h"
#include "config_pc.hpp"
#include "proof_utils.h"
#include "quantization.h"
#include "mimc.h"
#include <cassert>
#include <cstdio>
#include <vector>

extern void init_hash();
extern void init_SHA();

void init_test() {
  init_hash();
  init_SHA();
  current_randomness = F(0);
}

int main() {
  init_test();
  
  // Test with small values
  ProveParams p;
  p.quantization_bits = 4;
  
  // Create small RNN: n=2, m=2, k=2
  TimeStepInput in;
  in.x_t = {F(1), F(2)};
  in.h_prev = {F(3), F(4)};
  
  RNNWeights weights;
  // W_x: [2][2] - identity matrix
  weights.W_x = {{F(1), F(0)}, {F(0), F(1)}};
  // W_h: [2][2] - identity matrix
  weights.W_h = {{F(1), F(0)}, {F(0), F(1)}};
  // W_y: [2][2] - identity matrix
  weights.W_y = {{F(1), F(0)}, {F(0), F(1)}};
  // b1: [2] - zeros
  weights.b1 = {F(0), F(0)};
  // b2: [2] - zeros
  weights.b2 = {F(0), F(0)};
  
  // Compute expected outputs (simplified)
  // Wh = W_h * h_prev = [3, 4]
  // Wx = W_x * x_t = [1, 2]
  // a = Wh + Wx + b1 = [4, 6]
  // h = tanh(a) ≈ [4, 6] (simplified, should be tanh)
  // Zy = W_y * h = [4, 6]
  // z = Zy + b2 = [4, 6]
  // yHat = softmax(z) ≈ [0.5, 0.5] (simplified)
  
  TimeStepOutput out;
  out.a_t = {F(4), F(6)};
  out.h_t = {F(4), F(6)};  // Simplified - should be tanh(a_t)
  out.z_t = {F(4), F(6)};
  out.yHat_t = {F(2), F(3)};  // Simplified - should be softmax(z_t)
  for (size_t i = 0; i < out.a_t.size(); ++i) {
    out.exp_tanh.emit(out.a_t[i], out.h_t[i]);
  }
  for (size_t i = 0; i < out.z_t.size(); ++i) {
    out.exp_softmax.emit(out.z_t[i], out.yHat_t[i]);
  }
  
  // Generate proof
  LeafResult result = ProveLeaf(p, in, weights, out);
  
  // Verify proof structure
  assert(result.transcript.size() == 5);
  assert(result.transcript[0].type == MATMUL_PROOF);
  assert(result.transcript[1].type == MATMUL_PROOF);
  assert(result.transcript[2].type == MATMUL_PROOF);
  assert(result.transcript[3].type == RANGE_PROOF);
  assert(result.transcript[4].type == RANGE_PROOF);

  assert(result.logup_proofs.size() == 2);
  assert(result.logup_proofs[0].type == LOGUP_PROOF);
  assert(result.logup_proofs[1].type == LOGUP_PROOF);
  
  printf("  Transcript proofs: %zu\n", result.transcript.size());
  printf("  Logup proofs: %zu\n", result.logup_proofs.size());
  
  std::puts("test_prove_leaf OK");
  return 0;
}
