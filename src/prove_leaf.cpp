#include "prove_leaf.h"
#include "aux_witness.h"
#include "activation_gkr.h"
#include "sum_matmul_wrapper.h"
#include "proof_utils.h"
#include "quantization.h"
#include "utils.hpp"
#include "pol_verifier.h"
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <cstdlib>

using namespace std;

// Prove a single RNN timestep
// Combines linear operations (matrix-vector products) using ProveSumMatMul
// and non-linear operations (activations) using ProveActivationGKR
LeafResult ProveLeaf(const ProveParams &p,
                     const TimeStepInput &in,
                     const RNNWeights &weights,
                     const TimeStepOutput &out) {
  // Verify dimensions
  int m = weights.W_h.size();  // hidden size
  int n = in.x_t.size();       // input size
  int k = weights.W_y.size();  // output size
  
  // Check for empty matrices/vectors before accessing
  if (m == 0 || n == 0 || k == 0 || 
      weights.W_x.size() == 0 || weights.W_x[0].size() == 0 ||
      weights.W_h.size() == 0 || weights.W_h[0].size() == 0 ||
      weights.W_y.size() == 0 || weights.W_y[0].size() == 0) {
    cerr << "[ProveLeaf] Error: Empty matrices or vectors detected\n";
    exit(1);
  }

  assert((int)weights.W_x.size() == m);
  assert((int)weights.W_x[0].size() == n);
  assert((int)weights.W_h[0].size() == m);
  assert((int)weights.W_y[0].size() == m);
  assert((int)weights.b1.size() == m);
  assert((int)weights.b2.size() == k);
  assert((int)in.h_prev.size() == m);
  assert((int)out.h_t.size() == m);
  assert((int)out.a_t.size() == m);
  assert((int)out.z_t.size() == k);
  assert((int)out.yHat_t.size() == k);
  // Step 1: Prove linear operations using ProveSumMatMul (SUMMER strategy)
  // Linear operations in a timestep:
  // 1. Wh = W_h * h_prev
  // 2. Wx = W_x * x_t
  // 3. a = Wh + Wx + b1
  // 4. Zy = W_y * h_t
  // 5. z = Zy + b2
  
  // For now, we'll prove each matrix-vector product separately
  // In a full implementation, we could combine them into a single canonical form
  cout << "[ProveLeaf] Phase 1: proving linear relations (matmuls)\n";

  // Compute actual matrix-vector products to match what ProveSumMatMul expects
  // Compute Wh = W_h * h_prev
  FieldVector wh_result(m, F(0));
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      wh_result[i] = wh_result[i] + weights.W_h[i][j] * in.h_prev[j];
    }
  }
  
  // Compute Wx = W_x * x_t
  FieldVector wx_result(m, F(0));
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      wx_result[i] = wx_result[i] + weights.W_x[i][j] * in.x_t[j];
    }
  }
  
  // Compute Zy = W_y * h_t
  FieldVector zy_result(k, F(0));
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < m; ++j) {
      zy_result[i] = zy_result[i] + weights.W_y[i][j] * out.h_t[j];
    }
  }
  
  // Combine linear proofs using SUMMER-style block diagonal
  vector<vector<FieldVector>> matrices_combined{weights.W_h, weights.W_x, weights.W_y};
  vector<FieldVector> vectors_combined{in.h_prev, in.x_t, out.h_t};
  FieldVector combined_result;
  combined_result.reserve(wh_result.size() + wx_result.size() + zy_result.size());
  combined_result.insert(combined_result.end(), wh_result.begin(), wh_result.end());
  combined_result.insert(combined_result.end(), wx_result.begin(), wx_result.end());
  combined_result.insert(combined_result.end(), zy_result.begin(), zy_result.end());

  SumMatMulResult linear_proof = ProveSumMatMul(matrices_combined, vectors_combined, combined_result);
  
  // Step 2: Prove non-linear operations using ProveActivationGKR
  // 1. h = tanh(a) - tanh activation
  // 2. yHat = softmax(z) - softmax activation
  
  cout << "[ProveLeaf] Phase 2: proving activation ranges/logups\n";
  
  // Store quantization bits in local variable
  int num_bits = p.quantization_bits;
  if (num_bits <= 0) num_bits = Q;
  if (num_bits > Q) num_bits = Q;

  auto mask_signed = [&](const FieldVector &values) {
    FieldVector masked(values.size(), F_ZERO);
    const unsigned __int128 modulus =
        (num_bits >= 64) ? 0 : (static_cast<unsigned __int128>(1) << num_bits);
    const unsigned __int128 mask =
        (num_bits >= 64) ? (~static_cast<unsigned __int128>(0))
                         : (modulus - 1);
    for (size_t i = 0; i < values.size(); ++i) {
      __int128 raw = values[i].toint128();
      unsigned __int128 canonical = 0;
      if (raw >= 0) {
        canonical = static_cast<unsigned __int128>(raw) & mask;
      } else {
        unsigned __int128 abs = static_cast<unsigned __int128>(-raw) & mask;
        canonical = (modulus - abs) & mask;
        if (canonical == modulus) canonical = 0;
      }
      masked[i] = F(static_cast<long long>(canonical));
    }
    return masked;
  };
  
  // Build auxiliary witness for pre-activation values
  FieldVector masked_a = mask_signed(out.a_t);
  AuxWitness aux_a = BuildAuxWitness(masked_a, num_bits);
  
  // Prove tanh activation: h = tanh(a)
  ActivationProofs tanh_proofs = ProveActivationGKR(masked_a, aux_a, out.h_t,
                                                    out.exp_tanh,
                                                    ACTIVATION_TANH, num_bits);
  
  // Build auxiliary witness for pre-softmax values
  FieldVector masked_z = mask_signed(out.z_t);
  AuxWitness aux_z = BuildAuxWitness(masked_z, num_bits);
  
  // Prove softmax activation: yHat = softmax(z)
  ActivationProofs softmax_proofs = ProveActivationGKR(masked_z, aux_z, out.yHat_t,
                                                       out.exp_softmax,
                                                       ACTIVATION_SOFTMAX, num_bits);
  
  // Step 3: Combine all proofs into a single leaf proof
  // For now, we'll create a combined proof structure
  // In a full implementation, we would merge all proof fields
  
  LeafResult result;
  result.transcript.push_back(linear_proof.proof);
  result.transcript.push_back(tanh_proofs.range);
  result.transcript.push_back(softmax_proofs.range);

  result.logup_proofs.push_back(tanh_proofs.logup);
  result.logup_proofs.push_back(softmax_proofs.logup);

  std::cout << "[ProveLeaf] Phase summary: "
            << result.transcript.size() << " arithmetic proofs, "
            << result.logup_proofs.size() << " logup proofs\n";
  
  return result;
}
