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
  cout << "[ProveLeaf] Entered function\n";
  cout.flush();
  
  // Verify dimensions
  int m = weights.W_h.size();  // hidden size
  int n = in.x_t.size();       // input size
  int k = weights.W_y.size();  // output size
  
  cout << "[ProveLeaf] Dimensions: m=" << m << ", n=" << n << ", k=" << k << "\n";
  cout.flush();
  
  // Check for empty matrices/vectors before accessing
  if (m == 0 || n == 0 || k == 0 || 
      weights.W_x.size() == 0 || weights.W_x[0].size() == 0 ||
      weights.W_h.size() == 0 || weights.W_h[0].size() == 0 ||
      weights.W_y.size() == 0 || weights.W_y[0].size() == 0) {
    cerr << "[ProveLeaf] Error: Empty matrices or vectors detected\n";
    exit(1);
  }
  
  cout << "[ProveLeaf] Passed empty check\n";
  cout.flush();
  
  cout << "[ProveLeaf] Running assertions...\n";
  cout.flush();
  
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
  
  cout << "[ProveLeaf] All assertions passed\n";
  cout.flush();
  
  // Step 1: Prove linear operations using ProveSumMatMul (SUMMER strategy)
  // Linear operations in a timestep:
  // 1. Wh = W_h * h_prev
  // 2. Wx = W_x * x_t
  // 3. a = Wh + Wx + b1
  // 4. Zy = W_y * h_t
  // 5. z = Zy + b2
  
  // For now, we'll prove each matrix-vector product separately
  // In a full implementation, we could combine them into a single canonical form
  
  cout << "[ProveLeaf] Starting linear proofs...\n";
  
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
  
  // Prove Wh = W_h * h_prev
  vector<vector<FieldVector>> matrices_wh;
  vector<FieldVector> vectors_wh;
  matrices_wh.push_back(weights.W_h);
  vectors_wh.push_back(in.h_prev);
  cout << "[ProveLeaf] Calling ProveSumMatMul for Wh (matrix size=" << weights.W_h.size() 
       << "x" << (weights.W_h.empty() ? 0 : weights.W_h[0].size()) 
       << ", vector size=" << in.h_prev.size() << ")\n";
  SumMatMulResult wh_proof = ProveSumMatMul(matrices_wh, vectors_wh, wh_result);
  cout << "[ProveLeaf] ProveSumMatMul for Wh completed\n";
  
  // Prove Wx = W_x * x_t
  vector<vector<FieldVector>> matrices_wx;
  vector<FieldVector> vectors_wx;
  matrices_wx.push_back(weights.W_x);
  vectors_wx.push_back(in.x_t);
  cout << "[ProveLeaf] Calling ProveSumMatMul for Wx\n";
  SumMatMulResult wx_proof = ProveSumMatMul(matrices_wx, vectors_wx, wx_result);
  cout << "[ProveLeaf] ProveSumMatMul for Wx completed\n";
  
  // Prove Zy = W_y * h_t
  vector<vector<FieldVector>> matrices_zy;
  vector<FieldVector> vectors_zy;
  matrices_zy.push_back(weights.W_y);
  vectors_zy.push_back(out.h_t);
  cout << "[ProveLeaf] Calling ProveSumMatMul for Zy\n";
  SumMatMulResult zy_proof = ProveSumMatMul(matrices_zy, vectors_zy, zy_result);
  cout << "[ProveLeaf] ProveSumMatMul for Zy completed\n";
  
  // Step 2: Prove non-linear operations using ProveActivationGKR
  // 1. h = tanh(a) - tanh activation
  // 2. yHat = softmax(z) - softmax activation
  
  cout << "[ProveLeaf] Starting activation proofs...\n";
  cout.flush();
  
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
  
  cout << "[ProveLeaf] Building auxiliary witness for a_t (size=" << out.a_t.size() << ")\n";
  cout.flush();
  
  // Build auxiliary witness for pre-activation values
  FieldVector masked_a = mask_signed(out.a_t);
  AuxWitness aux_a = BuildAuxWitness(masked_a, num_bits);
  
  cout << "[ProveLeaf] BuildAuxWitness for a_t completed\n";
  cout.flush();
  
  cout << "[ProveLeaf] Calling ProveActivationGKR for tanh\n";
  cout.flush();
  
  // Prove tanh activation: h = tanh(a)
  ActivationProofs tanh_proofs = ProveActivationGKR(masked_a, aux_a, out.h_t,
                                                    out.exp_tanh,
                                                    ACTIVATION_TANH, num_bits);
  
  cout << "[ProveLeaf] ProveActivationGKR for tanh completed\n";
  cout.flush();
  
  // Build auxiliary witness for pre-softmax values
  cout << "[ProveLeaf] Building auxiliary witness for z_t (size=" << out.z_t.size() << ")\n";
  cout.flush();
  
  FieldVector masked_z = mask_signed(out.z_t);
  AuxWitness aux_z = BuildAuxWitness(masked_z, num_bits);
  
  cout << "[ProveLeaf] BuildAuxWitness for z_t completed\n";
  cout.flush();
  
  cout << "[ProveLeaf] Calling ProveActivationGKR for softmax\n";
  cout.flush();
  
  // Prove softmax activation: yHat = softmax(z)
  ActivationProofs softmax_proofs = ProveActivationGKR(masked_z, aux_z, out.yHat_t,
                                                       out.exp_softmax,
                                                       ACTIVATION_SOFTMAX, num_bits);
  
  cout << "[ProveLeaf] ProveActivationGKR for softmax completed\n";
  cout.flush();
  
  // Step 3: Combine all proofs into a single leaf proof
  // For now, we'll create a combined proof structure
  // In a full implementation, we would merge all proof fields
  
  cout << "[ProveLeaf] Starting proof combination...\n";
  cout.flush();
  
  LeafResult result;
  result.transcript.push_back(wh_proof.proof);
  result.transcript.push_back(wx_proof.proof);
  result.transcript.push_back(zy_proof.proof);
  result.transcript.push_back(tanh_proofs.range);
  result.transcript.push_back(softmax_proofs.range);

  result.logup_proofs.push_back(tanh_proofs.logup);
  result.logup_proofs.push_back(softmax_proofs.logup);

  if (!result.transcript.empty()) {
    result.step_proof = result.transcript.front();
  }
  result.commitments = { "W_COMMIT", "H_COMMIT" }; // Placeholder
  result.evals.push_back(std::make_tuple("SIG", "R", "V")); // Placeholder
  result.accumulator_out = "ACC_PLACEHOLDER"; // Placeholder

  std::cout << "[ProveLeaf] Built transcript with " << result.transcript.size()
            << " arithmetic proofs and " << result.logup_proofs.size()
            << " logup proofs (Q=" << num_bits << ")\n";
  
  return result;
}
