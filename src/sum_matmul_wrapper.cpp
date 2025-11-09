#include "sum_matmul_wrapper.h"
#include "config_pc.hpp"
#include "utils.hpp"
#include "mimc.h"
#include "polynomial.h"
#include "quantization.h"
#include "proof_utils.h"
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

using namespace std;

// Helper function to combine matrices and vectors into canonical form
// For SUMMER strategy: combine multiple M*v into single M_combined * x_combined
static vector<vector<F>> combine_matrices(const vector<vector<FieldVector>> &matrices) {
  if (matrices.empty()) return {};
  
  // Get dimensions
  size_t total_cols = 0;
  size_t rows = matrices[0].size();
  
  // Calculate total columns (sum of all matrix columns)
  for (const auto &M : matrices) {
    if (!M.empty()) {
      total_cols += M[0].size();
      assert(M.size() == rows); // All matrices must have same number of rows
    }
  }
  
  // Combine matrices by concatenating columns
  vector<vector<F>> combined(rows, vector<F>(total_cols));
  size_t col_offset = 0;
  
  for (const auto &M : matrices) {
    for (size_t i = 0; i < M.size(); ++i) {
      for (size_t j = 0; j < M[i].size(); ++j) {
        combined[i][col_offset + j] = M[i][j];
      }
    }
    col_offset += (M.empty() ? 0 : M[0].size());
  }
  
  return combined;
}

// Helper function to combine vectors into single vector
static vector<F> combine_vectors(const vector<FieldVector> &vectors) {
  vector<F> combined;
  for (const auto &v : vectors) {
    combined.insert(combined.end(), v.begin(), v.end());
  }
  return combined;
}

SumMatMulResult ProveSumMatMul(const vector<vector<FieldVector>> &matrices,
                               const vector<FieldVector> &vectors,
                               const FieldVector &result) {
  assert(matrices.size() == vectors.size());
  assert(!matrices.empty());
  
  // For Phase 1: Handle single matrix-vector product case
  // We'll extend to multiple products in later phases
  if (matrices.size() == 1) {
    // Single matrix-vector product: M * x = y
    vector<vector<F>> M = matrices[0];
    vector<F> x = vectors[0];
    
    // Verify dimensions
    assert(M.size() == result.size());
    assert(M[0].size() == x.size());
    
    // Compute actual result to verify correctness
    vector<F> computed_result(M.size(), F_ZERO);
    for (size_t i = 0; i < M.size(); ++i) {
      for (size_t j = 0; j < M[i].size(); ++j) {
        computed_result[i] = computed_result[i] + M[i][j] * x[j];
      }
    }
    
    // Verify correctness
    for (size_t i = 0; i < result.size(); ++i) {
      if (computed_result[i] != result[i]) {
        printf("[ProveSumMatMul] Error: computed result[%zu] != expected\n", i);
        exit(-1);
      }
    }
    
    // Pad to power of 2 for sumcheck protocol (minimum size 2 so verifier logic applies)
    size_t padded_rows = 1;
    int log_rows = 0;
    while (padded_rows < M.size()) {
      padded_rows <<= 1;
      log_rows += 1;
    }
    if (padded_rows < 2) {
      padded_rows = 2;
      log_rows = 1;
    }

    size_t padded_cols = 1;
    int log_cols = 0;
    while (padded_cols < M[0].size()) {
      padded_cols <<= 1;
      log_cols += 1;
    }
    if (padded_cols < 2) {
      padded_cols = 2;
      log_cols = 1;
    }

    M.resize(padded_rows);
    for (auto &row : M) {
      row.resize(padded_cols, F_ZERO);
    }
    x.resize(padded_cols, F_ZERO);
    
    // Create diagonal matrix from x for matrix multiplication proof
    // M * x is equivalent to M * diag(x) where diag(x) is diagonal matrix
    vector<vector<F>> x_diag(padded_cols, vector<F>(padded_cols, F_ZERO));
    for (size_t i = 0; i < x.size(); ++i) {
      x_diag[i][i] = x[i];
    }
    
    // Generate randomness for evaluation
    vector<F> r1 = generate_randomness(log_rows, current_randomness);
    vector<F> r2 = generate_randomness(log_cols, current_randomness);
    vector<F> r_combined = r1;
    r_combined.insert(r_combined.end(), r2.begin(), r2.end());
    
    // Compute the actual matrix product M * x_diag (after padding)
    vector<vector<F>> M_x_diag(padded_rows, vector<F>(padded_cols, F_ZERO));
    for (size_t i = 0; i < padded_rows; ++i) {
      for (size_t j = 0; j < padded_cols; ++j) {
        for (size_t k = 0; k < padded_cols; ++k) {
          M_x_diag[i][j] = M_x_diag[i][j] + M[i][k] * x_diag[k][j];
        }
      }
    }
    
    // Flatten the result matrix and evaluate it at r_combined
    // This matches how prove_rnn_dWx computes the sum: evaluate_vector(flat_dw, r)
    vector<F> flat_M_x_diag = convert2vector(M_x_diag);
    F evaluated_sum = evaluate_vector(flat_M_x_diag, r_combined);
    
    // Note: _prove_matrix2matrix checks: previous_sum * F(1ULL<<Q) == (Pr.q_poly[0].eval(0) + Pr.q_poly[0].eval(1))
    // The sumcheck computes sum_i V[0][i] * V[1][i] without Q scaling, so we need to divide by Q
    // so that (sum_to_prove * Q) equals the evaluated_sum
    F sum_to_prove = evaluated_sum / F(1ULL<<Q);
    
    // Use _prove_matrix2matrix to generate the proof
    // This proves the sum of all elements in M * x_diag (after padding)
    struct proof Pr = _prove_matrix2matrix(M, x_diag, r_combined, sum_to_prove);
    #ifdef DEBUG_SUM_MATMUL
    std::cout << "[ProveSumMatMul] Proof type=" << Pr.type
              << " randomness[0]_size=" << (Pr.randomness.empty() ? 0 : Pr.randomness[0].size())
              << " randomness_size=" << Pr.randomness.size() << "\n";
    #endif
    
    SumMatMulResult result_obj;
    result_obj.proof = Pr;
    result_obj.eval_point = r_combined.empty() ? F_ZERO : r_combined.back();
    result_obj.claimed_sum = sum_to_prove;
    
    return result_obj;
  } else {
    // Multiple matrix-vector products case
    // For Phase 1, we'll combine them into single canonical form
    // TODO: Implement proper SUMMER-style combination in later phase
    
    // Combine matrices and vectors
    vector<vector<F>> M_combined = combine_matrices(matrices);
    vector<F> x_combined = combine_vectors(vectors);
    
    // Verify dimensions
    assert(M_combined.size() == result.size());
    assert(M_combined[0].size() == x_combined.size());
    
    // For now, treat as single matrix-vector product
    return ProveSumMatMul({M_combined}, {x_combined}, result);
  }
}
