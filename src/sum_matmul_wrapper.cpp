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
    std::cout << "[ProveSumMatMul] Combining " << matrices.size()
              << " matrix-vector products into block-diagonal proof\n";

    // Multiple matrix-vector products: build a block-diagonal matrix so that
    // diag(M_0, ..., M_k) * concat(x_0, ..., x_k) = concat(M_0 x_0, ..., M_k x_k)
    size_t total_rows = 0;
    size_t total_cols = 0;
    vector<size_t> row_offsets;
    vector<size_t> col_offsets;
    row_offsets.reserve(matrices.size());
    col_offsets.reserve(matrices.size());

    for (size_t idx = 0; idx < matrices.size(); ++idx) {
      const auto &M = matrices[idx];
      assert(!M.empty());
      const size_t rows = M.size();
      const size_t cols = M[0].size();
      for (const auto &row : M) {
        assert(row.size() == cols);
      }
      assert(vectors[idx].size() == cols);

      row_offsets.push_back(total_rows);
      col_offsets.push_back(total_cols);
      total_rows += rows;
      total_cols += cols;
    }

    assert(result.size() == total_rows);

    vector<FieldVector> block(total_rows, FieldVector(total_cols, F_ZERO));
    FieldVector combined_vector(total_cols, F_ZERO);
    FieldVector computed(total_rows, F_ZERO);

    for (size_t idx = 0; idx < matrices.size(); ++idx) {
      const auto &M = matrices[idx];
      const auto &x = vectors[idx];
      const size_t row_off = row_offsets[idx];
      const size_t col_off = col_offsets[idx];
      for (size_t i = 0; i < M.size(); ++i) {
        F acc = F_ZERO;
        for (size_t j = 0; j < M[i].size(); ++j) {
          const F val = M[i][j];
          block[row_off + i][col_off + j] = val;
          acc = acc + val * x[j];
        }
        computed[row_off + i] = acc;
      }
      for (size_t j = 0; j < x.size(); ++j) {
        combined_vector[col_off + j] = x[j];
      }
    }

    for (size_t i = 0; i < total_rows; ++i) {
      if (computed[i] != result[i]) {
        printf("[ProveSumMatMul] Error: combined result mismatch at index %zu\n", i);
        exit(-1);
      }
    }

    return ProveSumMatMul({block}, {combined_vector}, result);
  }
}
