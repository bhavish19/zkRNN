#pragma once
#include "types.h"
#include "GKR.h"
#include <vector>
#include <tuple>

struct SumMatMulResult {
  struct proof proof;  // MATMUL_PROOF type
  F eval_point;        // Evaluation point for verification
  F claimed_sum;       // Claimed sum value
};

// Proves combined matrix-vector products: M_combined * x_combined = y_combined
// SUMMER strategy: combines multiple matrix-vector products into single proof
SumMatMulResult ProveSumMatMul(const std::vector<std::vector<FieldVector>> &matrices,
                               const std::vector<FieldVector> &vectors,
                               const FieldVector &result);
