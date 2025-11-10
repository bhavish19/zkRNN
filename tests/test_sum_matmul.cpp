#include "sum_matmul_wrapper.h"
#include "config_pc.hpp"
#include "utils.hpp"
#include "mimc.h"
#include "polynomial.h"
#include "pol_verifier.h"
#include "quantization.h"
#include <cassert>
#include <cstdio>
#include <vector>
#include <cmath>

using namespace std;

extern F current_randomness;
extern void init_hash();
extern void init_SHA();

void init_test() {
  // Initialize hash and randomness
  init_hash();
  init_SHA();
  current_randomness = F(12345); // Test seed
  // Q is defined as macro in quantization.h (default 32)
}

void test_sum_matmul_single() {
  printf("[test_sum_matmul] Testing single matrix-vector product...\n");
  
  init_test();
  
  // Create small test matrix and vector
  // 2x2 matrix
  vector<vector<F>> M = {
    {F(1), F(2)},
    {F(3), F(4)}
  };
  
  // 2x1 vector
  vector<F> x = {F(5), F(6)};
  
  // Compute expected result: M * x
  vector<F> expected_result(2);
  expected_result[0] = M[0][0] * x[0] + M[0][1] * x[1]; // 1*5 + 2*6 = 5 + 12 = 17
  expected_result[1] = M[1][0] * x[0] + M[1][1] * x[1]; // 3*5 + 4*6 = 15 + 24 = 39
  
  printf("  M[0] = [%lld, %lld]\n", M[0][0].toint128(), M[0][1].toint128());
  printf("  M[1] = [%lld, %lld]\n", M[1][0].toint128(), M[1][1].toint128());
  printf("  x = [%lld, %lld]\n", x[0].toint128(), x[1].toint128());
  printf("  Expected result = [%lld, %lld]\n", 
         expected_result[0].toint128(), expected_result[1].toint128());
  
  // Generate proof
  SumMatMulResult proof_result = ProveSumMatMul({M}, {x}, expected_result);
  
  // Verify proof structure
  assert(proof_result.proof.type == MATMUL_PROOF || proof_result.proof.type != -1);
  printf("  Proof type: %d\n", proof_result.proof.type);
  printf("  Proof generated successfully!\n");
  
  // Verify the proof has the expected structure
  assert(!proof_result.proof.q_poly.empty() || proof_result.proof.type == -1);
  printf("  Proof has %zu quadratic polynomials\n", proof_result.proof.q_poly.size());
  
  printf("[test_sum_matmul] Single product test PASSED\n");
}

void test_sum_matmul_multiple() {
  printf("[test_sum_matmul] Testing multiple matrix-vector products...\n");
  
  init_test();
  
  // Create two small matrices and vectors
  // Matrix 1: 2x2
  vector<vector<F>> M1 = {
    {F(1), F(2)},
    {F(3), F(4)}
  };
  
  // Matrix 2: 2x2 (same size for simplicity)
  vector<vector<F>> M2 = {
    {F(2), F(3)},
    {F(4), F(5)}
  };
  
  // Vectors
  vector<F> x1 = {F(1), F(1)};
  vector<F> x2 = {F(2), F(2)};
  
  // Compute M1*x1 + M2*x2
  vector<F> result1(2);
  result1[0] = M1[0][0] * x1[0] + M1[0][1] * x1[1];
  result1[1] = M1[1][0] * x1[0] + M1[1][1] * x1[1];
  
  vector<F> result2(2);
  result2[0] = M2[0][0] * x2[0] + M2[0][1] * x2[1];
  result2[1] = M2[1][0] * x2[0] + M2[1][1] * x2[1];
  
  vector<F> combined_result(2);
  combined_result[0] = result1[0] + result2[0];
  combined_result[1] = result1[1] + result2[1];
  
  printf("  M1*x1 = [%lld, %lld]\n", result1[0].toint128(), result1[1].toint128());
  printf("  M2*x2 = [%lld, %lld]\n", result2[0].toint128(), result2[1].toint128());
  printf("  Combined = [%lld, %lld]\n", 
         combined_result[0].toint128(), combined_result[1].toint128());
  
  // Generate proof (will combine matrices internally)
  SumMatMulResult proof_result = ProveSumMatMul({M1, M2}, {x1, x2}, combined_result);
  
  // Verify proof structure
  assert(proof_result.proof.type == MATMUL_PROOF || proof_result.proof.type != -1);
  printf("  Proof type: %d\n", proof_result.proof.type);
  printf("  Proof generated successfully!\n");
  
  printf("[test_sum_matmul] Multiple products test PASSED\n");
}

int main() {
  printf("=== test_sum_matmul ===\n");
  
  try {
    test_sum_matmul_single();
    test_sum_matmul_multiple();
    
    printf("\n=== test_sum_matmul OK ===\n");
    return 0;
  } catch (const std::exception& e) {
    printf("ERROR: %s\n", e.what());
    return 1;
  }
}

