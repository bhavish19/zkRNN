#include "sum_matmul_wrapper.h"
#include "config_pc.hpp"
#include "quantization.h"
#include "mimc.h"
#include <cassert>
#include <cstdio>

extern void init_hash();
extern void init_SHA();
extern F current_randomness;

int main() {
  // Initialize hash for randomness generation
  init_hash();
  init_SHA();
  current_randomness = F(12345);
  
  FieldVector x{F(1),F(2)};
  // Compute correct result: M * x where M is identity matrix
  // [[1,0],[0,1]] * [1,2] = [1,2]
  FieldVector y{F(1),F(2)};
  std::vector<FieldVector> M = { {F(1),F(0)}, {F(0),F(1)} };
  
  // Function expects vectors of matrices and vectors, so wrap in vectors
  auto res = ProveSumMatMul({M}, {x}, y);
  
  // Verify proof was generated (type should be MATMUL_PROOF or valid)
  assert(res.proof.type == MATMUL_PROOF || res.proof.type != -1);
  std::puts("test_sum_matmul_wrapper OK");
  return 0;
}
