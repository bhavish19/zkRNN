#include "ivc_adapter.h"
#include "prove_leaf.h"
#include "GKR.h"
#include "polynomial.h"
#include "proof_utils.h"
#include "config_pc.hpp"
#include "sum_matmul_wrapper.h"
#include "aux_witness.h"
#include "aux_sumcheck.h"
#include <cassert>
#include <cstdio>

extern void init_hash();
extern void init_SHA();

void init_test() {
  init_hash();
  init_SHA();
  current_randomness = F(0);
}

int main() {
  init_test();
  
  std::puts("Running test_ivc_adapter...");

  // Create placeholder leaf results with proof structures
  LeafResult a, b;

  // Matmul proof from identity matrix
  std::vector<std::vector<FieldVector>> matrices;
  std::vector<FieldVector> vectors;
  matrices.push_back({{F(1), F(0)}, {F(0), F(1)}});
  vectors.push_back({F(1), F(2)});
  FieldVector result_vec{F(1), F(2)};
  SumMatMulResult matmul_res = ProveSumMatMul(matrices, vectors, result_vec);
  a.transcript.push_back(matmul_res.proof);
  a.evals.push_back(std::make_tuple("SIG_A", "R_A", "V_A"));

  // Range proof using simple numbers
  FieldVector numbers{F(1), F(2)};
  AuxWitness aux_b = BuildAuxWitness(numbers, 4);
  proof range_res = ProveAuxConsistency(aux_b, numbers, 4);
  b.transcript.push_back(range_res);
  b.evals.push_back(std::make_tuple("SIG_B", "R_B", "V_B"));

  std::vector<LeafResult> leaves{a, b};

  AggregationResult res = FA_Aggregate(leaves, "ACC_PREV");

  // The result should be a non-empty serialized proof string
  assert(!res.serialized.empty());
  // The serialized string should contain the proof type and data
  // Format: proof_type|transcript_size|transcript_data|data_size|data
  assert(res.serialized.find("|") != std::string::npos); // Should contain separators

  std::puts("test_ivc_adapter OK");
  return 0;
}
