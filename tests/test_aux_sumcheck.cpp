#include "aux_sumcheck.h"
#include "aux_witness.h"
#include "config_pc.hpp"
#include "proof_utils.h"
#include "quantization.h"
#include "mimc.h"
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
  
  // Test with small values
  FieldVector pre{F(1), F(2)};
  int num_bits = 4;
  AuxWitness aux = BuildAuxWitness(pre, num_bits);
  
  // Verify auxiliary witness structure
  assert(aux.bits.size() == pre.size());
  for (size_t i = 0; i < aux.bits.size(); ++i) {
    assert((int)aux.bits[i].size() == num_bits);
  }
  assert(aux.flat_bits.size() == pre.size() * static_cast<size_t>(num_bits));
  assert(aux.scalars.size() == pre.size());
  for (size_t i = 0; i < pre.size(); ++i) {
    assert(aux.scalars[i] == pre[i]);
  }
  
  // Generate proof
  struct proof Pr = ProveAuxConsistency(aux, pre, num_bits);
  
  // Verify proof structure
  assert(Pr.type == RANGE_PROOF);
  assert(!Pr.q_poly.empty());
  assert(!Pr.c_poly.empty());
  assert(Pr.randomness.size() >= 2);
  assert(Pr.vr.size() >= 2);
  
  printf("  Proof type: %d (RANGE_PROOF)\n", Pr.type);
  printf("  Quadratic polynomials: %zu\n", Pr.q_poly.size());
  printf("  Cubic polynomials: %zu\n", Pr.c_poly.size());
  printf("  Randomness vectors: %zu\n", Pr.randomness.size());
  printf("  Verification values: %zu\n", Pr.vr.size());
  
  std::puts("test_aux_sumcheck OK");
  return 0;
}
