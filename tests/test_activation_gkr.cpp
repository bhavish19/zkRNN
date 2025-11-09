#include "activation_gkr.h"
#include "aux_witness.h"
#include "RNN.h"
#include "config_pc.hpp"
#include "proof_utils.h"
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
  FieldVector preact{F(1), F(2)};
  FieldVector h_out{F(1), F(2)};
  int num_bits = 4;
  AuxWitness aux = BuildAuxWitness(preact, num_bits);
  
  // Verify auxiliary witness structure
  assert(aux.bits.size() == preact.size());
  for (size_t i = 0; i < aux.bits.size(); ++i) {
    assert((int)aux.bits[i].size() == num_bits);
  }
  assert(aux.flat_bits.size() == preact.size() * static_cast<size_t>(num_bits));
  assert(aux.scalars.size() == preact.size());
  assert(aux.bit_planes.size() == static_cast<size_t>(num_bits));
  assert(aux.bit_width == num_bits);
  for (size_t i = 0; i < preact.size(); ++i) {
    assert(aux.scalars[i] == preact[i]);
    for (int j = 0; j < num_bits; ++j) {
      assert(aux.bit_planes[j][i] == F(aux.bits[i][j]));
    }
  }

  // Build simple ExpBatch witness (mirror inputs for test purposes)
  ExpBatch exp_batch;
  for (size_t i = 0; i < preact.size(); ++i) {
    exp_batch.emit(preact[i], h_out[i]);
  }

  // Generate activation proofs
  ActivationProofs proofs = ProveActivationGKR(preact, aux, h_out,
                                               exp_batch,
                                               ACTIVATION_TANH, num_bits);

  // Verify RANGE proof structure
  assert(proofs.range.type == RANGE_PROOF);
  assert(!proofs.range.q_poly.empty());
  assert(!proofs.range.c_poly.empty());
  assert(proofs.range.randomness.size() >= 2);

  // Verify LOGUP proof structure
  assert(proofs.logup.type == LOGUP_PROOF);
  assert(!proofs.logup.sig.empty());
  assert(!proofs.logup.final_claims_v.empty());
  assert(proofs.logup.randomness.size() >= 1);
  std::puts("test_activation_gkr OK");
  return 0;
}
