#include "activation_gkr.h"
#include "aux_sumcheck.h"
#include "proof_utils.h"
#include "logup.hpp"
#include "RNN.h"
#include <vector>
#include <cassert>
#include <stdexcept>

using namespace std;

namespace {

int DefaultLogupTableSize(ActivationType type) {
  switch (type) {
  case ACTIVATION_TANH:
  case ACTIVATION_SOFTMAX:
    return 15; // 2^15 entries, as in KAIZEN
  default:
    return 15;
  }
}

int DefaultLogupThreads() {
  return 32; // matches KAIZEN benchmark defaults
}

} // namespace

ActivationProofs ProveActivationGKR(const FieldVector &preact,
                                    const AuxWitness &aux,
                                    const FieldVector &h_out,
                                    const ExpBatch &exp_batch,
                                    ActivationType act_type,
                                    int num_bits) {
  // Verify dimensions match
  assert(preact.size() == h_out.size());
  assert(aux.bits.size() == preact.size());
  assert(aux.bit_width == num_bits);

  if (exp_batch.Xq.size() != exp_batch.Yq.size()) {
    throw std::invalid_argument("ExpBatch X/Y size mismatch for activation proof");
  }
  if (exp_batch.Xq.empty()) {
    throw std::invalid_argument("ExpBatch is empty; cannot generate logup proof");
  }
  
  ActivationProofs result;

  // Step 1: Prove auxiliary witness consistency (RANGE_PROOF)
  result.range = ProveAuxConsistency(aux, preact, num_bits);
  assert(result.range.type == RANGE_PROOF);

  // Step 2: Prove the lookup table relationships via KAIZEN logup protocol
  int table_logN = DefaultLogupTableSize(act_type);
  int threads = DefaultLogupThreads();
  result.logup = prove_logup(exp_batch, table_logN, threads);
  assert(result.logup.type == LOGUP_PROOF);

  return result;
}