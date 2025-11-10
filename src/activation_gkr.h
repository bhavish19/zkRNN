#pragma once
#include "types.h"
#include "aux_witness.h"

#include "activation_types.h"
#include "GKR.h"
#include <string>

struct ExpBatch;

struct ActivationProofs {
  struct proof logup;   // LOGUP proof for exponent table evaluation
  struct proof range;   // RANGE proof for auxiliary bits consistency
};

// Prove activation function using GKR circuit
// Returns a GKR_PROOF that includes:
// 1. Auxiliary witness consistency proof (RANGE_PROOF)
// 2. GKR proof for activation function computation
ActivationProofs ProveActivationGKR(const FieldVector &preact,
                                    const AuxWitness &aux,
                                    const FieldVector &h_out,
                                    const ExpBatch &exp_batch,
                                    ActivationType act_type,
                                    int num_bits);