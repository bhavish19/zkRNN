#pragma once
#include "aux_witness.h"
#include "GKR.h"
#include "types.h"
#include <vector>

// Prove auxiliary witness consistency - RANGE_PROOF for bit-decomposition
// This proves that the bits in aux correctly represent the quantized values in preact
struct proof ProveAuxConsistency(const AuxWitness &aux,
                                 const FieldVector &preact,
                                 int num_bits);