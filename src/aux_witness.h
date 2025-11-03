#pragma once

#include "types.h"
#include <vector>

using BitVec = std::vector<int>; // 0/1 per bit

struct AuxWitness {
  std::vector<BitVec> bits; // one BitVec per scalar in the input vector
};

AuxWitness BuildAuxWitness(const FieldVector &quantized, int Q);


