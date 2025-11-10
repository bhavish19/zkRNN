#pragma once

#include "types.h"
#include <vector>

using BitVec = std::vector<int>; // 0/1 per bit

struct AuxWitness {
  std::vector<BitVec> bits; // one BitVec per scalar in the input vector
  FieldVector flat_bits;    // flattened bits as field elements (LSB-first per scalar)
  FieldVector scalars;      // original quantised field elements corresponding to bits
  std::vector<FieldVector> bit_planes; // bit planes grouped by bit position
  int bit_width = 0;        // number of bits used for decomposition
};

AuxWitness BuildAuxWitness(const FieldVector &quantized, int Q);


