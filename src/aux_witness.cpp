#include "aux_witness.h"
#include <cassert>
#include <cstdint>

AuxWitness BuildAuxWitness(const FieldVector &quantized, int Q) {
  // Robustness: constrain Q to a sane range for downstream bit gadgets.
  // Extend as needed, but unsigned 128-bit covers plenty; we assert 1..64 per recommendation.
  assert(Q >= 1 && Q <= 64);
  AuxWitness a;
  a.bits.reserve(quantized.size());
  for (size_t i=0;i<quantized.size();++i) {
    BitVec b;
    // Extract via unsigned __int128 to avoid signedness issues before masking/shifting
    unsigned __int128 uv = static_cast<unsigned __int128>(quantized[i].toint128());
    unsigned __int128 mask = (Q >= 64) ? (~static_cast<unsigned __int128>(0)) :
                             ((static_cast<unsigned __int128>(1) << Q) - 1);
    unsigned __int128 val = uv & mask;
    for (int bit=0; bit<Q; ++bit) {
      b.push_back(static_cast<int>((val >> bit) & 1u));
    }
    a.bits.push_back(b);
  }
  return a;
}


