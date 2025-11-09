#include "aux_witness.h"
#include <cassert>
#include <cstdint>

AuxWitness BuildAuxWitness(const FieldVector &quantized, int Q) {
  // Robustness: constrain Q to a sane range for downstream bit gadgets.
  // Extend as needed, but unsigned 128-bit covers plenty; we assert 1..64 per recommendation.
  assert(Q >= 1 && Q <= 64);
  AuxWitness a;
  a.scalars = quantized;
  a.bit_width = Q;
  a.bit_planes.assign(static_cast<size_t>(Q), FieldVector());
  for (auto &plane : a.bit_planes) {
    plane.reserve(quantized.size());
  }

  if (quantized.empty()) {
    return a;
  }

  a.flat_bits.reserve(quantized.size() * static_cast<size_t>(Q));
  a.bits.reserve(quantized.size());
  for (size_t idx = 0; idx < quantized.size(); ++idx) {
    BitVec lane;
    lane.reserve(Q);
    unsigned __int128 uv = static_cast<unsigned __int128>(quantized[idx].toint128());
    unsigned __int128 mask = (Q >= 64) ? (~static_cast<unsigned __int128>(0))
                                       : ((static_cast<unsigned __int128>(1) << Q) - 1);
    unsigned __int128 val = uv & mask;
    for (int bit = 0; bit < Q; ++bit) {
      int bit_val = static_cast<int>((val >> bit) & 1u);
      F bit_f(bit_val);
      lane.push_back(bit_val);
      a.flat_bits.push_back(bit_f);
      a.bit_planes[bit].push_back(bit_f);
    }
    a.bits.push_back(std::move(lane));
  }

  return a;
}


