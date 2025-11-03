#include "aux_witness.h"
#include "config_pc.hpp" // F
#include <cassert>
#include <cstdio>

int main() {
  // Q=4; 3 -> 0011b, -1 -> 1111b (masked to 4 bits)
  FieldVector v; v.emplace_back(F(3)); v.emplace_back(F(-1));
  int Q = 4;
  AuxWitness a = BuildAuxWitness(v, Q);

  assert(a.bits.size() == 2);
  assert((int)a.bits[0].size() == Q);
  assert((int)a.bits[1].size() == Q);

  // 3: bits little-endian [1,1,0,0]
  assert(a.bits[0][0] == 1);
  assert(a.bits[0][1] == 1);
  assert(a.bits[0][2] == 0);
  assert(a.bits[0][3] == 0);

  // -1 mod 2^4 -> 1111
  for (int i=0;i<Q;++i) assert(a.bits[1][i] == 1);

  std::puts("test_aux_witness OK");
  return 0;
}


