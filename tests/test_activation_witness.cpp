#include "activation_witness.h"
#include "activation_circuit_loader.h"
#include "aux_witness.h"

#include <cassert>
#include <iostream>

namespace {

void ExpectLaneLayout(const FieldVector &layout,
                      std::size_t lanes,
                      int bit_width,
                      const FieldVector &preact,
                      const FieldVector &h_out,
                      const AuxWitness &aux) {
  const std::size_t stride = static_cast<std::size_t>(bit_width) + 2;
  assert(layout.size() == lanes * stride);

  for (std::size_t lane = 0; lane < lanes; ++lane) {
    const bool has_value = lane < preact.size();
    const F expected_preact = has_value ? preact[lane] : F_ZERO;
    const F expected_output = has_value ? h_out[lane] : F_ZERO;

    assert(layout[lane * stride + 0] == expected_preact);
    assert(layout[lane * stride + 1] == expected_output);

    for (int bit = 0; bit < bit_width; ++bit) {
      F expected_bit = F_ZERO;
      if (has_value) {
        expected_bit = F(aux.bits[lane][bit]);
      }
      assert(layout[lane * stride + 2 + bit] == expected_bit);
    }
  }
}

} // namespace

int main() {
  std::puts("Running test_activation_witness...");

  FieldVector preact{F(3), F(5)};
  FieldVector h_out{F(7), F(11)};
  const int bit_width = 4;

  AuxWitness aux = BuildAuxWitness(preact, bit_width);
  ActivationCircuitInfo info_exact = EnsureActivationCircuit(ACTIVATION_TANH, preact.size());

  FieldVector layout_exact = BuildActivationWitnessVector(preact, h_out, aux, info_exact);
  ExpectLaneLayout(layout_exact, info_exact.lanes, bit_width, preact, h_out, aux);

  // Request a larger circuit and ensure padding zeros are appended.
  ActivationCircuitInfo info_padded = EnsureActivationCircuit(ACTIVATION_TANH, 6);
  FieldVector layout_padded = BuildActivationWitnessVector(preact, h_out, aux, info_padded);
  ExpectLaneLayout(layout_padded, info_padded.lanes, bit_width, preact, h_out, aux);

  std::cout << "  stride = " << (bit_width + 2)
            << ", lanes_exact = " << info_exact.lanes
            << ", lanes_padded = " << info_padded.lanes << "\n";

  std::puts("test_activation_witness OK");
  return 0;
}

