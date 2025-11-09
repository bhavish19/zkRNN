#include "activation_witness.h"

#include <stdexcept>

FieldVector BuildActivationWitnessVector(const FieldVector &preact,
                                         const FieldVector &h_out,
                                         const AuxWitness &aux,
                                         const ActivationCircuitInfo &info) {
  if (preact.size() != h_out.size()) {
    throw std::invalid_argument("Activation witness mismatch: preact and output sizes differ");
  }
  if (aux.scalars.size() != preact.size()) {
    throw std::invalid_argument("AuxWitness scalars do not align with activation lanes");
  }
  if (aux.bit_width <= 0) {
    throw std::invalid_argument("AuxWitness bit width not initialised");
  }
  if (info.lanes == 0) {
    throw std::invalid_argument("Activation circuit info missing lane count");
  }

  const std::size_t stride = static_cast<std::size_t>(aux.bit_width) + 2; // preact + output + bits
  FieldVector layout;
  layout.reserve(info.lanes * stride);

  for (std::size_t lane = 0; lane < info.lanes; ++lane) {
    const bool has_value = lane < preact.size();

    layout.push_back(has_value ? preact[lane] : F_ZERO);
    layout.push_back(has_value ? h_out[lane] : F_ZERO);

    for (int bit = 0; bit < aux.bit_width; ++bit) {
      F bit_val = F_ZERO;
      if (has_value) {
        const int b = aux.bits[lane][bit];
        bit_val = F(b);
      }
      layout.push_back(bit_val);
    }
  }

  assert(layout.size() == info.lanes * stride);

  return layout;
}


