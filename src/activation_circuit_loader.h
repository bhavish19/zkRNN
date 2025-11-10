#pragma once

#include <cstddef>
#include <string>

#include "activation_types.h"

struct ActivationCircuitInfo {
  std::string circuit_path;   // Path to the circuit artefact
  std::string input_label;    // Human-friendly tag (e.g., "tanh")
  std::size_t lanes = 0;      // Number of activation lanes handled by the circuit
};

// Ensure that we have a KAIZEN circuit description for the requested activation and dimension.
// The loader selects the smallest available artefact that can accommodate the requested lanes.
ActivationCircuitInfo EnsureActivationCircuit(ActivationType type,
                                              std::size_t dimension);


