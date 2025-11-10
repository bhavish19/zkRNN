#pragma once

#include "types.h"
#include "aux_witness.h"
#include "activation_circuit_loader.h"

// Assemble the per-gate activation witness vector expected by the activation
// circuit. Layout is interleaved on a per-lane basis:
//   [preact_0, h_out_0, bits_0..., preact_1, h_out_1, bits_1..., ...]
// The returned vector length equals `info.lanes * (aux.bit_width + 2)`. Lanes
// beyond the timestep count are zero padded; overflow raises.
FieldVector BuildActivationWitnessVector(const FieldVector &preact,
                                         const FieldVector &h_out,
                                         const AuxWitness &aux,
                                         const ActivationCircuitInfo &info);


