#pragma once
#include "finalise_proof.h"
#include <string>

// Phase 8: Real cryptographic verification
// Verifies the final proof cryptographically using the appropriate verification functions
// Returns true if verification succeeds, false otherwise
bool VerifyProof(const FinalProof &proof);
