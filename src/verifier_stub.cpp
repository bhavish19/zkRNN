#include "verifier_stub.h"
#include "proof_serialization.h"
#include "pol_verifier.h"
#include "proof_utils.h"
#include <iostream>
#include <stdexcept>

// Verification functions are declared in pol_verifier.h

bool VerifyProof(const FinalProof &proof) {
  // Reset MiMC transcript state to match proof generation environment
  x_transcript.clear();
  y_transcript.clear();
  current_randomness = F_ZERO;

  std::cout << "[VerifyProof] checking final proof (type=" << proof.final_proof.type << ")\n";
  std::cout << "  root_acc size=" << proof.root_acc.size() << "\n";
  
  // Step 1: Prefer the supplied proof structure; fall back to deserializing root_acc only
  // if we lack essential data (common for legacy or placeholder accumulators).
  struct proof p = proof.final_proof;
  auto HasEssentialData = [](const struct proof &candidate) {
    return candidate.type != 0 &&
           (!candidate.randomness.empty() || !candidate.q_poly.empty() ||
            !candidate.sig.empty());
  };

  bool used_supplied_struct = HasEssentialData(p);
  if (!used_supplied_struct &&
      !proof.root_acc.empty() &&
      proof.root_acc.find("ACC_FOLDED_") != 0 &&
      proof.root_acc.find("FINAL_PROOF_") != 0) {
    try {
      p = DeserializeProofFromString(proof.root_acc);
      used_supplied_struct = HasEssentialData(p);
      std::cout << "  Deserialized proof from root_acc (randomness="
                << (p.randomness.empty() ? 0 : static_cast<int>(p.randomness.size()))
                << ")\n";
    } catch (...) {
      std::cout << "  Could not deserialize root_acc, using proof structure directly\n";
    }
  }

  if (!used_supplied_struct) {
    std::cout << "  Warning: proof structure lacks essential data; cannot verify\n";
    return false;
  }
  
  // Step 2: Check that proof structure is valid
  if (p.type != GKR_PROOF && 
      p.type != MATMUL_PROOF && 
      p.type != RANGE_PROOF &&
      p.type != RANGE_PROOF_OPT &&
      p.type != RANGE_PROOF_LOOKUP) {
    std::cout << "  Invalid proof type: " << p.type << "\n";
    return false;
  }
  
  // Step 3: Perform cryptographic verification based on proof type
  try {
    if (p.type == GKR_PROOF) {
      std::cout << "  Verifying GKR proof cryptographically...\n";
      if (p.randomness.empty()) {
        std::cout << "  Error: GKR proof missing randomness transcript\n";
        return false;
      }
      if (p.q_poly.empty()) {
        std::cout << "  Error: GKR proof missing q_polys\n";
        return false;
      }
      if (p.sig.empty() || p.final_claims_v.empty()) {
        std::cout << "  Error: GKR proof missing sigma/final claims data\n";
        return false;
      }
      // verify_gkr mutates global MiMC transcript state; ensure vectors have expected shape
      if (p.gr.size() < p.vr.size()) {
        std::cout << "  Padding gr to match vr size (" << p.vr.size() << ")\n";
        auto gr_copy = p.gr;
        gr_copy.resize(p.vr.size(), F_ZERO);
        p.gr = gr_copy;
      }
      verify_gkr(p);
      std::cout << "  GKR proof verification: SUCCESS\n";
    } else if (p.type == MATMUL_PROOF) {
      // Verify matrix multiplication proof cryptographically
      std::cout << "  Verifying MATMUL proof cryptographically...\n";
      
      // Check that proof has required fields
      if (p.q_poly.empty()) {
        std::cout << "  Error: MATMUL proof has no q_poly\n";
        return false;
      }
      if (p.randomness.empty() || p.randomness[0].empty()) {
        std::cout << "  Error: MATMUL proof has invalid randomness\n";
        return false;
      }
      
      // Call cryptographic verification
      // Note: verify_matrix2matrix will exit(-1) if verification fails
      verify_matrix2matrix(p);
      std::cout << "  MATMUL proof verification: SUCCESS\n";
      
    } else if (p.type == RANGE_PROOF || 
               p.type == RANGE_PROOF_OPT || 
               p.type == RANGE_PROOF_LOOKUP) {
      // Verify bit decomposition proof cryptographically
      std::cout << "  Verifying RANGE proof cryptographically...\n";
      
      // Check that proof has required fields
      if (p.q_poly.empty() && p.c_poly.empty()) {
        std::cout << "  Error: RANGE proof has no polynomials\n";
        return false;
      }
      
      // Call cryptographic verification
      // Note: verify_bit_decomposition will exit(-1) if verification fails
      verify_bit_decomposition(p);
      std::cout << "  RANGE proof verification: SUCCESS\n";
    }
    
    std::cout << "  All cryptographic verifications passed\n";
    return true;
    
  } catch (const std::exception& e) {
    std::cout << "  Verification error: " << e.what() << "\n";
    return false;
  } catch (...) {
    // verify_* functions may exit(-1) on failure, which we can't catch
    // If we reach here, verification likely failed
    std::cout << "  Verification failed (caught exception)\n";
    return false;
  }
}
