#include "ivc_adapter.h"
#include "prove_leaf.h"
#include "proof_serialization.h"
#include "pol_verifier.h"
#include "GKR.h"
#include "proof_utils.h"
#include "polynomial.h"
#include "utils.hpp"
#include <iostream>
#include <vector>
#include <string> 
#include <sstream>
#include <cassert>
#include <stdexcept>
#include <algorithm>

using namespace std;

namespace {

proof PadProofForEncoding(const proof &original) {
  proof padded = original;

  if (padded.q_poly.empty()) {
    padded.q_poly.push_back(quadratic_poly(F(0), F(0), F(0)));
  }

  if (padded.randomness.empty()) {
    padded.randomness.push_back(vector<F>());
  }
  if (padded.randomness.size() < 2) {
    padded.randomness.resize(2);
  }

  if (padded.sig.empty()) {
    padded.sig.push_back(vector<F>{F(0)});
  }

  if (padded.final_claims_v.size() < padded.sig.size()) {
    padded.final_claims_v.resize(padded.sig.size());
  }
  for (size_t i = 0; i < padded.sig.size(); ++i) {
    if (padded.sig[i].empty()) {
      padded.sig[i].push_back(F(0));
    }
    if (padded.final_claims_v[i].size() < padded.sig[i].size()) {
      padded.final_claims_v[i].resize(padded.sig[i].size(), F(0));
    }
  }

  if (padded.vr.size() < padded.sig.size()) {
    padded.vr.resize(padded.sig.size(), F(0));
  }
  if (padded.gr.size() < padded.sig.size()) {
    padded.gr.resize(padded.sig.size(), F(0));
  }

  return padded;
}

}

AggregationResult FA_Aggregate(const std::vector<LeafResult>& children,
                               const std::string& prev_acc) {
  AggregationResult result;
  result.serialized = prev_acc;
  result.proof_struct = proof{};
  result.ok = false;
  if (children.empty()) {
    return result; // Return previous accumulator if no children
  }
  
  std::cout << "[FA_Aggregate] folding " << children.size()
            << " leaves with prev_acc=" << prev_acc << "\n";
  
  // Step 1: Collect all proofs from children (make mutable copies)
  vector<struct proof> proofs;
  for (size_t i = 0; i < children.size(); ++i) {
    if (!children[i].transcript.empty()) {
      proofs.insert(proofs.end(), children[i].transcript.begin(), children[i].transcript.end());
      std::cout << "  child" << i << ": transcript proofs=" << children[i].transcript.size() << "\n";
    } else if (children[i].step_proof.type != 0) {
      proofs.push_back(children[i].step_proof);
      std::cout << "  child" << i << ": proof_type=" << children[i].step_proof.type << " (fallback)\n";
    } else {
      std::cout << "  child" << i << ": no transcript available" << "\n";
    }

    if (!children[i].logup_proofs.empty()) {
      proofs.insert(proofs.end(), children[i].logup_proofs.begin(), children[i].logup_proofs.end());
      std::cout << "  child" << i << ": logup proofs=" << children[i].logup_proofs.size() << "\n";
    }
  }

  if (proofs.empty()) {
    std::cerr << "[FA_Aggregate] No proofs collected; returning previous accumulator" << std::endl;
    return result;
  }
  
  // Ensure proofs are properly initialized before encoding
  for (size_t i = 0; i < proofs.size(); ++i) {
    if (proofs[i].type == GKR_PROOF) {
      // Ensure all required vectors are initialized
      if (proofs[i].randomness.empty()) {
        proofs[i].randomness.push_back(vector<F>());
      }
      if (proofs[i].randomness.size() < 2) {
        while (proofs[i].randomness.size() < 2) {
          proofs[i].randomness.push_back(vector<F>());
        }
      }
      if (proofs[i].gr.empty() && !proofs[i].vr.empty()) {
        proofs[i].gr.resize(proofs[i].vr.size(), F(0));
      }
      if (proofs[i].final_claims_v.empty() && !proofs[i].sig.empty()) {
        proofs[i].final_claims_v.resize(proofs[i].sig.size());
        for (size_t j = 0; j < proofs[i].sig.size(); ++j) {
          proofs[i].final_claims_v[j].resize(proofs[i].sig[j].size(), F(0));
        }
      }
    }
  }
  
  // Step 2: If there's a previous accumulator, deserialize it and add to proofs
  if (!prev_acc.empty() && prev_acc != "ACC_INIT") {
    try {
      proof prev_proof = DeserializeProofFromString(prev_acc);
      proofs.push_back(prev_proof);
      std::cout << "[FA_Aggregate] included previous accumulator proof\n";
    } catch (const std::exception &e) {
      std::cerr << "[FA_Aggregate] Failed to deserialize previous accumulator: "
                << e.what() << "\n";
    } catch (...) {
      std::cerr << "[FA_Aggregate] Unknown error deserializing previous accumulator\n";
    }
  }
  
  try {
    vector<proof> normalized;
    normalized.reserve(proofs.size());
    for (const auto &p : proofs) {
      normalized.push_back(PadProofForEncoding(p));
    }
    VerifyBundle bundle = BuildVerificationBundle(normalized);
    std::vector<F> gkr_data;
    for (const auto &chunk : bundle.data) {
      gkr_data.insert(gkr_data.end(), chunk.begin(), chunk.end());
    }
    gkr_data.push_back(F(1));

    std::vector<F> randomness(10);
    for (auto &r : randomness) {
      r = random();
    }

    struct proof acc_proof = prove_verification(gkr_data, randomness, bundle.transcript);
    std::string accumulator = SerializeProofToString(acc_proof);
    std::cout << " => new_acc (size=" << accumulator.size() << ")\n";
    result.serialized = accumulator;
    result.proof_struct = acc_proof;
    result.ok = true;
    return result;
  } catch (const std::exception &e) {
    std::cerr << "[FA_Aggregate] Error generating accumulator: " << e.what() << "\n";
  } catch (...) {
    std::cerr << "[FA_Aggregate] Unknown error generating accumulator\n";
  }

  std::cerr << "[FA_Aggregate] Falling back to previous accumulator" << std::endl;
  return result;
}
