#include "witness_snapshot.h"
#include "aux_witness.h"
#include "aux_sumcheck.h"
#include "activation_gkr.h"
#include "proof_utils.h"
#include "mimc.h"
#include "config_pc.hpp"
#include "quantization.h"
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <stdexcept>

extern void init_hash();
extern void init_SHA();

static void init_env() {
  init_hash();
  init_SHA();
  current_randomness = F_ZERO;
  x_transcript.clear();
  y_transcript.clear();
}

static int RequiredBits(const FieldVector &vec) {
  int max_bits = 1;
  for (const auto &elem : vec) {
    __int128 raw = elem.toint128();
    unsigned __int128 v = raw < 0
                              ? static_cast<unsigned __int128>(-raw)
                              : static_cast<unsigned __int128>(raw);
    int bits = 0;
    while (v > 0) {
      ++bits;
      v >>= 1;
    }
    if (bits > max_bits) max_bits = bits;
  }
  if (max_bits > 60) max_bits = 60; // field modulus is ~2^61, keep margin
  if (max_bits < 1) max_bits = 1;
  return max_bits;
}

static FieldVector MaskToWidthSigned(const FieldVector &vec, int bits) {
  const unsigned __int128 modulus =
      (bits >= 64) ? 0 : (static_cast<unsigned __int128>(1) << bits);
  const unsigned __int128 mask =
      (bits >= 64) ? (~static_cast<unsigned __int128>(0))
                   : (modulus - 1);

  FieldVector out(vec.size(), F_ZERO);
  for (size_t i = 0; i < vec.size(); ++i) {
    __int128 raw = vec[i].toint128();
    unsigned __int128 canonical = 0;
    if (raw >= 0) {
      canonical = static_cast<unsigned __int128>(raw) & mask;
    } else {
      unsigned __int128 abs = static_cast<unsigned __int128>(-raw) & mask;
      canonical = (modulus - abs) & mask;
      if (canonical == modulus) canonical = 0;
    }
    out[i] = F(static_cast<long long>(canonical));
  }
  return out;
}

static std::string ResolveSnapshotPath() {
  const char *env = std::getenv("WITNESS_SNAPSHOT_PATH");
  if (env && *env) {
    std::ifstream test(env);
    if (test.good())
      return std::string(env);
  }

  const char *candidates[] = {
      "../bench/witness_dump.json",
      "../../bench/witness_dump.json",
      "../../../bench/witness_dump.json",
      "bench/witness_dump.json"
  };

  for (const char *candidate : candidates) {
    std::ifstream test(candidate);
    if (test.good())
      return std::string(candidate);
  }
  throw std::runtime_error("Unable to locate witness_dump.json. Set WITNESS_SNAPSHOT_PATH.");
}

int main() {
  init_env();

  const std::string snapshot_path = ResolveSnapshotPath();
  WitnessSnapshot snapshot = LoadWitnessSnapshot(snapshot_path);
  assert(snapshot.hidden_size > 0);

  const int num_bits = Q;

  FieldVector masked_a = MaskToWidthSigned(snapshot.a_t0, num_bits);
  FieldVector masked_z = MaskToWidthSigned(snapshot.z_t0, num_bits);

  AuxWitness aux_a = BuildAuxWitness(masked_a, num_bits);
  AuxWitness aux_z = BuildAuxWitness(masked_z, num_bits);

  aux_a.scalars = masked_a;
  aux_z.scalars = masked_z;

  std::cerr << "[test] num_bits=" << num_bits << "\n";
  auto range_a = ProveAuxConsistency(aux_a, masked_a, num_bits);
  assert(range_a.type == RANGE_PROOF);
  auto range_z = ProveAuxConsistency(aux_z, masked_z, num_bits);
  assert(range_z.type == RANGE_PROOF);

  ActivationProofs tanh_proof = ProveActivationGKR(masked_a, aux_a,
                                                   snapshot.h_t0,
                                                   snapshot.exp_tanh,
                                                   ACTIVATION_TANH, num_bits);
  assert(tanh_proof.range.type == RANGE_PROOF);
  assert(tanh_proof.logup.type == LOGUP_PROOF);

  ActivationProofs softmax_proof = ProveActivationGKR(masked_z, aux_z,
                                                      snapshot.y_hat_t0,
                                                      snapshot.exp_softmax,
                                                      ACTIVATION_SOFTMAX, num_bits);
  assert(softmax_proof.range.type == RANGE_PROOF);
  assert(softmax_proof.logup.type == LOGUP_PROOF);

  std::puts("test_witness_snapshot OK");
  return 0;
}


