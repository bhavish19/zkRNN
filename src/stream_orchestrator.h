#pragma once

#include "prove_leaf.h"
#include "ivc_adapter.h"
#include "finalise_proof.h"
#include "verifier_stub.h"
#include <vector>
#include <string>

class StreamOrchestrator {
public:
  StreamOrchestrator(const ProveParams &params, int arity_k);
  void OnInput(const TimeStepInput &in);
  std::string FinaliseAndOpen();

private:
  ProveParams params_;
  int arity_k_;
  std::vector<LeafResult> buffer_;
  std::string acc_;
  struct proof acc_proof_;
  bool acc_valid_ = false;
  RNNWeights weights_;
  bool weights_initialized_ = false;
  int input_dim_ = 0;
  int hidden_dim_ = 0;
  int output_dim_ = 0;

  void EnsureWeightsInitialized(int input_dim, int hidden_dim, int output_dim);
};