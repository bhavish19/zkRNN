#include "stream_orchestrator.h"
#include "config_pc.hpp" // for F
#include "proof_utils.h"
#include "quantization.h"
#include <cassert>
#include <cstdio>

extern void init_hash();
extern void init_SHA();

void init_test() {
  init_hash();
  init_SHA();
  current_randomness = F(0);
}

int main() {
  std::puts("Running test_stream_orchestrator...");
  
  init_test();

  ProveParams p;
  p.quantization_bits = Q + 2;  // allow headroom for range proof
  p.pc_type = 1;

  StreamOrchestrator orchestrator(p, /*arity_k=*/2);

  const F q0 = quantize(0.0f);
  const F q1 = quantize(0.01f);
  const F q2 = quantize(0.02f);
  const F q3 = quantize(0.03f);

  TimeStepInput in1;
  in1.x_t = {q1};
  in1.h_prev = {q0};

  TimeStepInput in2;
  in2.x_t = {q2};
  in2.h_prev = {q1};

  // Feed two inputs -> triggers aggregation (k=2)
  orchestrator.OnInput(in1);
  orchestrator.OnInput(in2);

  // Add one more input, then finalize to force final aggregation
  TimeStepInput in3;
  in3.x_t = {q3};
  in3.h_prev = {q2};
  orchestrator.OnInput(in3);

  std::string final_acc = orchestrator.FinaliseAndOpen();
  assert(!final_acc.empty());
  std::puts("test_stream_orchestrator OK");
  return 0;
}
