#include "activation_circuit_loader.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

bool FileHasContent(const std::string &path) {
  std::ifstream in(path);
  return in.good() && in.peek() != std::ifstream::traits_type::eof();
}

void ExpectCircuit(ActivationType type,
                   std::size_t requested,
                   std::size_t expected_dim) {
  ActivationCircuitInfo info = EnsureActivationCircuit(type, requested);
  assert(info.lanes == expected_dim);
  assert(FileHasContent(info.circuit_path));
  std::cout << "  located " << info.input_label
            << " circuit -> " << info.circuit_path
            << " (lanes=" << info.lanes << ")\n";
}

} // namespace

int main() {
  std::puts("Running test_activation_loader...");

  ExpectCircuit(ACTIVATION_TANH, 2, 2);
  ExpectCircuit(ACTIVATION_TANH, 5, 6);      // rounds up
  ExpectCircuit(ACTIVATION_SOFTMAX, 8, 12);  // rounds up

  bool threw = false;
  try {
    (void)EnsureActivationCircuit(ACTIVATION_SOFTMAX, 32);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  assert(threw && "Expected oversized activation request to throw");

  std::puts("test_activation_loader OK");
  return 0;
}

