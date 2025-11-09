#include "aux_sumcheck.h"
#include "fieldElement.hpp"
#include "aux_witness.h"
#include "proof_utils.h"
#include "utils.hpp"
#include "quantization.h"
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

using namespace std;

// Prove auxiliary witness consistency - RANGE_PROOF for bit-decomposition
// This proves that the bits in aux correctly represent the quantized values in preact
struct proof ProveAuxConsistency(const AuxWitness &aux,
                                 const FieldVector &preact,
                                 int num_bits) {
  // Verify dimensions match
  assert(aux.bits.size() == preact.size());
  assert(aux.scalars.size() == preact.size());
  assert(aux.bit_width == num_bits);
  assert(aux.bit_planes.size() == static_cast<size_t>(num_bits));
  for (size_t i = 0; i < preact.size(); ++i) {
    assert(aux.scalars[i] == preact[i]);
    assert(static_cast<int>(aux.bits[i].size()) == num_bits);
    for (int j = 0; j < num_bits; ++j) {
      assert(aux.bit_planes[j].size() == preact.size());
      assert(aux.bit_planes[j][i] == F(aux.bits[i][j]));
    }
  }
  
  // Flattened bit view already prepared by BuildAuxWitness
  vector<F> bits = aux.flat_bits;
  
  // Pad bits and numbers to power of 2 for sumcheck protocol
  // Ensure minimum size of 2 for sumcheck to work correctly
  pad_vector(bits);
  if (bits.size() < 2) bits.resize(2, F_ZERO);
  
  vector<F> numbers = aux.scalars;
  pad_vector(numbers);
  if (numbers.size() < 2) numbers.resize(2, F_ZERO);
  
  // Generate randomness for evaluation
  // r length = log2(numbers.size()) after padding
  // This matches how _prove_bit_decomposition is called in the original code:
  // r = generate_randomness((int)log2(numbers.size()), F(32));
  int log_numbers = (int)ceil(log2(numbers.size()));
  if (log_numbers == 0) log_numbers = 1;  // Ensure at least 1 for size-1 vectors
  vector<F> r = generate_randomness(log_numbers, current_randomness);
  
  // Compute the sum we're proving: evaluate the numbers vector at randomness r
  // This matches: evaluate_vector(numbers, r)
  F sum_to_prove = evaluate_vector(numbers, r);
  
  // Verify bit-decomposition correctness locally
  for (size_t i = 0; i < preact.size(); ++i) {
    F reconstructed = F_ZERO;
    F power = F(1);
    for (int j = 0; j < num_bits; ++j) {
      reconstructed = reconstructed + aux.bit_planes[j][i] * power;
      power = power * F(2);
    }
    assert(reconstructed == aux.scalars[i]);
  }
  
  // Generate the bit-decomposition proof
  // Note: r is used for both evaluating numbers and preparing the bit matrix
#ifdef DEBUG_AUX_SUMCHECK
  auto print_int128 = [](std::int64_t hi, std::uint64_t lo) {
    if (hi == 0) {
      std::cerr << static_cast<long long>(lo);
    } else {
      __int128 full = (static_cast<__int128>(hi) << 64) |
                      static_cast<__int128>(lo);
      long double v = static_cast<long double>(full);
      std::cerr << static_cast<long double>(v);
    }
  };

  std::cerr << "[ProveAuxConsistency] num_bits=" << num_bits
            << " vector_size=" << numbers.size()
            << " sum=";
  print_int128(sum_to_prove.real, sum_to_prove.img);
  std::cerr << "\n";
  for (size_t idx = 0; idx < numbers.size(); ++idx) {
    std::cerr << "  number[" << idx << "]=";
    print_int128(numbers[idx].real, numbers[idx].img);
    std::cerr << " bits:";
    for (size_t b = 0; b < aux.bits[idx].size(); ++b)
      std::cerr << aux.bits[idx][b];
    std::cerr << "\n";
  }
#endif
  struct proof Pr = _prove_bit_decomposition(bits, r, sum_to_prove, num_bits);
  
  return Pr;
}
