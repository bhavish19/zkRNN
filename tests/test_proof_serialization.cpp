#include "proof_serialization.h"
#include "sum_matmul_wrapper.h"
#include "aux_sumcheck.h"
#include "aux_witness.h"
#include "activation_gkr.h"
#include "RNN.h"
#include "prove_leaf.h"
#include "proof_utils.h"
#include "config_pc.hpp"
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
  init_test();
  
  std::puts("Running test_proof_serialization...");
  
  // Test 1: Serialize/Deserialize MATMUL_PROOF
  {
    std::puts("  Test 1: MATMUL_PROOF serialization");
    
    // Create a simple matrix-vector product proof
    std::vector<std::vector<FieldVector>> matrices;
    std::vector<FieldVector> vectors;
    FieldVector result;
    
    // Simple 2x2 identity matrix
    matrices.push_back({{F(1), F(0)}, {F(0), F(1)}});
    vectors.push_back({F(1), F(2)});
    result = {F(1), F(2)}; // Identity * [1,2] = [1,2]
    
    SumMatMulResult matmul_result = ProveSumMatMul(matrices, vectors, result);
    assert(matmul_result.proof.type == MATMUL_PROOF);
    
    // Serialize
    SerializedProof serialized = SerializeProof(matmul_result.proof);
    assert(serialized.proof_type == MATMUL_PROOF);
    assert(!serialized.data.empty());
    assert(!serialized.transcript.empty());
    
    // Deserialize
    struct proof deserialized = DeserializeProof(serialized);
    assert(deserialized.type == MATMUL_PROOF);
    assert(deserialized.q_poly.size() == matmul_result.proof.q_poly.size());
    
    std::puts("    MATMUL_PROOF serialization OK");
  }
  
  // Test 2: Serialize/Deserialize RANGE_PROOF
  {
    std::puts("  Test 2: RANGE_PROOF serialization");
    
    // Create a simple bit-decomposition proof
    FieldVector numbers{F(1), F(2), F(3)};
    int num_bits = 4;
    AuxWitness aux = BuildAuxWitness(numbers, num_bits);
    
    struct proof range_proof = ProveAuxConsistency(aux, numbers, num_bits);
    assert(range_proof.type == RANGE_PROOF);
    
    // Serialize
    SerializedProof serialized = SerializeProof(range_proof);
    assert(serialized.proof_type == RANGE_PROOF);
    assert(!serialized.data.empty());
    assert(!serialized.transcript.empty());
    
    // Deserialize
    struct proof deserialized = DeserializeProof(serialized);
    assert(deserialized.type == RANGE_PROOF);
    assert(deserialized.q_poly.size() == range_proof.q_poly.size());
    assert(deserialized.c_poly.size() == range_proof.c_poly.size());
    
    std::puts("    RANGE_PROOF serialization OK");
  }
  
  // Test 3: Serialize/Deserialize GKR_PROOF (from activation)
  {
    std::puts("  Test 3: GKR_PROOF serialization");
    
    FieldVector preact{F(1), F(2)};
    FieldVector h_out{F(1), F(2)};
    int num_bits = 4;
    AuxWitness aux = BuildAuxWitness(preact, num_bits);
    
    ExpBatch exp_batch;
    for (size_t i = 0; i < preact.size(); ++i) {
      exp_batch.emit(preact[i], h_out[i]);
    }

    ActivationProofs activation = ProveActivationGKR(preact, aux, h_out,
                                                     exp_batch,
                                                     ACTIVATION_TANH, num_bits);
    assert(activation.logup.type == LOGUP_PROOF);
    assert(activation.range.type == RANGE_PROOF);
 
    // Serialize
    SerializedProof serialized_logup = SerializeProof(activation.logup);
    SerializedProof serialized_range = SerializeProof(activation.range);
    assert(serialized_logup.proof_type == LOGUP_PROOF);
    assert(serialized_range.proof_type == RANGE_PROOF);
    assert(!serialized_logup.data.empty());
    assert(!serialized_logup.transcript.empty());
    assert(!serialized_range.data.empty());
    assert(!serialized_range.transcript.empty());
 
    // Deserialize
    struct proof deserialized_logup = DeserializeProof(serialized_logup);
    struct proof deserialized_range = DeserializeProof(serialized_range);
    assert(deserialized_logup.type == LOGUP_PROOF);
    assert(deserialized_range.type == RANGE_PROOF);
 
    std::puts("    GKR_PROOF serialization OK");
  }
  
  // Test 4: String serialization/deserialization
  {
    std::puts("  Test 4: String serialization");
    
    // Create a simple proof
    FieldVector numbers{F(5), F(10)};
    int num_bits = 4;
    AuxWitness aux = BuildAuxWitness(numbers, num_bits);
    struct proof original = ProveAuxConsistency(aux, numbers, num_bits);
    
    // Serialize to string
    std::string serialized_str = SerializeProofToString(original);
    assert(!serialized_str.empty());
    
    // Deserialize from string
    struct proof deserialized = DeserializeProofFromString(serialized_str);
    assert(deserialized.type == original.type);
    assert(deserialized.q_poly.size() == original.q_poly.size());
    
    std::puts("    String serialization OK");
  }
  
  // Test 5: Field element serialization
  {
    std::puts("  Test 5: Field element serialization");
    
    // fieldElement constructor takes (long long x, long long y)
    F elem1(123, 456);
    std::string str = FieldElementToString(elem1);
    F elem2 = StringToFieldElement(str);
    
    assert(elem1.real == elem2.real);
    assert(elem1.img == elem2.img);
    
    std::puts("    Field element serialization OK");
  }
  
  std::puts("test_proof_serialization OK");
  return 0;
}

