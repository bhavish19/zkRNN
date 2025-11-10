# zkRNN Implementation Plan - Runnable Phases

This document outlines the implementation phases, ensuring each phase produces runnable, testable code.

## Phase 0: Foundation Setup (Current State)
**Status**: ✅ Complete
- Placeholder structures created
- Basic test framework in place
- Code compiles and runs with dummy implementations

---

## Phase 1: Single Timestep Linear Proof (SUMMER Strategy)
**Goal**: Implement SumMatMul proof for linear operations in one RNN timestep
**Runnable Target**: Generate and verify MATMUL_PROOF for matrix-vector products

### Implementation Steps:

1. **Implement `ProveSumMatMul` in `sum_matmul_wrapper.cpp`**
   - Accept: weight matrices `W_h`, `W_x`, vectors `h_prev`, `x_t`
   - Compute: `Wh = W_h * h_prev`, `Wx = W_x * x_t`, `a_t = Wh + Wx + b1`
   - Combine into single matrix multiplication form
   - Use `_prove_matrix2matrix` from `main.cpp` to generate `MATMUL_PROOF`
   - Return `SumMatMulResult` with actual proof structure

2. **Update `sum_matmul_wrapper.h`**
   - Change return type to include `struct proof` instead of strings
   - Add field element types

3. **Create test: `test_sum_matmul.cpp`**
   - Generate small test matrices (e.g., 4x4, 8x8)
   - Call `ProveSumMatMul`
   - Verify proof structure is valid
   - Test compiles and runs

**Deliverables**:
- ✅ `ProveSumMatMul` generates real `MATMUL_PROOF`
- ✅ Test passes with small matrices
- ✅ Code compiles and runs

**Testing**:
```cpp
// test_sum_matmul.cpp
void test_sum_matmul_basic() {
  // Create test matrices and vectors
  std::vector<FieldVector> M = {W_h, W_x};
  FieldVector x = h_prev;
  FieldVector y = computed_output;
  
  SumMatMulResult result = ProveSumMatMul(M, x, y);
  assert(result.proof.type == MATMUL_PROOF);
}
```

---

## Phase 2: Auxiliary Witness Consistency Proof
**Goal**: Prove bit-decomposition consistency for auxiliary witness
**Runnable Target**: Generate and verify bit-decomposition proof

### Implementation Steps:

1. **Enhance `aux_witness.cpp`**
   - Add function: `struct proof ProveAuxWitnessConsistency(const AuxWitness &aux, const FieldVector &quantized)`
   - Use `_prove_bit_decomposition` logic from `main.cpp`
   - Generate `RANGE_PROOF` type proof
   - Verify bits match quantized values

2. **Update `aux_witness.h`**
   - Add proof function declaration
   - Include necessary headers (`GKR.h`)

3. **Create test: `test_aux_witness.cpp`**
   - Build `AuxWitness` for test vector
   - Generate consistency proof
   - Verify proof structure
   - Test compiles and runs

**Deliverables**:
- ✅ Auxiliary witness consistency proof generation
- ✅ Test passes
- ✅ Code compiles and runs

**Testing**:
```cpp
// test_aux_witness.cpp
void test_aux_witness_consistency() {
  FieldVector test_vec = {F(100), F(200), F(300)};
  AuxWitness aux = BuildAuxWitness(test_vec, 32);
  struct proof P = ProveAuxWitnessConsistency(aux, test_vec);
  assert(P.type == RANGE_PROOF);
}
```

---

## Phase 3: Activation GKR Proof (KAIZEN Strategy)
**Goal**: Implement GKR proof for non-linear activation functions
**Runnable Target**: Generate GKR_PROOF for tanh activation

### Implementation Steps:

1. **Implement `ProveActivationGKR` in `activation_gkr.cpp`**
   - Accept: `preact` (a_t), `aux` (AuxWitness), `h_out` (h_t)
   - Prepare GKR circuit data from inputs and bit-decompositions
   - Generate GKR circuit description for tanh function
   - Call `generate_GKR_proof` or similar
   - Generate auxiliary witness consistency proof
   - Return `struct proof` with type `GKR_PROOF`

2. **Update `activation_gkr.h`**
   - Change return type from `std::string` to `struct proof`
   - Include `GKR.h`

3. **Create test: `test_activation_gkr.cpp`**
   - Create test pre-activation vector
   - Build auxiliary witness
   - Compute tanh output
   - Generate GKR proof
   - Verify proof structure
   - Test compiles and runs

**Deliverables**:
- ✅ `ProveActivationGKR` generates real `GKR_PROOF`
- ✅ Test passes with small vectors
- ✅ Code compiles and runs

**Testing**:
```cpp
// test_activation_gkr.cpp
void test_activation_gkr_tanh() {
  FieldVector preact = {F(100), F(200)};
  AuxWitness aux = BuildAuxWitness(preact, 32);
  FieldVector h_out = compute_tanh(preact);  // actual computation
  
  struct proof P = ProveActivationGKR(preact, aux, h_out);
  assert(P.type == GKR_PROOF);
}
```

---

## Phase 4: Complete Single Timestep Proof (ProveLeaf)
**Goal**: Combine linear and non-linear proofs into complete per-timestep proof
**Runnable Target**: Generate complete proof for one RNN timestep

### Implementation Steps:

1. **Implement `ProveLeaf` in `prove_leaf.cpp`**
   - Accept: `TimeStepInput` (x_t, h_prev)
   - Compute linear operations: `a_t = W_h*h_prev + W_x*x_t + b1`
   - Generate SumMatMul proof for linear part
   - Compute non-linear: `h_t = tanh(a_t)`
   - Generate GKR proof for activation
   - Optionally: output layer `z_t = W_y*h_t + b2`, `yHat_t = softmax(z_t)`
   - Assemble `LeafResult` with actual proof structures

2. **Update `prove_leaf.h`**
   - Change `LeafResult` to include `struct proof` instead of strings
   - Add necessary headers

3. **Update `stream_orchestrator.cpp`**
   - Modify to work with actual proof structures (temporary: keep string-based accumulator for now)
   - Call updated `ProveLeaf`

4. **Create test: `test_prove_leaf.cpp`**
   - Create test timestep input
   - Call `ProveLeaf`
   - Verify `LeafResult` contains valid proofs
   - Test compiles and runs

**Deliverables**:
- ✅ `ProveLeaf` generates complete per-timestep proof
- ✅ Test passes
- ✅ Code compiles and runs

**Testing**:
```cpp
// test_prove_leaf.cpp
void test_prove_leaf_complete() {
  TimeStepInput input;
  input.x_t = {F(1), F(2)};
  input.h_prev = {F(0), F(0)};
  
  ProveParams params;
  params.Q = 32;
  
  LeafResult result = ProveLeaf(params, input);
  assert(!result.proofs.empty());
  // Verify proof structures
}
```

---

## Phase 5: Proof Serialization
**Goal**: Serialize proof structures for storage/transmission
**Runnable Target**: Serialize/deserialize proofs correctly

### Implementation Steps:

1. **Create `proof_serialize.h` and `.cpp`**
   - `serialize_proof(const struct proof&)` → `std::vector<uint8_t>`
   - `deserialize_proof(const std::vector<uint8_t>&)` → `struct proof`
   - Use existing `encode_*_proof` functions from `pol_verifier.cpp`
   - Handle all proof types (GKR, MATMUL, RANGE)

2. **Update `LeafResult` structure**
   - Add serialized proof field
   - Keep proof structures for verification

3. **Create test: `test_proof_serialize.cpp`**
   - Generate proof
   - Serialize
   - Deserialize
   - Verify deserialized proof matches original
   - Test compiles and runs

**Deliverables**:
- ✅ Proof serialization/deserialization works
- ✅ Test passes
- ✅ Code compiles and runs

---

## Phase 6: Recursive Aggregation (FA_Aggregate)
**Goal**: Implement KAIZEN-style recursive folding of k leaf proofs
**Runnable Target**: Aggregate multiple leaf proofs into one

### Implementation Steps:

1. **Implement `FA_Aggregate` in `ivc_adapter.cpp`**
   - Accept: `vector<LeafResult>` (k leaf proofs)
   - Extract proofs from each leaf
   - Use polynomial commitments for aggregation (reuse `prove_aggr` logic from `main.cpp`)
   - Generate aggregated proof
   - Compute new accumulator root
   - Return accumulator string (for now, can be commitment hash)

2. **Update `ivc_adapter.h`**
   - Change accumulator type to commitment root/hash
   - Include necessary headers

3. **Create test: `test_ivc_aggregate.cpp`**
   - Generate multiple leaf proofs
   - Call `FA_Aggregate`
   - Verify aggregated proof structure
   - Test compiles and runs

**Deliverables**:
- ✅ `FA_Aggregate` aggregates k leaf proofs
- ✅ Test passes
- ✅ Code compiles and runs

**Testing**:
```cpp
// test_ivc_aggregate.cpp
void test_aggregate_leafs() {
  std::vector<LeafResult> leaves;
  for (int i = 0; i < 3; i++) {
    leaves.push_back(ProveLeaf(params, inputs[i]));
  }
  
  std::string prev_acc = "ACC_INIT";
  std::string new_acc = FA_Aggregate(leaves, prev_acc);
  assert(new_acc != prev_acc);
}
```

---

## Phase 7: Final Proof Finalization
**Goal**: Implement `FinaliseProof` with actual proof structures
**Runnable Target**: Create final proof from aggregated proofs

### Implementation Steps:

1. **Update `FinalProof` structure in `finalise_proof.h`**
   - Change from strings to actual proof structures
   - Include `vector<struct proof>` for round transcripts
   - Add commitment root field

2. **Implement `FinaliseProof` in `finalise_proof.cpp`**
   - Accept: `vector<struct proof>` instead of strings
   - Aggregate all proofs using `verify_proof` logic
   - Generate final verification proof
   - Compute final accumulator root
   - Return `FinalProof` with actual proofs

3. **Update test: `test_finalise_verifier.cpp`**
   - Generate actual proofs (not dummy strings)
   - Call `FinaliseProof`
   - Verify structure
   - Test compiles and runs

**Deliverables**:
- ✅ `FinaliseProof` works with actual proofs
- ✅ Test passes
- ✅ Code compiles and runs

---

## Phase 8: Real Verification (VerifyProof)
**Goal**: Implement actual cryptographic verification
**Runnable Target**: Verify final proof correctly

### Implementation Steps:

1. **Implement `VerifyProof` in `verifier_stub.cpp`**
   - Accept: `FinalProof` with actual proof structures
   - Verify each proof in transcripts using `verify_proof()` from `pol_verifier.cpp`
   - Verify accumulator chain integrity
   - Check polynomial commitments
   - Return true/false based on actual verification

2. **Update `verifier_stub.h`**
   - Include necessary verification headers

3. **Update test: `test_finalise_verifier.cpp`**
   - Generate valid proof
   - Verify it passes
   - Generate invalid proof
   - Verify it fails
   - Test compiles and runs

**Deliverables**:
- ✅ `VerifyProof` performs real verification
- ✅ Test passes for valid proofs
- ✅ Test fails for invalid proofs
- ✅ Code compiles and runs

---

## Phase 9: Stream Orchestrator Integration
**Goal**: Integrate all components in StreamOrchestrator
**Runnable Target**: Complete end-to-end proof generation

### Implementation Steps:

1. **Update `StreamOrchestrator` in `stream_orchestrator.cpp`**
   - Use actual `ProveLeaf` with real proofs
   - Use actual `FA_Aggregate` with real aggregation
   - Use actual `FinaliseProof` and `VerifyProof`
   - Handle actual proof structures throughout

2. **Create integration test**
   - Test multiple timesteps
   - Test aggregation
   - Test finalization
   - Test verification
   - Test compiles and runs

**Deliverables**:
- ✅ Complete end-to-end flow works
- ✅ Integration test passes
- ✅ Code compiles and runs

---

## Phase 10: Optimization and Polish
**Goal**: Optimize and finalize implementation
**Runnable Target**: Optimized, production-ready code

### Implementation Steps:

1. **Memory optimization**
   - Proper cleanup of proof structures
   - Efficient serialization

2. **Error handling**
   - Add proper error checks
   - Meaningful error messages

3. **Documentation**
   - Code comments
   - API documentation

4. **Performance testing**
   - Benchmark with various sizes
   - Profile bottlenecks

**Deliverables**:
- ✅ Optimized code
- ✅ Proper error handling
- ✅ Documentation
- ✅ Performance benchmarks

---

## Implementation Order Summary

1. **Phase 1**: SumMatMul (linear) - Foundation for linear ops
2. **Phase 2**: Aux witness consistency - Foundation for non-linear
3. **Phase 3**: Activation GKR (non-linear) - Foundation for activations
4. **Phase 4**: ProveLeaf (complete timestep) - Combines linear + non-linear
5. **Phase 5**: Serialization - Needed for aggregation
6. **Phase 6**: FA_Aggregate (recursion) - Aggregates multiple timesteps
7. **Phase 7**: FinaliseProof - Final proof creation
8. **Phase 8**: VerifyProof - Real verification
9. **Phase 9**: StreamOrchestrator - End-to-end integration
10. **Phase 10**: Optimization - Polish

## Testing Strategy

Each phase will have:
- Unit tests for the specific component
- Integration tests with previous phases
- Compilation and runtime verification
- Small test cases that run quickly

## Dependencies Between Phases

```
Phase 1 (SumMatMul)
    ↓
Phase 2 (Aux Witness)
    ↓
Phase 3 (Activation GKR)
    ↓
Phase 4 (ProveLeaf) ──→ Phase 5 (Serialization)
    ↓                         ↓
Phase 6 (FA_Aggregate) ←──────┘
    ↓
Phase 7 (FinaliseProof)
    ↓
Phase 8 (VerifyProof)
    ↓
Phase 9 (StreamOrchestrator)
    ↓
Phase 10 (Optimization)
```

Each phase builds on previous phases while maintaining runnable code.

