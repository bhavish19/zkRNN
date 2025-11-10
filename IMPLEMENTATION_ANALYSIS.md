# zkRNN Implementation Analysis

## Overview
This document analyzes the existing zkRNN codebase and explains what needs to be implemented next for the placeholder code to work with the actual proof system.

## Final Implementation Summary

After completing implementation, the zkRNN prover now:

- Uses KAIZEN’s real cryptographic pipeline end-to-end while plugging in SUMMER-style linear proofs per timestep.
- Builds authentic leaf proofs (MATMUL + RANGE + LOGUP) for each RNN timestep via a rewritten `ProveLeaf`.
- Aggregates leaves with KAIZEN’s `prove_verification`, serialises/deserialises `struct proof`, and verifies with genuine `verify_matrix2matrix`, `verify_bit_decomposition`, and `verify_gkr`.
- Loads activation circuits on demand, aligns witnesses to KAIZEN layouts, and constructs auxiliary bit gadgets with two’s-complement masking.
- Restores KAIZEN’s quantised RNN forward helpers (tanh/softmax, exp batches), provides JSON snapshot dump/load, and allows inference replay via `--snapshot=…`.
- Instruments inference runs with prover/commit/logup/verifier timings plus peak memory, appending results to `bench_results_infer.csv`.
- Extends the build to pull in KAIZEN’s logup/poly-commit/expander code for both app and tests, which now cover snapshot loading, activation wiring, IVCAggregate, and orchestrator flows.

The sections below retain the historical analysis for reference.

## Protocol Design: Hybrid KAIZEN + SUMMER Approach

The implementation follows a **hybrid protocol** that combines:

1. **KAIZEN's Recursive Framework**: Tree-based IVC structure applied at **per-timestep granularity**
   - Each RNN timestep (OnInput(x_i)) generates a recursive proof
   - Contrasts with KAIZEN/SUMMER which apply recursion to full training iterations

2. **KAIZEN's Non-Linear Function Handling**: 
   - Quantize all real numbers to fixed-point field elements
   - Use GKR circuits for non-linear activations (tanh, softmax)
   - Provide bit-decompositions via auxiliary witness (AUX_i)
   - Use KAIZEN's sumcheck protocols for auxiliary witness consistency

3. **SUMMER's Linear Function Strategy**:
   - Reduce all linear operations to canonical matrix multiplication
   - Use SumMatMul protocol for efficient linear proofs
   - Combine multiple matrix-vector products into single proof

## Key Distinction: Per-Timestep vs Full-Iteration Recursion

**Current Original Code**: Proves full training iterations (forward + backward + update)

**Target Hybrid Protocol**: Proves individual RNN timesteps
- `ProveLeaf` should prove ONE timestep: x_t → a_t → h_t → z_t → yHat_t
- Recursion happens at timestep level, not training iteration level

## Current State: Placeholders Created

### 1. **Final Proof Finalization** (`finalise_proof.h` / `finalise_proof.cpp`)
**Current Implementation:**
- Simple placeholder that takes a final accumulator string and round transcripts
- Returns a `FinalProof` struct with `root_acc` and `transcripts` fields
- Currently just concatenates strings: `"FINAL_PROOF_OF_" + final_acc`

**What It Should Do:**
- Take the actual proof structures (`vector<struct proof>`) from multiple training rounds
- Convert them to a unified final proof format
- Serialize the proofs properly (the original code uses field elements `F`, not strings)
- Create a final accumulator/root that represents the entire training proof

### 2. **Verifier Stub** (`verifier_stub.h` / `verifier_stub.cpp`)
**Current Implementation:**
- Simple string-based verification: checks if `root_acc` contains `"FINAL_PROOF_OF_"`
- Returns true/false based on string pattern matching

**What It Should Do:**
- Actually verify the `FinalProof` structure
- Verify all the proof components (GKR proofs, polynomial commitments, etc.)
- Check the integrity of the accumulator chain
- Perform cryptographic verification of the proof commitments

### 3. **Stream Orchestrator** (`stream_orchestrator.cpp`)
**Current Implementation:**
- Uses placeholder `ProveLeaf` and `FA_Aggregate` functions
- Works with string-based accumulators
- Calls `FinaliseProof` and `VerifyProof` with string inputs

**What It Should Do:**
- Integrate with the actual RNN proof generation functions
- Handle real proof structures instead of strings
- Properly aggregate proofs using the actual aggregation logic
- **Per-timestep recursion**: Each `OnInput(x_i)` should prove ONE RNN timestep

### 4. **ProveLeaf** (`prove_leaf.cpp`)
**Current Implementation:**
- Placeholder that calls dummy `ProveActivationGKR` and returns string results

**What It Should Do (Per Single Timestep):**
1. **Linear Operations** (SUMMER approach):
   - Compute: `Wh = W_h * h_prev`, `Wx = W_x * x[t]`, `Zy = W_y * h[t]`
   - Combine all linear operations into single SumMatMul proof
   - Use `_prove_matrix2matrix` or similar for combined proof

2. **Non-Linear Operations** (KAIZEN approach):
   - For tanh activation: `a_t → h_t = tanh(a_t)`
   - For softmax: `z_t → yHat_t = softmax(z_t)`
   - Quantize inputs: Build `AuxWitness` with bit-decompositions
   - Generate GKR proof for activation circuit (NOT LogUp)
   - Use sumcheck protocols for auxiliary witness consistency

3. **Return LeafResult** with:
   - Actual proof structures (not strings)
   - Commitments for the timestep
   - Evaluation points for verification

### 5. **SumMatMul Wrapper** (`sum_matmul_wrapper.cpp`)
**Current Implementation:**
- Placeholder returning dummy strings

**What It Should Do:**
- Accept multiple matrix-vector products: `[W_h*h_prev, W_x*x_t, W_y*h_t]`
- Combine them into canonical matrix multiplication form
- Generate single `MATMUL_PROOF` using `_prove_matrix2matrix` logic
- Return actual proof structure with commitments

### 6. **Activation GKR** (`activation_gkr.cpp`)
**Current Implementation:**
- Placeholder returning dummy string

**What It Should Do:**
- Accept pre-activation vector `a_t` and auxiliary witness `AuxWitness`
- Generate GKR circuit for tanh/softmax function
- Use bit-decomposition from auxiliary witness
- Prove using `generate_GKR_proof` or similar
- Verify auxiliary witness consistency using sumcheck protocols
- Return actual `struct proof` with type `GKR_PROOF`

## Original Code Architecture

### Proof System Components

1. **Proof Types** (`GKR.h`):
   - `GKR_PROOF` - GKR-based proofs for circuit verification
   - `MATMUL_PROOF` - Matrix multiplication proofs
   - `LOOKUP_PROOF` - Lookup table proofs
   - `HASH_SUMCHECK` - Hash sumcheck proofs
   - `RANGE_PROOF` - Range/decomposition proofs
   - `LOGUP_PROOF` - LogUp proofs for exponentials

2. **Main Proof Structure** (`struct proof` in `GKR.h`):
   ```cpp
   struct proof {
       int type;
       vector<F> partial_eval, output, global_randomness;
       vector<vector<F>> randomness;
       vector<quadratic_poly> q_poly;
       vector<cubic_poly> c_poly;
       vector<layer_proof> proofs;
       // ... many more fields
   };
   ```

3. **Proof Workflow in `main.cpp`**:
   - **Forward Pass**: Generates proofs for RNN forward propagation
   - **Backward Pass**: Generates proofs for gradient computation
   - **Parameter Update**: Generates proofs for weight updates
   - **Aggregation**: Aggregates multiple proofs into one using polynomial commitments
   - **Recursive Verification**: Uses recursive proof verification at multiple levels

4. **Key Functions**:
   - `verify_proof(vector<proof> proofs)` - Verifies a batch of proofs and creates a verification proof
   - `prove_aggr()` - Aggregates proofs using polynomial commitments
   - `mimc_sumcheck()` - Creates sumcheck proofs for hash values

### RNN Training Proof Flow

The original code follows this pattern (from `main.cpp`):

```
For each training round:
  1. Forward pass → generate proofs → add to Transcript
  2. Backward pass → generate proofs → add to Transcript  
  3. Parameter update → generate proofs → add to Transcript
  4. Aggregate proofs at this level
  5. Create verification proof for aggregated batch
  6. Store for next-level aggregation

After all rounds:
  1. Final aggregation of all round proofs
  2. Final verification proof
  3. Verify the entire proof chain
```

## What Needs to Be Implemented

### Priority 1: Implement Per-Timestep Proof Generation

#### A. **Implement `ProveLeaf` for Single RNN Timestep**

**Signature:**
```cpp
LeafResult ProveLeaf(const ProveParams &p, const TimeStepInput &in);
```

**Input:**
- `in.x_t`: Input vector at timestep t (FieldVector)
- `in.h_prev`: Previous hidden state (FieldVector)

**Steps:**
1. **Linear Part** (SUMMER strategy):
   ```cpp
   // Compute pre-activation: a_t = W_h * h_prev + W_x * x_t + b1
   FieldVector Wh = mat_vec_mul(W_h, h_prev);
   FieldVector Wx = mat_vec_mul(W_x, x_t);
   FieldVector a_t = add(add(Wh, Wx), b1);
   
   // Combine into SumMatMul proof
   SumMatMulResult linear_proof = ProveSumMatMul({W_h, W_x}, {h_prev, x_t}, a_t);
   ```

2. **Non-Linear Part** (KAIZEN strategy):
   ```cpp
   // Quantize and build auxiliary witness
   AuxWitness aux = BuildAuxWitness(a_t, p.Q);
   
   // Compute activation: h_t = tanh(a_t)
   FieldVector h_t = compute_tanh(a_t);  // actual computation
   
   // Prove activation using GKR
   struct proof activation_proof = ProveActivationGKR(a_t, aux, h_t);
   ```

3. **Output Layer** (if needed for this timestep):
   ```cpp
   FieldVector Zy = mat_vec_mul(W_y, h_t);
   FieldVector z_t = add(Zy, b2);
   FieldVector yHat_t = compute_softmax(z_t);
   ```

4. **Assemble LeafResult**:
   ```cpp
   LeafResult result;
   result.step_proof_serialized = serialize_proofs({linear_proof, activation_proof});
   result.commitments = {get_commitment(linear_proof), get_commitment(activation_proof)};
   result.evals = {get_evals(linear_proof), get_evals(activation_proof)};
   result.accumulator_out = compute_accumulator(h_t);  // output state for next timestep
   ```

#### B. **Implement `ProveSumMatMul` (SUMMER Strategy)**

**Signature:**
```cpp
SumMatMulResult ProveSumMatMul(const std::vector<FieldVector> &M,
                                const FieldVector &x,
                                const FieldVector &y);
```

**Implementation:**
- Accept multiple matrix-vector products: `[W_h*h_prev, W_x*x_t]`
- Combine into canonical form: `M_combined * x_combined = y_combined`
- Use existing `_prove_matrix2matrix` from `main.cpp` (lines 483-520)
- Generate `MATMUL_PROOF` type proof
- Return commitments and evaluation points

**Key Function to Use:**
- `_prove_matrix2matrix(vector<vector<F>> M1, vector<vector<F>> M2, vector<F> r_eval, F previous_sum)` from `main.cpp`

#### C. **Implement `ProveActivationGKR` (KAIZEN Strategy)**

**Signature:**
```cpp
struct proof ProveActivationGKR(const FieldVector &preact,
                                 const AuxWitness &aux,
                                 const FieldVector &h_out);
```

**Implementation:**
1. **Build GKR Circuit**:
   - Create circuit for tanh/softmax function
   - Use bit-decomposition from `aux.bits`
   
2. **Generate GKR Proof**:
   ```cpp
   // Convert aux witness to format expected by GKR
   vector<F> gkr_data = prepare_gkr_data(preact, aux, h_out);
   vector<F> randomness = generate_randomness(...);
   
   // Generate proof using existing GKR infrastructure
   struct proof Pr = generate_GKR_proof(circuit_name, input_file, gkr_data, randomness, false);
   Pr.type = GKR_PROOF;
   ```

3. **Verify Auxiliary Witness Consistency**:
   - Use sumcheck protocols to verify bit-decomposition matches quantized values
   - Similar to `_prove_bit_decomposition` logic in `main.cpp` (lines 534-576)

**Key Functions to Use:**
- `generate_GKR_proof()` from `GKR.cpp`
- `_prove_bit_decomposition()` from `main.cpp` for auxiliary witness consistency
- Circuit generation for tanh/softmax functions

### Priority 2: Connect Placeholders to Real Proof System

#### A. **Implement `FinaliseProof` Properly**

**Current Signature:**
```cpp
FinalProof FinaliseProof(const std::string &final_acc,
                         const std::vector<std::string> &round_transcripts);
```

**Should Change To:**
```cpp
FinalProof FinaliseProof(const std::string &final_acc,
                         const std::vector<std::vector<struct proof>> &round_transcripts);
// OR convert to a serialized format that can be stored/transmitted
```

**Implementation Steps:**
1. Accept actual `vector<struct proof>` from each round (not strings)
2. Serialize proofs to a format suitable for final proof (could use field element encoding)
3. Create a final aggregated proof structure
4. Compute the final accumulator root (likely a Merkle root or hash of all proofs)
5. Return a `FinalProof` with the actual proof data

**Key Integration Points:**
- Use `verify_proof()` from `pol_verifier.cpp` to verify each round's proofs
- Use `prove_aggr()` logic to aggregate proofs if needed
- Use `mimc_sumcheck()` for hash verification

#### B. **Implement `VerifyProof` Properly**

**Current Signature:**
```cpp
bool VerifyProof(const FinalProof &proof);
```

**Implementation Steps:**
1. Deserialize the proof components from `FinalProof`
2. Verify each proof in `proof.transcripts` using the appropriate verifier
3. Verify the accumulator chain integrity
4. Check polynomial commitments if applicable
5. Verify hash commitments using MIMC

**Key Functions to Use:**
- `verify_proof()` from `pol_verifier.cpp` - for individual proof verification
- `verify_gkr()` from `pol_verifier.cpp` - for GKR proofs
- `verify_matrix2matrix()` - for matrix multiplication proofs
- Merkle tree verification for commitment verification

#### C. **Update `FinalProof` Structure**

**Current Structure:**
```cpp
struct FinalProof {
  std::string root_acc;
  std::vector<std::string> transcripts;
};
```

**Should Include:**
```cpp
struct FinalProof {
  std::string root_acc;  // Final accumulator/root hash
  std::vector<struct proof> round_proofs;  // Actual proof structures
  // OR serialized format
  commitment final_commitment;  // Final polynomial commitment
  std::vector<F> final_witness;  // Final witness data
};
```

### Priority 3: Integrate with Stream Orchestrator

The `StreamOrchestrator` needs to:
1. Use actual `ProveLeaf` implementation (per-timestep proof)
2. Use actual `FA_Aggregate` implementation (KAIZEN-style recursive folding)
3. Work with real proof structures, not strings

**Key Integration Points:**
- `ProveLeaf` should prove ONE RNN timestep (not full iteration)
- `FA_Aggregate` should fold k leaf proofs using KAIZEN's recursive framework
- Accumulator should be actual commitment roots from polynomial commitments
- Use `prove_aggr()` or similar for aggregation

### Priority 4: Proof Serialization

Since proofs contain complex structures (polynomials, vectors of field elements), you need:

1. **Serialization Functions:**
   - `serialize_proof(const struct proof&)` → `std::vector<uint8_t>` or `std::string`
   - `deserialize_proof(const std::vector<uint8_t>&)` → `struct proof`

2. **Encoding Functions** (already exist in `pol_verifier.cpp`):
   - `encode_gkr_proof()` - for GKR proofs
   - `encode_m2m_proof()` - for matrix multiplication proofs
   - `encode_hash_proof()` - for hash proofs
   - These return `vector<F>` which can be serialized

### Priority 5: Test Integration

The test file `test_finalise_verifier.cpp` currently:
- Creates dummy string inputs
- Calls `FinaliseProof` and `VerifyProof`
- Asserts verification succeeds

**To Make It Work:**
1. Generate actual proofs using the RNN proof functions
2. Convert to `FinalProof` format
3. Verify using real verification logic

## Implementation Roadmap

### Phase 1: Per-Timestep Proof Generation (Core Protocol)
1. ✅ Create placeholder structures (DONE)
2. ⏳ Implement `ProveSumMatMul` - Combine linear operations into single proof
3. ⏳ Implement `ProveActivationGKR` - GKR proof for non-linear activations
4. ⏳ Implement `BuildAuxWitness` consistency proofs (sumcheck for bits)
5. ⏳ Implement `ProveLeaf` - Complete per-timestep proof generation
6. ⏳ Test single timestep proof generation

### Phase 2: Recursive Aggregation (KAIZEN Framework)
1. ⏳ Implement `FA_Aggregate` - Recursive folding of k leaf proofs
2. ⏳ Integrate with polynomial commitments for accumulator
3. ⏳ Implement tree-based IVC structure
4. ⏳ Test multi-timestep aggregation

### Phase 3: Finalization and Verification
1. ⏳ Update `FinalProof` to accept `vector<struct proof>` instead of strings
2. ⏳ Implement basic serialization for proof structures
3. ⏳ Update `FinaliseProof` to work with actual proofs
4. ⏳ Update `VerifyProof` to perform real verification
5. ⏳ Integrate with `verify_proof()` for verification

### Phase 4: Optimization and Integration
1. ⏳ Optimize proof serialization
2. ⏳ Batch verification where possible
3. ⏳ Memory management for large proofs
4. ⏳ Handle multiple proof types (GKR, MATMUL, etc.)

## Key Files to Implement/Modify

### Core Protocol Implementation (Priority 1)

1. **`src/prove_leaf.cpp`**:
   - Implement per-timestep proof generation
   - Call `ProveSumMatMul` for linear operations
   - Call `ProveActivationGKR` for non-linear operations
   - Assemble `LeafResult` with actual proof structures

2. **`src/sum_matmul_wrapper.cpp`**:
   - Implement SUMMER-style SumMatMul protocol
   - Combine multiple matrix-vector products
   - Use `_prove_matrix2matrix` from `main.cpp`
   - Return actual `MATMUL_PROOF` structure

3. **`src/activation_gkr.cpp`**:
   - Implement KAIZEN-style GKR proof for activations
   - Use `generate_GKR_proof` for tanh/softmax circuits
   - Verify auxiliary witness consistency
   - Return actual `GKR_PROOF` structure

4. **`src/aux_witness.cpp`** (if needed):
   - Add sumcheck proof for bit-decomposition consistency
   - Similar to `_prove_bit_decomposition` in `main.cpp`

### Recursive Framework (Priority 2)

5. **`src/ivc_adapter.cpp`**:
   - Implement KAIZEN-style recursive folding
   - Use polynomial commitments for aggregation
   - Implement tree-based IVC structure

### Finalization (Priority 3)

6. **`src/finalise_proof.h` and `.cpp`**:
   - Change signatures to accept `vector<struct proof>`
   - Implement actual proof aggregation/finalization
   - Add serialization logic

7. **`src/verifier_stub.h` and `.cpp`**:
   - Implement real verification logic
   - Call `verify_proof()`, `verify_gkr()`, etc.
   - Verify commitments and hash chains

8. **`src/stream_orchestrator.cpp`**:
   - Update to work with real proof structures
   - Integrate with actual `ProveLeaf` and `FA_Aggregate`

### Testing (Priority 4)

9. **`tests/test_finalise_verifier.cpp`**:
   - Generate real proofs for single timestep
   - Test actual proof verification
   - Test recursive aggregation

## Dependencies

The implementation will need:
- `pol_verifier.h/cpp` - For `verify_proof()` and related functions
- `GKR.h/cpp` - For proof structures and GKR proof generation
- `poly_commit.h/cpp` - For polynomial commitments
- `mimc.h/cpp` - For hash verification
- `config_pc.hpp` - For field element type `F`

## Important Implementation Notes

### Protocol-Specific Requirements

1. **Do NOT use LogUp for non-linear functions**:
   - The original code uses `exp_tanh` and `exp_softmax` with LogUp lookups
   - **Our hybrid protocol requires GKR + auxiliary witness instead**
   - Remove LogUp dependencies from activation proofs

2. **Per-Timestep Granularity**:
   - `ProveLeaf` proves ONE timestep: `x_t → a_t → h_t → z_t → yHat_t`
   - NOT a full training iteration (forward + backward + update)
   - This is the key distinction from KAIZEN/SUMMER

3. **Linear Operations**:
   - Combine all matrix-vector products: `W_h*h_prev`, `W_x*x_t`, `W_y*h_t`
   - Use single SumMatMul proof (SUMMER strategy)
   - Reference: `_prove_matrix2matrix` in `main.cpp`

4. **Non-Linear Operations**:
   - Use GKR circuits (KAIZEN strategy)
   - Provide bit-decomposition via `AuxWitness`
   - Prove auxiliary witness consistency with sumcheck
   - Reference: `_prove_bit_decomposition` in `main.cpp`

### Technical Details

- The original code uses field elements (`F`) throughout, not strings
- Proofs are complex structures with polynomials, commitments, and randomness
- The verification process involves multiple rounds of sumcheck protocols
- Proof aggregation uses polynomial commitments (likely Orion or similar PCS)
- The system supports recursive proof verification for efficiency
- Quantization uses `Q` bits (default 32, configurable via `ProveParams`)

### Existing Code to Reuse

- `_prove_matrix2matrix()` - For linear operations (lines 483-520 in `main.cpp`)
- `_prove_bit_decomposition()` - For auxiliary witness consistency (lines 534-576)
- `generate_GKR_proof()` - For GKR circuit proofs (in `GKR.cpp`)
- `prove_aggr()` - For proof aggregation (in `main.cpp`)
- `verify_proof()` - For proof verification (in `pol_verifier.cpp`)

### Code to Avoid/Modify

- **LogUp proofs**: The original code uses `prove_logup()` for exponentials
  - **We should NOT use this** - use GKR instead
- **Full iteration proofs**: Original `prove_rnn_forward()` proves entire sequence
  - **We need per-timestep proofs** instead

