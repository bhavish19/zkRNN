# zkRNN Testing Guide

This guide explains how to test the zkRNN implementation.

## Available Tests

### Phase 1: Basic Components

1. **test_aux_witness** - Tests auxiliary witness construction (bit-decomposition)
   ```bash
   ./build/tests/test_aux_witness
   ```

2. **test_sum_matmul_wrapper** - Tests SumMatMul proof generation (linear operations)
   ```bash
   ./build/tests/test_sum_matmul_wrapper
   ```

3. **test_sum_matmul** - Comprehensive SumMatMul tests
   ```bash
   ./build/tests/test_sum_matmul
   ```

4. **test_aux_sumcheck** - Tests auxiliary witness consistency proof (RANGE_PROOF)
   ```bash
   ./build/tests/test_aux_sumcheck
   ```

### Phase 2: Combined Components

5. **test_activation_gkr** - Tests activation function GKR proof (tanh/softmax)
   ```bash
   ./build/tests/test_activation_gkr
   ```

6. **test_prove_leaf** - Tests complete leaf proof (linear + non-linear)
   ```bash
   ./build/tests/test_prove_leaf
   ```

### Phase 3: Aggregation and Finalization

7. **test_ivc_adapter** - Tests FA_Aggregate (recursive folding of k leaf proofs)
   ```bash
   ./build/tests/test_ivc_adapter
   ```

8. **test_proof_serialization** - Tests proof serialization/deserialization
   ```bash
   ./build/tests/test_proof_serialization
   ```

9. **test_finalise_verifier** - Tests FinaliseProof and VerifyProof
   ```bash
   ./build/tests/test_finalise_verifier
   ```

### Phase 4: End-to-End

10. **test_stream_orchestrator** - Tests complete RNN inference flow
    ```bash
    ./build/tests/test_stream_orchestrator
    ```

## Running All Tests

### Option 1: Use the test runner script
```bash
./tests/run_all_tests.sh
```

### Option 2: Run tests individually
```bash
# Build all tests
./build.sh

# Run each test
./build/tests/test_aux_witness
./build/tests/test_sum_matmul_wrapper
./build/tests/test_sum_matmul
./build/tests/test_aux_sumcheck
./build/tests/test_activation_gkr
./build/tests/test_prove_leaf
./build/tests/test_ivc_adapter
./build/tests/test_proof_serialization
./build/tests/test_finalise_verifier
./build/tests/test_stream_orchestrator
```

### Option 3: Run tests with CMake
```bash
cd build
ctest --output-on-failure
```

## Expected Output

Each test should print:
- Test name and status
- Any debug output from the implementation
- "test_<name> OK" on success

## Troubleshooting

### Segmentation Fault
If you see a segmentation fault:
1. Check that all dependencies are linked correctly in `tests/CMakeLists.txt`
2. Ensure all required initialization functions are called (e.g., `init_hash()`, `init_SHA()`)
3. Check that proof structures are properly initialized before encoding

### Linker Errors
If you see undefined reference errors:
1. Check `tests/CMakeLists.txt` to ensure all required `.cpp` files are included
2. Rebuild: `./build.sh`

### Test Failures
If a test fails:
1. Check the test output for error messages
2. Verify that the test inputs are correct
3. Check that the implementation matches the expected behavior

## Testing Strategy

1. **Unit Tests**: Run individual component tests first (Phase 1)
2. **Integration Tests**: Run combined component tests (Phase 2)
3. **System Tests**: Run aggregation and finalization tests (Phase 3)
4. **End-to-End Tests**: Run the complete orchestrator test (Phase 4)

## Next Steps

After all tests pass:
- Phase 8: Implement real VerifyProof with cryptographic verification
- Phase 9: Integrate all components in StreamOrchestrator
- Phase 10: Optimization and polish

