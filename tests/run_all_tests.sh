#!/bin/bash
# Test runner script for zkRNN implementation
# Run all tests and report results

set -e  # Exit on error

echo "=========================================="
echo "zkRNN Test Suite"
echo "=========================================="
echo ""

# Build all tests first
echo "Building all tests..."
cd "$(dirname "$0")/.."
./build.sh
echo ""

# Test counter
PASSED=0
FAILED=0

# Function to run a test
run_test() {
    local test_name=$1
    local test_executable="./build/tests/$test_name"
    
    echo "----------------------------------------"
    echo "Running: $test_name"
    echo "----------------------------------------"
    
    if [ -f "$test_executable" ]; then
        if "$test_executable" 2>&1; then
            echo "✓ $test_name PASSED"
            ((PASSED++))
        else
            echo "✗ $test_name FAILED"
            ((FAILED++))
        fi
    else
        echo "✗ $test_name NOT FOUND (build failed?)"
        ((FAILED++))
    fi
    echo ""
}

# Run tests in order of dependency
# Phase 1: Basic components
run_test "test_aux_witness"
run_test "test_sum_matmul_wrapper"
run_test "test_sum_matmul"
run_test "test_aux_sumcheck"

# Phase 2: Combined components
run_test "test_activation_gkr"
run_test "test_prove_leaf"

# Phase 3: Aggregation and finalization
run_test "test_ivc_adapter"
run_test "test_proof_serialization"
run_test "test_finalise_verifier"

# Phase 4: End-to-end
run_test "test_stream_orchestrator"

# Summary
echo "=========================================="
echo "Test Summary"
echo "=========================================="
echo "Passed: $PASSED"
echo "Failed: $FAILED"
echo "Total:  $((PASSED + FAILED))"
echo ""

if [ $FAILED -eq 0 ]; then
    echo "✓ All tests passed!"
    exit 0
else
    echo "✗ Some tests failed"
    exit 1
fi

