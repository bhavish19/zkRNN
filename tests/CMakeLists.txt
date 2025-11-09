cmake_minimum_required(VERSION 3.10)

find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)

set(LOGUP_SOURCES
  ${CMAKE_SOURCE_DIR}/src/logup.cpp
  ${CMAKE_SOURCE_DIR}/src/poly_commit.cpp
  ${CMAKE_SOURCE_DIR}/src/linear_code_encode.cpp
  ${CMAKE_SOURCE_DIR}/src/merkle_tree.cpp
  ${CMAKE_SOURCE_DIR}/src/expander.cpp
  ${CMAKE_SOURCE_DIR}/tests/logup_globals.cpp
)

set(COMMON_PROOF_LIBS
  ${GMP_LIBRARIES}
  ${MPFR_LIBRARIES}
  OpenMP::OpenMP_CXX
  Threads::Threads
)

add_executable(test_aux_witness
  ${CMAKE_SOURCE_DIR}/src/aux_witness.cpp
  ${CMAKE_SOURCE_DIR}/src/fieldElement.cpp
  ${CMAKE_SOURCE_DIR}/tests/test_aux_witness.cpp)

target_include_directories(test_aux_witness PRIVATE
  ${CMAKE_SOURCE_DIR}/src)

  add_executable(test_prove_leaf
  ${CMAKE_SOURCE_DIR}/src/prove_leaf.cpp
  ${CMAKE_SOURCE_DIR}/src/sum_matmul_wrapper.cpp
  ${CMAKE_SOURCE_DIR}/src/activation_gkr.cpp
  ${LOGUP_SOURCES}
  ${CMAKE_SOURCE_DIR}/src/activation_witness.cpp
  ${CMAKE_SOURCE_DIR}/src/activation_circuit_loader.cpp
  ${CMAKE_SOURCE_DIR}/src/aux_sumcheck.cpp
  ${CMAKE_SOURCE_DIR}/src/aux_witness.cpp
  ${CMAKE_SOURCE_DIR}/src/GKR.cpp
  ${CMAKE_SOURCE_DIR}/src/prover.cpp
  ${CMAKE_SOURCE_DIR}/src/verifier.cpp
  ${CMAKE_SOURCE_DIR}/src/circuit.cpp
  ${CMAKE_SOURCE_DIR}/src/timer.cpp
  ${CMAKE_SOURCE_DIR}/src/pol_verifier.cpp
  ${CMAKE_SOURCE_DIR}/src/fieldElement.cpp
  ${CMAKE_SOURCE_DIR}/src/proof_utils.cpp
  ${CMAKE_SOURCE_DIR}/src/utils.cpp
  ${CMAKE_SOURCE_DIR}/src/mimc.cpp
  ${CMAKE_SOURCE_DIR}/src/polynomial.cpp
  ${CMAKE_SOURCE_DIR}/src/quantization.cpp
  ${CMAKE_SOURCE_DIR}/tests/test_prove_leaf.cpp)

target_include_directories(test_prove_leaf PRIVATE
  ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(test_prove_leaf
  PRIVATE
    ${COMMON_PROOF_LIBS}
    XKCP)

add_executable(test_activation_gkr
  ${CMAKE_SOURCE_DIR}/src/activation_gkr.cpp
  ${LOGUP_SOURCES}
  ${CMAKE_SOURCE_DIR}/src/activation_witness.cpp
  ${CMAKE_SOURCE_DIR}/src/activation_circuit_loader.cpp
  ${CMAKE_SOURCE_DIR}/src/aux_sumcheck.cpp
  ${CMAKE_SOURCE_DIR}/src/aux_witness.cpp
  ${CMAKE_SOURCE_DIR}/src/GKR.cpp
  ${CMAKE_SOURCE_DIR}/src/prover.cpp
  ${CMAKE_SOURCE_DIR}/src/verifier.cpp
  ${CMAKE_SOURCE_DIR}/src/circuit.cpp
  ${CMAKE_SOURCE_DIR}/src/timer.cpp
  ${CMAKE_SOURCE_DIR}/src/fieldElement.cpp
  ${CMAKE_SOURCE_DIR}/src/proof_utils.cpp
  ${CMAKE_SOURCE_DIR}/src/utils.cpp
  ${CMAKE_SOURCE_DIR}/src/mimc.cpp
  ${CMAKE_SOURCE_DIR}/src/polynomial.cpp
  ${CMAKE_SOURCE_DIR}/src/quantization.cpp
  ${CMAKE_SOURCE_DIR}/tests/test_activation_gkr.cpp)

target_include_directories(test_activation_gkr PRIVATE
  ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(test_activation_gkr
  PRIVATE
    ${COMMON_PROOF_LIBS}
    XKCP)

add_executable(test_activation_loader
  ${CMAKE_SOURCE_DIR}/src/activation_circuit_loader.cpp
  ${CMAKE_SOURCE_DIR}/tests/test_activation_loader.cpp)

target_include_directories(test_activation_loader PRIVATE
  ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(test_activation_loader
  PRIVATE
    ${COMMON_PROOF_LIBS})

add_executable(test_activation_witness
  ${CMAKE_SOURCE_DIR}/src/activation_witness.cpp
  ${CMAKE_SOURCE_DIR}/src/activation_circuit_loader.cpp
  ${CMAKE_SOURCE_DIR}/src/aux_witness.cpp
  ${CMAKE_SOURCE_DIR}/src/fieldElement.cpp
  ${CMAKE_SOURCE_DIR}/src/proof_utils.cpp
  ${CMAKE_SOURCE_DIR}/src/utils.cpp
  ${CMAKE_SOURCE_DIR}/src/mimc.cpp
  ${CMAKE_SOURCE_DIR}/src/polynomial.cpp
  ${CMAKE_SOURCE_DIR}/src/quantization.cpp
  ${CMAKE_SOURCE_DIR}/tests/test_activation_witness.cpp)

target_include_directories(test_activation_witness PRIVATE
  ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(test_activation_witness
  PRIVATE
    ${COMMON_PROOF_LIBS})

  add_executable(test_ivc_adapter
  ${CMAKE_SOURCE_DIR}/src/ivc_adapter.cpp
  ${CMAKE_SOURCE_DIR}/src/proof_serialization.cpp
  ${CMAKE_SOURCE_DIR}/src/pol_verifier.cpp
  ${CMAKE_SOURCE_DIR}/src/GKR.cpp
  ${CMAKE_SOURCE_DIR}/src/prover.cpp
  ${CMAKE_SOURCE_DIR}/src/verifier.cpp
  ${CMAKE_SOURCE_DIR}/src/circuit.cpp
  ${CMAKE_SOURCE_DIR}/src/timer.cpp
  ${CMAKE_SOURCE_DIR}/src/prove_leaf.cpp
  ${CMAKE_SOURCE_DIR}/src/sum_matmul_wrapper.cpp
  ${CMAKE_SOURCE_DIR}/src/activation_gkr.cpp
  ${LOGUP_SOURCES}
  ${CMAKE_SOURCE_DIR}/src/activation_witness.cpp
  ${CMAKE_SOURCE_DIR}/src/activation_circuit_loader.cpp
  ${CMAKE_SOURCE_DIR}/src/aux_sumcheck.cpp
  ${CMAKE_SOURCE_DIR}/src/aux_witness.cpp
  ${CMAKE_SOURCE_DIR}/src/fieldElement.cpp
  ${CMAKE_SOURCE_DIR}/src/proof_utils.cpp
  ${CMAKE_SOURCE_DIR}/src/utils.cpp
  ${CMAKE_SOURCE_DIR}/src/mimc.cpp
  ${CMAKE_SOURCE_DIR}/src/polynomial.cpp
  ${CMAKE_SOURCE_DIR}/src/quantization.cpp
  ${CMAKE_SOURCE_DIR}/tests/test_ivc_adapter.cpp)

target_include_directories(test_ivc_adapter PRIVATE
  ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(test_ivc_adapter
  PRIVATE
    ${COMMON_PROOF_LIBS}
    XKCP)

add_executable(test_stream_orchestrator
  ${CMAKE_SOURCE_DIR}/src/stream_orchestrator.cpp
  ${CMAKE_SOURCE_DIR}/src/ivc_adapter.cpp
  ${CMAKE_SOURCE_DIR}/src/proof_serialization.cpp
  ${CMAKE_SOURCE_DIR}/src/pol_verifier.cpp
  ${CMAKE_SOURCE_DIR}/src/GKR.cpp
  ${CMAKE_SOURCE_DIR}/src/prover.cpp
  ${CMAKE_SOURCE_DIR}/src/verifier.cpp
  ${CMAKE_SOURCE_DIR}/src/circuit.cpp
  ${CMAKE_SOURCE_DIR}/src/timer.cpp
  ${CMAKE_SOURCE_DIR}/src/prove_leaf.cpp
  ${CMAKE_SOURCE_DIR}/src/sum_matmul_wrapper.cpp
  ${CMAKE_SOURCE_DIR}/src/activation_gkr.cpp
  ${LOGUP_SOURCES}
  ${CMAKE_SOURCE_DIR}/src/RNN.cpp
  ${CMAKE_SOURCE_DIR}/src/activation_witness.cpp
  ${CMAKE_SOURCE_DIR}/src/activation_circuit_loader.cpp
  ${CMAKE_SOURCE_DIR}/src/aux_sumcheck.cpp
  ${CMAKE_SOURCE_DIR}/src/aux_witness.cpp
  ${CMAKE_SOURCE_DIR}/src/fieldElement.cpp
  ${CMAKE_SOURCE_DIR}/src/finalise_proof.cpp
  ${CMAKE_SOURCE_DIR}/src/verifier_stub.cpp
  ${CMAKE_SOURCE_DIR}/src/proof_utils.cpp
  ${CMAKE_SOURCE_DIR}/src/utils.cpp
  ${CMAKE_SOURCE_DIR}/src/mimc.cpp
  ${CMAKE_SOURCE_DIR}/src/polynomial.cpp
  ${CMAKE_SOURCE_DIR}/src/quantization.cpp
  ${CMAKE_SOURCE_DIR}/tests/test_stream_orchestrator.cpp
)

target_include_directories(test_stream_orchestrator PRIVATE
  ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(test_stream_orchestrator
  PRIVATE
    ${COMMON_PROOF_LIBS}
    XKCP)

target_compile_definitions(test_stream_orchestrator PRIVATE DEBUG_LOGUP=1)

add_executable(test_sum_matmul_wrapper
  ${CMAKE_SOURCE_DIR}/src/sum_matmul_wrapper.cpp
  ${CMAKE_SOURCE_DIR}/src/proof_utils.cpp
  ${CMAKE_SOURCE_DIR}/src/fieldElement.cpp
  ${CMAKE_SOURCE_DIR}/src/utils.cpp
  ${CMAKE_SOURCE_DIR}/src/mimc.cpp
  ${CMAKE_SOURCE_DIR}/src/polynomial.cpp
  ${CMAKE_SOURCE_DIR}/src/quantization.cpp
  ${CMAKE_SOURCE_DIR}/tests/test_sum_matmul_wrapper.cpp)

target_include_directories(test_sum_matmul_wrapper PRIVATE
  ${CMAKE_SOURCE_DIR}/src)

add_executable(test_sum_matmul
  ${CMAKE_SOURCE_DIR}/src/sum_matmul_wrapper.cpp
  ${CMAKE_SOURCE_DIR}/src/proof_utils.cpp
  ${CMAKE_SOURCE_DIR}/src/fieldElement.cpp
  ${CMAKE_SOURCE_DIR}/src/utils.cpp
  ${CMAKE_SOURCE_DIR}/src/mimc.cpp
  ${CMAKE_SOURCE_DIR}/src/polynomial.cpp
  ${CMAKE_SOURCE_DIR}/src/quantization.cpp
  ${CMAKE_SOURCE_DIR}/tests/test_sum_matmul.cpp)

target_include_directories(test_sum_matmul PRIVATE
  ${CMAKE_SOURCE_DIR}/src)

target_include_directories(test_sum_matmul_wrapper PRIVATE
  ${CMAKE_SOURCE_DIR}/src)

add_executable(test_aux_sumcheck
  ${CMAKE_SOURCE_DIR}/src/aux_sumcheck.cpp
  ${CMAKE_SOURCE_DIR}/src/fieldElement.cpp
  ${CMAKE_SOURCE_DIR}/src/aux_witness.cpp
  ${CMAKE_SOURCE_DIR}/src/proof_utils.cpp
  ${CMAKE_SOURCE_DIR}/src/utils.cpp
  ${CMAKE_SOURCE_DIR}/src/mimc.cpp
  ${CMAKE_SOURCE_DIR}/src/polynomial.cpp
  ${CMAKE_SOURCE_DIR}/src/quantization.cpp
  ${CMAKE_SOURCE_DIR}/tests/test_aux_sumcheck.cpp)

target_include_directories(test_aux_sumcheck PRIVATE
  ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(test_aux_sumcheck
  ${GMP_LIBRARIES}
  ${MPFR_LIBRARIES})

add_executable(test_finalise_verifier
  ${CMAKE_SOURCE_DIR}/src/finalise_proof.cpp
  ${CMAKE_SOURCE_DIR}/src/verifier_stub.cpp
  ${CMAKE_SOURCE_DIR}/src/proof_serialization.cpp
  ${CMAKE_SOURCE_DIR}/src/pol_verifier.cpp
  ${CMAKE_SOURCE_DIR}/src/GKR.cpp
  ${CMAKE_SOURCE_DIR}/src/prover.cpp
  ${CMAKE_SOURCE_DIR}/src/verifier.cpp
  ${CMAKE_SOURCE_DIR}/src/circuit.cpp
  ${CMAKE_SOURCE_DIR}/src/timer.cpp
  ${CMAKE_SOURCE_DIR}/src/proof_utils.cpp
  ${CMAKE_SOURCE_DIR}/src/fieldElement.cpp
  ${CMAKE_SOURCE_DIR}/src/polynomial.cpp
  ${CMAKE_SOURCE_DIR}/src/mimc.cpp
  ${CMAKE_SOURCE_DIR}/src/utils.cpp
  ${CMAKE_SOURCE_DIR}/tests/test_finalise_verifier.cpp)

target_include_directories(test_finalise_verifier PRIVATE
  ${CMAKE_SOURCE_DIR}/src)

add_executable(test_proof_serialization
  ${CMAKE_SOURCE_DIR}/src/proof_serialization.cpp
  ${CMAKE_SOURCE_DIR}/src/pol_verifier.cpp
  ${CMAKE_SOURCE_DIR}/src/GKR.cpp
  ${CMAKE_SOURCE_DIR}/src/prover.cpp
  ${CMAKE_SOURCE_DIR}/src/verifier.cpp
  ${CMAKE_SOURCE_DIR}/src/circuit.cpp
  ${CMAKE_SOURCE_DIR}/src/timer.cpp
  ${CMAKE_SOURCE_DIR}/src/sum_matmul_wrapper.cpp
  ${CMAKE_SOURCE_DIR}/src/activation_gkr.cpp
  ${LOGUP_SOURCES}
  ${CMAKE_SOURCE_DIR}/src/activation_witness.cpp
  ${CMAKE_SOURCE_DIR}/src/activation_circuit_loader.cpp
  ${CMAKE_SOURCE_DIR}/src/aux_sumcheck.cpp
  ${CMAKE_SOURCE_DIR}/src/aux_witness.cpp
  ${CMAKE_SOURCE_DIR}/src/proof_utils.cpp
  ${CMAKE_SOURCE_DIR}/src/utils.cpp
  ${CMAKE_SOURCE_DIR}/src/mimc.cpp
  ${CMAKE_SOURCE_DIR}/src/polynomial.cpp
  ${CMAKE_SOURCE_DIR}/src/quantization.cpp
  ${CMAKE_SOURCE_DIR}/src/fieldElement.cpp
  ${CMAKE_SOURCE_DIR}/tests/test_proof_serialization.cpp
)

target_include_directories(test_proof_serialization PRIVATE
  ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(test_proof_serialization
  PRIVATE
    ${COMMON_PROOF_LIBS}
    XKCP)

target_compile_definitions(test_proof_serialization PRIVATE DEBUG_LOGUP=1)

add_executable(test_witness_snapshot
  ${CMAKE_SOURCE_DIR}/src/witness_snapshot.cpp
  ${CMAKE_SOURCE_DIR}/src/proof_serialization.cpp
  ${CMAKE_SOURCE_DIR}/src/pol_verifier.cpp
  ${CMAKE_SOURCE_DIR}/src/activation_gkr.cpp
  ${LOGUP_SOURCES}
  ${CMAKE_SOURCE_DIR}/src/activation_witness.cpp
  ${CMAKE_SOURCE_DIR}/src/activation_circuit_loader.cpp
  ${CMAKE_SOURCE_DIR}/src/aux_sumcheck.cpp
  ${CMAKE_SOURCE_DIR}/src/aux_witness.cpp
  ${CMAKE_SOURCE_DIR}/src/GKR.cpp
  ${CMAKE_SOURCE_DIR}/src/prover.cpp
  ${CMAKE_SOURCE_DIR}/src/verifier.cpp
  ${CMAKE_SOURCE_DIR}/src/circuit.cpp
  ${CMAKE_SOURCE_DIR}/src/timer.cpp
  ${CMAKE_SOURCE_DIR}/src/fieldElement.cpp
  ${CMAKE_SOURCE_DIR}/src/proof_utils.cpp
  ${CMAKE_SOURCE_DIR}/src/utils.cpp
  ${CMAKE_SOURCE_DIR}/src/mimc.cpp
  ${CMAKE_SOURCE_DIR}/src/polynomial.cpp
  ${CMAKE_SOURCE_DIR}/src/quantization.cpp
  ${CMAKE_SOURCE_DIR}/tests/test_witness_snapshot.cpp)

target_include_directories(test_witness_snapshot PRIVATE
  ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(test_witness_snapshot
  PRIVATE
    ${COMMON_PROOF_LIBS}
    XKCP)

target_compile_definitions(test_witness_snapshot PRIVATE DEBUG_AUX_SUMCHECK=1)