cmake_minimum_required(VERSION 3.10)

add_executable(test_aux_witness
  ${CMAKE_SOURCE_DIR}/src/aux_witness.cpp
  ${CMAKE_SOURCE_DIR}/src/fieldElement.cpp
  ${CMAKE_SOURCE_DIR}/tests/test_aux_witness.cpp)

target_include_directories(test_aux_witness PRIVATE
  ${CMAKE_SOURCE_DIR}/src)


