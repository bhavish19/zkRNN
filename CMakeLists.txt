cmake_minimum_required(VERSION 3.10)
project(ProofOfLearning)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-w -march=native -mavx -lgmp -msha")
set(CMAKE_C_FLAGS "-w -march=native -mavx -lgmp -msha")


link_directories(lib)

include_directories(src)
include_directories(lib)

add_compile_definitions(PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

add_library(XKCP STATIC IMPORTED GLOBAL)
set_target_properties(XKCP PROPERTIES
  IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/lib/libXKCP.a"
  IMPORTED_LINK_INTERFACE_LANGUAGES "C"
  INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/lib/libXKCP.a.headers"
)

add_subdirectory(src)
add_subdirectory(tests)