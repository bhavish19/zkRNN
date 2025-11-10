find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)

# Get all source files except proof_utils.cpp (it's only for tests to avoid conflicts with main.cpp)
file(GLOB SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
list(REMOVE_ITEM SRC "${CMAKE_CURRENT_SOURCE_DIR}/proof_utils.cpp")

include_directories(.)
add_executable(zkpot ${SRC})
target_link_libraries(zkpot PRIVATE OpenMP::OpenMP_CXX Threads::Threads gmp XKCP)
