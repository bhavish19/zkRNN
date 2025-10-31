find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)

aux_source_directory(. SRC)
include_directories(.)
add_executable(zkpot ${SRC})
target_link_libraries(zkpot PRIVATE OpenMP::OpenMP_CXX Threads::Threads gmp XKCP)
