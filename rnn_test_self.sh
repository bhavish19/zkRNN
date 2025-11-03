#!/usr/bin/env bash
export BENCH_CSV=../bench/zkpot.csv
mkdir -p ../bench
cd build/src

./zkpot 64 128 128 128 1 10
# ./zkpot 256 128 128 128 1 10

# ./zkpot 64 256 512 256 1 10
# ./zkpot 256 256 512 256 1 10

# ./zkpot 64 256 1024 256 1 10
# ./zkpot 256 256 1024 256 1 10

# ./zkpot 64 1024 1024 1024 1 10
# ./zkpot 256 1024 1024 1024 1 10

# ./zkpot 64 2048 2048 2048 1 10
# ./zkpot 256 2048 2048 2048 1 10









