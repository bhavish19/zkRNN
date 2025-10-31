#!/usr/bin/env bash
export BENCH_CSV=../bench/zkpot.csv
mkdir -p ../bench
cd build/src


./zkpot 16 256 512 256 1 10
./zkpot 32 256 512 256 1 10
./zkpot 64 256 512 256 1 10
./zkpot 128 256 512 256 1 10
./zkpot 256 256 512 256 1 10
./zkpot 512 256 512 256 1 10

