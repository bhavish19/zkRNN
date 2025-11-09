#!/usr/bin/env bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}" 2>/dev/null)" && pwd)"
REPO_ROOT="$SCRIPT_DIR"
cd "$REPO_ROOT/build/src" || exit 1

mkdir -p "$REPO_ROOT/bench"
export BENCH_CSV="$REPO_ROOT/bench/zkpot.csv"
export BENCH_DUMP="$REPO_ROOT/bench/witness_dump.json"

./zkpot 8 32 32 32 1 3
# ./zkpot 64 128 128 128 1 10
# ./zkpot 256 128 128 128 1 10

# ./zkpot 64 256 512 256 1 10
# ./zkpot 256 256 512 256 1 10

# ./zkpot 64 256 1024 256 1 10
# ./zkpot 256 256 1024 256 1 10

# ./zkpot 64 1024 1024 1024 1 10
# ./zkpot 256 1024 1024 1024 1 10

# ./zkpot 64 2048 2048 2048 1 10
# ./zkpot 256 2048 2048 2048 1 10









