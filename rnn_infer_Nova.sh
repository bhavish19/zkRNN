#!/usr/bin/env bash
set -euo pipefail

export BENCH_CSV=../bench/zkpot.csv
mkdir -p ../bench
cd build/src

T=${1:-64}
N=${2:-256}
H=${3:-512}
O=${4:-256}
MODEL_PATH=${5:-}
INPUT_PATH=${6:-}

cmd=("./zkpot" "$T" "$N" "$H" "$O" 1 1 "--mode=infer")

if [[ -n "${MODEL_PATH}" ]]; then
  cmd+=("--model=${MODEL_PATH}")
fi

if [[ -n "${INPUT_PATH}" ]]; then
  cmd+=("--inputs=${INPUT_PATH}")
fi

"${cmd[@]}"

