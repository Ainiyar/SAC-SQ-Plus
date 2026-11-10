#!/bin/sh
set -eu
cd "$(dirname "$0")"
g++ -O3 -std=c++17 src/aug_tropical_multiwitness.cpp -o baseline16
g++ -O3 -std=c++17 src/exact_interval_cover_hull.cpp -o exact_interval
g++ -O3 -std=c++17 src/exact_best_scenario_support.cpp -o exact_support
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 run_r399.py
