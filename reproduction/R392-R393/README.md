# R392–R393 reproduction bundle

This subdirectory preserves the execution layout for the archive-independent 11-box certification and the analytical constructions used by the SAC-SQ+ manuscript.

## Environment

- Linux x86_64
- Python 3.10 or newer
- NumPy 1.24 or newer
- GCC or Clang with C++17 support (`__int128` is used)

Install the Python dependency and build the two exact-certification programs:

```bash
python3 -m pip install -r requirements.txt
mkdir -p bin output

g++ -O3 -std=c++17 R392/source/exact_interval_cover_hull.cpp -o bin/exact_interval
g++ -O3 -std=c++17 R392/source/exact_best_scenario_support.cpp -o bin/exact_support
```

## Reproduce the 11 frozen cost-box certifications

```bash
python3 R392/source/r392_full_reproduction.py \
  --inputs real_boxes/inputs \
  --out output/real_boxes \
  --cert bin/exact_interval \
  --support bin/exact_support
```

The program creates witnesses from the L/M/U cost scenarios and midpoint exact-cardinality dynamic programming. The archived path identities are consulted only after certification for comparison. Use `--boxes R372_baseline` to run a single small control box.

## Reproduce analytical constructions

```bash
python3 run_singleton_transition.py
python3 run_same_count_stress.py
python3 run_h8_asymmetric_stress.py
```

These scripts write generated inputs and results under `inputs/` and `results/`. They use the exact same R392 certification source and the analytical ground truth encoded in the generators.

## Exhaustive small-instance validation

```bash
python3 sources/validation/random_exhaustive_validation_fast.py
```

This is the 1,000-case positive/negative-penalty exhaustive validation driver. It can take longer than the deterministic examples.

## Data format

The `real_boxes/inputs/*_L.bin` and `*_U.bin` files are signed 64-bit little-endian integer matrices in row-major `(n+1) x (n+1)` layout, with cost scale `Q=10^12`. The associated witness text files are retained as post-run comparison references; they are not used to generate candidates or certify prefixes.

Original UCI observations are not redistributed. Their provenance and dataset DOIs are documented in the repository-level README. The binary matrices are frozen, model-derived certification inputs.
