# SAC-SQ+: Reproducibility Artifact

**Manuscript:** *Collective Exclusion Certificates for Complete Regularization Paths under Interval-Valued Segment Costs*  
**Authors:** Ainiyaer Xilali; Yutian Chen  
**Affiliations:** Fudan University  
**Repository:** `https://github.com/Ainiyar/SAC-SQ-Plus` (planned public repository)

This repository contains the code and data needed to inspect and reproduce the main computational claims of the accompanying R400 manuscript. See `PUBLIC_RELEASE_CONTENTS.md` for the full claim-to-artifact inventory and the deliberate data-redistribution boundaries. It includes the complete-path certification algorithms, exact support filtering, independent witness construction, the earlier fixed-cell comparator, pairwise-only ablations, analytical constructions, exhaustive small-instance tests, the R399 paired benchmark, and the upstream fitted-cost acquisition pipeline used for the RQ3 stagewise comparison.

## What is included

- `reproduction/R392-R393/` — archive-independent witness construction, exact continuous-interval certification, exact closed-support filtering, 11 frozen real fitted-cost boxes, deterministic same-count and singleton-transition constructions, and random validation.
- `reproduction/R399/` — independent witness initialization, both certifiers, the R354 support filter, all 135 paired repetitions (270 executions), raw logs, exact outputs, and summary calculations. It uses the shared 11 frozen L/U boxes under `reproduction/R392-R393/real_boxes/inputs/` to avoid duplicating large binary matrices.
- `reproduction/R390/` — the earlier 15-repeat n=256 same-host CPU comparison, archived 11-box comparison summary, independent-seed parity summaries, random-test results, and SKAB stress logs.
- `reproduction/R381-exhaustive/` — the independent 250-random-DAG exact identity/support comparison against exhaustive enumeration.
- `reproduction/R382-pairwise/` — the six-arm pairwise-only ablation using the same R357 executable with an empty witness set.
- `reproduction/R385-uniform-interpolation/` — the executable m=8 separation between uniform interpolation and the full independent interval-cost box.
- `acquisition/` — the R333 fitted-cost pipeline snapshot, the R343 updated n=1024 snapshot, the R372 certificate-bank-pricing pair, all source dependencies, the official CROPS socket bridge, UCI data preparation scripts, frozen X/y windows, and archived original acquisition outputs.
- `reproduction/R392-R393/figures/` and `figure_data/` — the five figures used by the manuscript and the numerical data underlying the quantitative panels.
- `docs/REPRODUCTION_MAP.md` — maps each major manuscript claim to its implementation, inputs, and recorded outputs.

The repository is deliberately organized around the paper's evidence chain. It does not include obsolete exploratory branches or internal project-management records.

## Quick setup: certification experiments

The core certification experiments require Python 3.10+ (NumPy) and a C++17 compiler (`g++`). From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

### Reproduce all 11 frozen fitted-cost boxes (R392)

```bash
cd reproduction/R392-R393
mkdir -p bin output

g++ -O3 -std=c++17 R392/source/exact_interval_cover_hull.cpp -o bin/exact_interval
g++ -O3 -std=c++17 R392/source/exact_best_scenario_support.cpp -o bin/exact_support

python3 R392/source/r392_full_reproduction.py \
  --inputs real_boxes/inputs \
  --out output/real_boxes \
  --cert bin/exact_interval \
  --support bin/exact_support
```

For a quick single-box run, add `--boxes R372_baseline`.

### Reproduce the 250 exact small-DAG checks (R381)

```bash
python3 reproduction/R381-exhaustive/exhaustive_250.py
```

This independently enumerates all legal paths, calculates exact `Fraction` supports, then compares the result with R357 candidate enumeration and R354 exact support filtering. It reports the 646 possible path-instance pairs, fractional endpoints, and singleton supports.

### Reproduce the six-arm pairwise-only ablation (R382)

```bash
python3 reproduction/R382-pairwise/run_pairwise_ablation.py
```

The pairwise-only arm uses the R357 binary with an empty witness file, retaining same-node strict robust dominance and 16 closed beta cells. Three SAC arms reach the 2,000,000-label cap; the three strict-reference arms complete and are checked against the R392 exact supports. This is a controlled ablation of the shared-beta search skeleton, not a reproduction of the original Tropical Dijkstra executable.

### Reproduce the uniform-interpolation separation (R385)

```bash
cd reproduction/R385-uniform-interpolation
python3 reproduce.py
```

This rebuilds both C++ programs from source and verifies the 256 full-box possible paths, their exact closed supports, and the unique all-A path under uniform interpolation.

### Reproduce the 1,000 random small-DAG candidate-inclusion tests (R390)

Build the exact interval certifier as above, then run:

```bash
cd reproduction/R392-R393
python3 sources/validation/random_exhaustive_validation_fast.py
```

This test checks that every independently enumerated box-possible path is retained in the exact certifier's safe candidate superset. It is separate from the 250-case test, which also checks exact final supports.

### Reproduce the earlier R390 15-repeat comparison

```bash
python3 reproduction/R390/repeated_same_host_benchmark.py
```

This runs 15 fresh-process repetitions on the six n=256 arms and records enumerator plus R354 support CPU. The original single-run 11-box comparison is preserved as `reproduction/R390/evidence/11_box_summary.json`; the newer R399 paired benchmark is the preferred full 11-box rerun. Archived R390 summary data are included; new CPU values will depend on the host.

### Reproduce the R399 paired benchmark

```bash
cd reproduction/R399
bash rebuild_and_run.sh
```

This performs 135 paired repetitions (270 certifier executions) across 11 workload arms. It may take several minutes. The original logs and outputs are already included; rerunning is optional for inspecting the published results.

## Reproduce the fitted-cost acquisition experiments (RQ3)

The manuscript's stagewise acquisition comparison uses the original R333 pipeline for three matched n=256 workload pairs and the R343 updated pipeline for the Occupancy n=1024 SAC arm. The strict n=1024 CROPS+PELT reference is retained as an archived original run because rerunning it requires a long sequence of direct segment fits.

The exact acquisition source snapshots and dependencies are included under `acquisition/`. The original UCI observations are not duplicated in this repository; the preprocessing windows, derived fitted-cost boxes, and original run outputs are included. Download the source datasets and recreate the windows with:

```bash
cd acquisition
python3 scripts/download_uci_data.py
python3 scripts/prepare_uci_windows.py
```

The preparation script calls the archived R333 `full_pipeline.preprocess` implementation and compares the generated X/y windows with the frozen references.

For a fresh reproduction of the three n=256 matched SAC/CROPS+PELT pairs:

```bash
python3 scripts/run_rq3_n256_pairs.py
```

This launches six full acquisition runs and can take substantial CPU time. It uses the original R333 source snapshot and the original experiment configuration.

For the two R372 Occupancy n=256 workload arms (persistent baseline and certificate-bank pricing):

```bash
python3 scripts/run_r372_pair.py
```

The baseline and certificate-bank arms are launched in separate fresh processes; the second arm uses the baseline run as its cost-model calibration.

For the updated Occupancy n=1024 SAC arm:

```bash
python3 scripts/run_r343_n1024.py
```

This requires R 4.1.2, `crops` 1.0.5, and `jsonlite`. The bridge is `acquisition/crops_bridge.R`. The RQ3 stagewise CPU figures are historical same-ECS measurements; reproducing them on different hardware is not expected to produce identical wall-clock or CPU values. The archived original runs and timing records are included under `acquisition/data/archived_runs/`.

See `acquisition/docs/ACQUISITION_REPRODUCTION.md` for the exact run contract and the distinction between online acquisition CPU, fixed-cell certification replay, and R399's frozen-box-to-support CPU.

## Data and interpretation boundaries

- The 11 R399 workload arms are configurations across Occupancy, EEG, SECOM, and R372 controls; they are not 11 independent datasets.
- R399 begins from frozen L/U cost boxes and excludes upstream GLM fitting. The RQ3 stagewise comparison starts earlier and includes fitted-cost acquisition plus the original fixed-cell certification replay. These CPU measures must not be combined.
- The 2.633x R399 certification-plus-support speedup is an internal controlled comparison between the 16-cell and exact-interval certifiers under the same independent witness pools and output contract. It is not an official third-party Tropical Dijkstra benchmark.
- The compact SECOM strict arm is slightly slower under exact-interval certification; this result is retained.
- SKAB n=768 is an unresolved stress case at the 16-million-label budget. Its upstream dataset is not redistributed here; see `docs/SKAB_REPRODUCTION.md` for provenance and the original source terms.

## Licenses

- Authors' source code and scripts: MIT License (`LICENSE`).
- Authors' original and derived data, frozen cost boxes, experiment tables, figures, and documentation: CC BY 4.0 (`DATA_LICENSE.md`), with dataset attribution in `DATA_SOURCES.md`.
- The original UCI datasets are not bundled; download them from UCI under their published CC BY 4.0 terms.
- Third-party data and software retain their original licenses. The external R `crops` package is not bundled and remains under its own GPL license.

## Citation

Please cite the accompanying manuscript and this exact software/data release. `CITATION.cff` contains the software citation metadata. The Zenodo DOI will be added to the citation metadata and Data and Code Availability statement after the GitHub Release has been archived. See `docs/ZENODO_PUBLICATION_GUIDE.md` for the mixed-license release procedure.
