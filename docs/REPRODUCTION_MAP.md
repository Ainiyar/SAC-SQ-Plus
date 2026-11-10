# Manuscript-to-artifact map

This map links the R400 paper's principal claims to the source, inputs, and archived outputs in this repository.

| Manuscript component / claim | Implementation | Inputs and archived evidence |
|---|---|---|
| Best-scenario exact closed-support characterization | `reproduction/R392-R393/R392/source/exact_best_scenario_support.cpp` | 11 frozen L/U boxes; exact-support outputs under `reproduction/R392-R393/reference_results/R392/` |
| Collective multi-witness prefix exclusion on fixed beta cells | `reproduction/R399/src/aug_tropical_multiwitness.cpp` | R357 source, same-binary pairwise ablation under `reproduction/R382-pairwise/` |
| Exact continuous witness-interval cover | `reproduction/R392-R393/R392/source/exact_interval_cover_hull.cpp` | R390 correctness runs, 1,000 random small-DAG test driver, and 11-box complete results |
| Archive-independent candidate/witness generation | `reproduction/R392-R393/R392/source/r392_full_reproduction.py` | L/M/U 17-point seeds, midpoint exact-k DP, 11-box result package |
| 250 independent exact identity/support checks | `reproduction/R381-exhaustive/exhaustive_250.py` | Independent path enumeration and Fraction support oracle |
| Same-count exponential identity family and isolated transition identity | `reproduction/R392-R393/` scripts and results | Deterministic h=8 constructions, exact supports, and singleton transition at beta=0.665 |
| Uniform-interpolation separation | `reproduction/R385-uniform-interpolation/` | m=8 L/U matrices, 256 candidates, R357/R354 logs, exact interpolation checks |
| Pairwise-only ablation | `reproduction/R382-pairwise/run_pairwise_ablation.py` | Six n=256 L/U arms from R399 inputs; exact-support comparison against R392 |
| Earlier 15-repeat and 11-box CPU comparisons | `reproduction/R390/` | Original 15-repeat n=256 summary, 11-box result summary, and portable rerun scripts |
| Independent-witness paired comparison (R399) | `reproduction/R399/run_r399.py` | 11 boxes, 135 paired repetitions / 270 executions, raw logs and exact outputs |
| Fitted-cost source acquisition and PrefixBank construction | `acquisition/r333/`, `acquisition/r343/`, `acquisition/r372/`, `acquisition/common/` | Original UCI X/y windows, eight R333 acquisition outputs, R343 full JSON/NPZ and strict reference, R372 baseline/certificate-bank pair |
| RQ3 stagewise CPU comparison | `acquisition/data/archived_runs/`, `reproduction/R392-R393/figure_data/` | R380 stagewise table, R333/R343 source runs, original CPU ledgers |
| Manuscript figures | `reproduction/R392-R393/figures/` | Exact R395 PNGs used by R400; quantitative source tables in sibling `figure_data/` |

## Important scope distinctions

1. The 11 R399 arms are workload configurations, not 11 independent datasets.
2. R399 begins from frozen L/U boxes and excludes upstream GLM fitting. RQ3 stagewise acquisition starts earlier and includes source fitting, PrefixBank construction/query, and the original certification replay.
3. The R399 comparator is the paper's controlled 16-cell versus exact-interval implementation comparison under the same witness pool and output contract. It is not an official third-party Tropical Dijkstra implementation.
4. The exact-box guarantee concerns the Cartesian product of certified edge-cost intervals. It does not assert that every combination in this conservative outer box is realizable by the shared upstream model parameters.
5. SKAB n=768 remains an unresolved stress case at the reported label budget; its upstream dataset is not redistributed here.
