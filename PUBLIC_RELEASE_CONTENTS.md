# Public release contents: R400 evidence coverage

This release candidate is intentionally broader than the last few certification experiments. It includes the older source-acquisition and baseline components that generate the fitted-cost inputs used throughout R400.

## Complete evidence chain included

1. **Upstream fitted-cost acquisition**
   - R333 original pipeline source and its complete imported dependency closure, including the R326 shared DP/hull and R325 interval arithmetic modules.
   - R343 updated source snapshot used for the Occupancy n=1024 SAC arm.
   - R372 certificate-bank-pricing source snapshot used for two additional Occupancy n=256 workload arms.
   - CROPS 1.0.5 R bridge, Python/R environment records, UCI download and exact window-preparation scripts.
   - Frozen X/y windows, eight original R333 run outputs, the full R343 JSON/NPZ output and strict reference, and both R372 paired outputs with traces.

2. **Core path-certification algorithms**
   - R357 fixed-cell multi-witness enumerator.
   - R390 exact continuous witness-interval certifier.
   - R354 exact closed-support filter.
   - R392 independent L/M/U and exact-k witness construction, including the 11 frozen real cost boxes.

3. **Correctness and structural experiments**
   - R381 independent exact exhaustive validation (250 randomized DAGs; 646 possible path-instance pairs; 27 fractional-endpoint instances; 76 singleton supports).
   - R385 m=8 uniform-interpolation separation (256 full-box possible identities).
   - R393 same-count and isolated-transition constructions.
   - R390 1,000 random candidate-inclusion tests and the original 11-box / 15-repeat summaries.

4. **Comparative experiments**
   - R382 six-arm pairwise-only ablation, with original archived outputs and a portable rerun driver.
   - R399 independent-witness paired benchmark: 11 workload arms, 135 pairs, 270 executions, raw logs, exact outputs, and summary calculations.

5. **Paper outputs**
   - The five exact R395 PNG figures used in R400.
   - Numerical backing tables and stagewise CPU ledgers.
   - Manuscript-to-artifact mapping and scope notes.

## Data included and not included

- Included: authors' processed UCI X/y windows, all 11 frozen L/U boxes, derived experiment outputs, support records, raw benchmark logs, and numerical figure data.
- Not duplicated: original UCI raw archives. They are publicly available under CC BY 4.0 and are fetched by `acquisition/scripts/download_uci_data.py`.
- Not redistributed: SKAB raw data and SKAB-specific derived matrices, because the Kaggle data record carries AGPL-3.0 terms. The R390 stress logs and source citation are retained; see `docs/SKAB_REPRODUCTION.md`.
- Not included: the rejected manuscript, internal editorial analysis, the internal research-state document, and unrelated exploratory R333 branches that are not dependencies of the R400 results.

The 11 R399 arms are workload configurations, not 11 independent datasets. R399's CPU starts from frozen L/U boxes and excludes upstream fitting; RQ3's stagewise CPU includes acquisition and must remain a separate comparison.
