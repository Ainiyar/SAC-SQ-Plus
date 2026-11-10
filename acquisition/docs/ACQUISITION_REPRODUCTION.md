# RQ3 fitted-cost acquisition reproduction

## Source snapshots and their role

- `acquisition/r333/` is the archived R333 pipeline snapshot used for the three n=256 matched Occupancy, EEG, and SECOM workload pairs.
- `acquisition/r343/` is the updated pipeline snapshot used for the final Occupancy n=1024 SAC arm. Its `proof_relevance.py` implements the cached exact window-partition update used in that run.
- `acquisition/r372/` is the patched certificate-bank-pricing snapshot used for the two additional Occupancy n=256 arms. The portable runner launches the baseline and certificate-bank arms in separate fresh processes and passes the baseline JSON as the new arm's cost calibration.
- `acquisition/common/` contains the shared source dependencies imported by both snapshots, including `safe_archive.py`, the PrefixBank interval arithmetic module, the R326 shared DP/hull primitives, and the source-fit / relevance / closure modules.
- `acquisition/crops_bridge.R` calls the external `crops` R package 1.0.5. The package itself is not bundled.

The only portability edits to the archived R333/R343/R372 source snapshots change an unused Occupancy n=512 cached-input branch from an original-workspace absolute path to the repository-relative `data/processed/frozen_reference/occupancy_768_n512_Xy.npz`. The active RQ3 windows and algorithm code are otherwise the archived source snapshots.

## UCI windows

The exact source preprocessing is called from `r333/full_pipeline.py::preprocess`:

| Dataset | Zero-based window | Rows | Predictors |
|---|---:|---:|---:|
| Occupancy | [1280, 1536) | 256 | 5 original variables + intercept + linear trend |
| EEG Eye State | [1408, 1664) | 256 | 14 EEG channels + intercept + linear trend |
| SECOM | [512, 768) | 256 | 8 selected sensor channels + intercept + linear trend |
| Occupancy | [1536, 2560) | 1024 | 5 original variables + intercept + linear trend |
| Occupancy (R372 pair) | [1536, 1792) | 256 | 5 original variables + intercept + linear trend |

Occupancy rows are stable-sorted by timestamp. EEG preserves the original ARFF row order. SECOM preserves the original source sensor ordering after selecting the first eight non-constant columns with at most 5% missingness; the selected window is median-imputed. For each window, predictor columns are standardized using the population standard deviation (values below `1e-12` are replaced by 1); the intercept and `linspace(-1,1,n)` trend are appended without standardization. The target labels are retained as integers.

Run from `acquisition/`:

```bash
python3 scripts/download_uci_data.py
python3 scripts/prepare_uci_windows.py
```

The generated arrays are written to `data/processed/generated/` and compared against the frozen original X/y arrays in `data/processed/frozen_reference/`. The R372 pair uses the separate [1536,1792) Occupancy window.

## Re-running online acquisition

The three n=256 pairs are reproduced with:

```bash
python3 scripts/run_rq3_n256_pairs.py
```

This executes the original six arms (SAC lookahead and strict CROPS+PELT for each dataset) with the archived R333 configuration: persistent pricing, entropy target priority, debt temperature, sparse4 initialization, source pool 1, incremental lower closure, proposal-fit reuse, no bank cap, and a 400-action budget. The lookahead arm charges round overhead. The script writes each result JSON, log, and a summary under `data/fresh_runs/rq3_n256/`.

The R372 pair is reproduced with:

```bash
python3 scripts/run_r372_pair.py
```

It uses the frozen Occupancy [1536,1792) X/y window, the common R372 run configuration, and the certificate-bank-pricing flag only for the second arm. Both archived outputs and action/proof traces are included.

The updated Occupancy n=1024 SAC arm is reproduced with:

```bash
python3 scripts/run_r343_n1024.py
```

It uses the frozen preprocessed n=1024 X/y window and the R343 source snapshot. It does not rerun the historical strict CROPS+PELT n=1024 reference, whose archived run is retained under `data/archived_runs/r333/` and `data/archived_runs/r343/`.

## CPU interpretation

The R380 stagewise table is the sum of archived online acquisition CPU and the same-host R357/R354 verification replay CPU. It is not a single continuous stopwatch. The R399 paired comparison starts from frozen L/U matrices and does not include upstream GLM fitting. These scopes are deliberately kept separate in the manuscript and must remain separate in any re-analysis.
