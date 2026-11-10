# Fitted-cost acquisition source

This directory contains the source dependencies needed to inspect and rerun the RQ3 fitted-cost acquisition pipeline.

- `common/`: shared R333 source modules and the R326/R325 numerical primitives imported by `safe_archive.py`.
- `r333/`: the original pipeline snapshot used for the three n=256 matched acquisition pairs.
- `r343/`: the updated pipeline snapshot used for the Occupancy n=1024 SAC arm.
- `r372/`: the certificate-bank-pricing snapshot used for the two additional Occupancy n=256 control arms.
- `crops_bridge.R`: socket bridge to the external R `crops` package.
- `scripts/`: UCI download/preprocessing and portable experiment runners, including the R372 two-arm pair.
- `data/processed/frozen_reference/`: frozen X/y windows extracted from the original R333 package; these are comparison references, not raw observations.
- `data/archived_runs/`: original acquisition JSON/NPZ/log/time outputs, including the R343 updated arm and the historical strict n=1024 reference.

## Environment

The recorded environment was Python 3.10.12, NumPy 1.21.5, SciPy 1.8.0, Numba 0.55.1, pandas 1.5.3, R 4.1.2, CROPS 1.0.5, and jsonlite. Install the Python environment from the repository root with `pip install -r requirements-acquisition.txt`. Install the R packages with:

```r
install.packages("jsonlite", repos="https://cloud.r-project.org")
install.packages("https://cran.r-project.org/src/contrib/crops_1.0.5.tar.gz", repos=NULL, type="source")
```

Set `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and `NUMBA_NUM_THREADS=1` for the historical single-threaded comparison protocol.
