# R399 — Same-Box Independent Witness Paired Certification

Run `bash rebuild_and_run.sh` on Linux with Python3, NumPy and g++.

Source inputs: `../R392-R393/real_boxes/inputs/*_L.bin` and `*_U.bin` (the shared frozen completed model-fit interval boxes). Independent witnesses are generated from these cost boxes using R392 `r392/source/r392_full_reproduction.py` functions, with no archive candidate input. `r392/results` supports are consulted only as post-run reference values.

Run logs (`runs/<case>/...`) include each `COMPLETE` enumerator stdout, exact R354 support file and independently regenerated witnesses. `runs/R399_RAW_RUNS.json` contains every timing observation, `runs/R399_SUMMARY.json` the per-arm medians and aggregates; `R399_PAIRED_SUMMARY.csv` is the manuscript table source.

Scientific interpretation and limitations are in `R399_SCIENTIFIC_REPORT.md`.
