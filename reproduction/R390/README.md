# R390 exact-interval benchmark and validation assets

This folder preserves the R390 source drivers and archived evidence that feed the R400 quantitative figures. `evidence/n256_15rep_pair_summary.json` is the original six-arm, 15-fresh-process same-host comparison. `evidence/11_box_summary.json` is the original single-run fixed-cell versus exact-interval comparison over all 11 frozen boxes. `evidence/random_1000_summary.json` is the 1,000-case safe-candidate-inclusion experiment; the executable driver is under `reproduction/R392-R393/sources/validation/`. The SKAB label-budget logs are retained under `evidence/skab/`; the SKAB data and derived matrices are not redistributed.

The portable `repeated_same_host_benchmark.py` reruns the six n=256 15-repeat comparison. The archived `11_box_summary.json` is the single-run reference used in the manuscript; the newer R399 paired benchmark is the preferred full 11-box rerun. CPU measurements will vary by host.
