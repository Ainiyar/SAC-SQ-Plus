# R381 independent exhaustive validation

The portable runner `exhaustive_250.py` independently enumerates every legal segmentation path for 250 randomized small DAGs (n in {7,8,9,10}, h in {1,2,3}), constructs integer lower/upper edge-cost boxes, and computes exact rational support endpoints using Python `Fraction`. It then checks candidate inclusion from the R357 enumerator and exact identity/support equality after the R354 filter.

Recorded result: 250/250 instances passed; 646 path-instance pairs were box-possible; 27 had at least one non-integer rational support endpoint; 76 had singleton supports. The per-instance summary is `results/summary_250.json`. These are implementation cross-checks; the manuscript's analytic theorems provide the proof.
