# R382 pairwise-only ablation

`run_pairwise_ablation.py` rebuilds the R357 binary and runs six n=256 arms with an empty witness file and 16 closed beta cells. The three SAC arms reach the 2,000,000-label cap; the strict-reference arms complete and their R354 exact supports are compared with the R392 reference outputs.

This is a controlled ablation of the same search skeleton with multi-witness pruning disabled. It is not the original Tropical Dijkstra implementation and should not be described as an official third-party baseline.
