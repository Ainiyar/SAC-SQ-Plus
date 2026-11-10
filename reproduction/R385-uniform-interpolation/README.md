# R385 m=8 uniform-interpolation separation

Run `python3 reproduce.py`. Both C++ programs are compiled from source. The script recreates the L/U matrices and witness/candidate files, runs R357 with 16 beta cells, applies R354 to the returned paths and all 256 blockwise A/B candidates, and checks selected exact-rational interpolation cases.

The finite interpolation checks are implementation cross-checks. The all-t, all-beta conclusion follows from the proposition proof in the manuscript, not from finite sampling.
