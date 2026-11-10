# R385 m=8 uniform-interpolation separation: executable validation

## Construction

Eight blocks of length `3h` (`h=1`) each contain two two-edge routes. Route A edges have intervals `[0,1]`; route B edges have `[0.5,1.5]`; all other legal edges have fixed cost `H=1000`. The shared penalty is `beta in [0,1.2]`, and costs/penalties are represented at scale `Q=10^12`. The sufficient high-edge condition `H > 2m(beta_1+epsilon)` holds.

## R357/R354 result

The actual R357 executable returned `COMPLETE`, created 1,021 labels, pruned 212 prefixes using three feasible witnesses and returned exactly the 256 blockwise A/B paths. The R354 exact support executable processed all 256 candidates and marked all 256 possible, each with exact closed support `[0,1.2]`. The returned R357 path set was compared with an independently generated set of all `2^8` route choices.

## Uniform interpolation check

An independent exact-rational DAG dynamic program was run at `t in {0, 1/2, 1}` and `beta in {0, 3/5, 6/5}`. All nine checks returned the unique all-A path. The proposition in the manuscript establishes the same result for every real `t in [0,1]` and every `beta in [0,1.2]`; the finite checks are implementation cross-checks, not a replacement for the proof.

See `R357.log`, `R354_R357.txt`, `R354_all256.txt`, `L.bin`, `U.bin`, `witnesses.txt`, and `all_256_candidates.log` for the complete raw output and inputs. `verify_results.py` reruns the output comparison and fixed-interpolation checks; `build_and_run.py` rebuilds R354 from the archived source and reruns both algorithms (its binary paths are local to the original workspace).
