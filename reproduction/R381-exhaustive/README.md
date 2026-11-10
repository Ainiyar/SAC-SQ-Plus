# R381 independent exhaustive validation

Run `python3 exhaustive_250.py` from any working directory. The script compiles the R357 multi-witness enumerator and R354 exact-support filter from the R399 source tree, generates 250 small random DAG instances, exhaustively enumerates all legal paths, computes exact Fraction supports, and checks exact identity/support equality.

The original R381 validation report is retained. The new runner replaces the original workspace-specific binary paths with repository-relative source builds; the mathematical oracle and random instance family are unchanged in intent.
