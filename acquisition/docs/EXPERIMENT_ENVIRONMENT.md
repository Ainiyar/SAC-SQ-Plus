# Runtime environment and measurement

ECS: AMD EPYC under KVM,8 virtual CPUs(4 cores,2 threads/core),15321MiB RAM,
Linux5.15.0-191-generic. Experiments run sequentially with
OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1.

Python3.10.12; NumPy1.21.5, SciPy1.8, Numba0.55.1, pandas1.5.3;
R4.1.2 with official crops1.0.5 and jsonlite. Strict arithmetic is70-digit
directed Decimal with outward integer gridQ=1e12.

Numerical pipeline CPU begins at raw-data preprocessing after common data/graph
and selected pricing-kernel warmup. It includes all initial ordinary proposals,
directed evidence, source fits, bank construction/queries, every guidance and
closure operation, archive/support management and final phase proof, plus the
official R child CPU. Fresh-process CPU measured by GNU time additionally
includes Python/R startup, imports, warmup and process-side output work.
Independent processes reuse the Numba code cache prepared during shared setup.
This fresh-process metric starts a new interpreter with that compiled-code
cache available; the numerical pipeline boundary follows common kernel warmup.

Each cohort has its own frozen protocol and timing boundary. Repeated studies
report medians of within-pair ratios; single prospective workload pairs report
their observed full costs. Raw JSON/gz, NPZ, R bridge traces and GNU time files
are retained. The independent graph implementation checks the final edge-cost
intervals using a different count recurrence, forward trie-product and exact
Fraction support intersections.
