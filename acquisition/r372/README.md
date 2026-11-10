# R372 certificate-bank-pricing source snapshot

This directory contains the patched R372 `full_pipeline.py`, `proof_lookahead.py`, and `certificate_bank_pricing.py` used for the two additional Occupancy n=256 workload arms. The portable runner is `../scripts/run_r372_pair.py`; it launches the baseline and certificate-bank-pricing arms in separate fresh processes. The second arm consumes the baseline JSON as its cost-model calibration.

The original R372 compatibility issue was limited to an optional `retain_edge_probabilities` keyword passed to `proof_relevance.prepare_context`; the public source snapshot uses the compatible call signature. The exact experiment parameters and outputs are preserved under `../data/archived_runs/r372/`.
