# R372 configuration and compatibility

Both arms use the same Occupancy window [1536,1792), n=256, h=8, beta=[0.1,1.2], and a 400-action budget. Shared settings: lookahead, sparse4 initialization, persistent pricing, witness source priority, entropy target priority, debt temperature, source pool 1, uniform obligation weight, bank cap 0, incremental closure, tiling closure geometry, full bank geometry, target source policy, proposal-fit reuse, and round-overhead charging.

The second arm additionally enables certificate-bank pricing and uses the baseline JSON as its calibration source. The public runner `acquisition/scripts/run_r372_pair.py` launches the two arms in separate processes.
