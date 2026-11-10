#!/usr/bin/env python3
"""Portable runner for the R343 updated Occupancy n=1024 acquisition arm."""
from pathlib import Path
import importlib.util, json, os, sys
import numpy as np

ACQ = Path(__file__).resolve().parents[1]
REPO = ACQ.parent
COMMON = ACQ / "common"
R343 = ACQ / "r343"
DEFAULT_INPUT = ACQ / "data" / "processed" / "frozen_reference" / "occupancy_1536_n1024_Xy.npz"
DEFAULT_OUT = ACQ / "data" / "fresh_runs" / "r343_n1024" / "r343_window.json"


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    a = p.parse_args()
    if not a.input.exists():
        raise SystemExit(f"Missing X/y input {a.input}; run the UCI download and preprocessing steps first")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    for k in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"]:
        env[k] = "1"
    sys.path.insert(0, str(COMMON))
    sys.path.insert(0, str(R343))
    os.chdir(ACQ)
    with np.load(a.input) as z:
        X = z["X"].copy(); y = z["y"].copy()
    spec = importlib.util.spec_from_file_location("full_pipeline", R343 / "full_pipeline.py")
    fp = importlib.util.module_from_spec(spec); sys.modules["full_pipeline"] = fp; spec.loader.exec_module(fp)
    def frozen_preprocess(start, dataset="occupancy", n=128):
        if (dataset, start, n) != ("occupancy", 1536, 1024):
            raise ValueError(f"Unexpected R343 input request {(dataset, start, n)}")
        return X.copy(), y.copy()
    fp.preprocess = frozen_preprocess
    out = a.out.resolve()
    sys.argv = [str(R343 / "full_pipeline.py"), "--dataset", "occupancy", "--start", "1536", "--n", "1024",
                "--method", "lookahead", "--pricing", "persistent", "--target-priority", "entropy",
                "--entropy-temperature", "debt", "--initial-policy", "sparse4", "--source-pool", "1",
                "--closure", "--closure-kernel", "incremental", "--reuse-proposal-fits", "--bank-cap", "0",
                "--budget", "400", "--charge-round-overhead", "--bank-geometry", "union", "--out", str(out)]
    fp.main()
    result = json.loads(out.read_text())
    result["release_input_path"] = str(a.input)
    out.write_text(json.dumps(result, indent=2))
    print(json.dumps({"status": result.get("status"), "total_cpu": result.get("total_cpu"),
                      "actions": result.get("actions"), "archive_paths": len(result.get("archive_paths", [])),
                      "complete_box_support_certificate": result.get("complete_box_support_certificate")}, indent=2))


if __name__ == "__main__":
    main()
