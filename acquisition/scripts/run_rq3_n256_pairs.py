#!/usr/bin/env python3
"""Re-run the three original n=256 matched SAC-vs-CROPS+PELT acquisitions.

This is an expensive optional reproduction: six full acquisition runs. It uses
R333 source snapshots and the original protocol from run_whole_round_holdout.py.
"""
from pathlib import Path
import json, os, resource, subprocess, sys, time

ACQ = Path(__file__).resolve().parents[1]
REPO = ACQ.parent
R333 = ACQ / "r333"
COMMON = ACQ / "common"
OUT = ACQ / "data" / "fresh_runs" / "rq3_n256"
CASES = [("occupancy", 1280), ("eeg", 1408), ("secom", 512)]


def child_cpu():
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    return r.ru_utime + r.ru_stime


def main():
    for name in ["occupancy/uci_occupancy.zip", "eeg/uci_eeg.zip", "secom/uci_secom.zip"]:
        if not (ACQ / name).exists():
            raise SystemExit(f"Missing {name}; run python scripts/download_uci_data.py first")
    OUT.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update({"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "1"})
    env["PYTHONPATH"] = os.pathsep.join([str(R333), str(COMMON), env.get("PYTHONPATH", "")]).rstrip(os.pathsep)
    rows = []
    os.chdir(ACQ)
    for dataset, start in CASES:
        for method in ["lookahead", "crops_pelt"]:
            stem = f"{dataset}_{start}_{method}"
            output = OUT / f"{stem}.json"
            command = [sys.executable, str(R333 / "full_pipeline.py"),
                       "--dataset", dataset, "--start", str(start), "--n", "256",
                       "--method", method, "--pricing", "persistent", "--target-priority", "entropy",
                       "--entropy-temperature", "debt", "--initial-policy", "sparse4", "--source-pool", "1",
                       "--closure", "--closure-kernel", "incremental", "--reuse-proposal-fits",
                       "--bank-cap", "0", "--budget", "400", "--out", str(output)]
            if method == "lookahead":
                command.append("--charge-round-overhead")
            before = child_cpu(); wall0 = time.perf_counter()
            log = OUT / f"{stem}.log"
            with log.open("w", encoding="utf-8") as f:
                p = subprocess.run(command, cwd=ACQ, env=env, stdout=f, stderr=subprocess.STDOUT)
            cpu = child_cpu() - before; wall = time.perf_counter() - wall0
            if p.returncode:
                raise RuntimeError(f"{stem} failed; inspect {log}")
            result = json.loads(output.read_text())
            row = {"dataset": dataset, "start": start, "n": 256, "method": method,
                   "returncode": p.returncode, "status": result.get("status"),
                   "pipeline_cpu": result.get("total_cpu"), "child_process_cpu": cpu,
                   "wall_seconds": wall, "actions": result.get("actions"),
                   "archive_paths": len(result.get("archive_paths", [])),
                   "result_json": output.name, "log": log.name}
            rows.append(row)
            (OUT / "summary.json").write_text(json.dumps({"runs": rows}, indent=2))
            print(json.dumps(row), flush=True)
    print(f"Completed {len(rows)} n=256 acquisition runs. Results: {OUT}", flush=True)


if __name__ == "__main__":
    main()
