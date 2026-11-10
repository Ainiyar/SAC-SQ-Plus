#!/usr/bin/env python3
"""Reproduce the six-arm R382 pairwise-only vs multi-witness ablation.

The pairwise-only arm is the same R357 binary with an empty witness file and
16 closed beta cells. It is a controlled ablation of the shared-beta search
skeleton, not a re-execution of the original Tropical Dijkstra implementation.
"""
from pathlib import Path
import csv, json, resource, subprocess, tempfile, time

ROOT = Path(__file__).resolve().parents[2]
R399 = ROOT / "reproduction" / "R399"
INPUT = ROOT / "reproduction" / "R392-R393" / "real_boxes" / "inputs"
R392 = ROOT / "reproduction" / "R392-R393" / "reference_results" / "R392"
SRC = R399 / "src"
OUT = Path(__file__).resolve().parent / "results"
B0, B1 = 100_000_000_000, 1_200_000_000_000
ARMS = [
    "whole_round_holdout_occupancy_1280_lookahead",
    "whole_round_holdout_occupancy_1280_crops_pelt",
    "whole_round_holdout_eeg_1408_lookahead",
    "whole_round_holdout_eeg_1408_crops_pelt",
    "whole_round_holdout_secom_512_lookahead",
    "whole_round_holdout_secom_512_crops_pelt",
]


def child_cpu():
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    return r.ru_utime + r.ru_stime


def parse_paths(text):
    return [line[5:].strip() for line in text.splitlines() if line.startswith("PATH ")]


def parse_exact_support(log_text, support_text):
    paths = parse_paths(log_text); out = {}
    for line in support_text.splitlines():
        if not line.startswith("INDEX "): continue
        z = line.split()
        if z[3] == "1": out[paths[int(z[1])]] = (z[5], z[7], z[9])
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="sac_r382_") as td:
        td = Path(td); r357 = td / "r357"; r354 = td / "r354"; empty = td / "empty_witnesses.txt"
        empty.write_text("")
        subprocess.run(["g++", "-O3", "-std=c++17", str(SRC / "aug_tropical_multiwitness.cpp"), "-o", str(r357)], check=True)
        subprocess.run(["g++", "-O3", "-std=c++17", str(SRC / "exact_best_scenario_support.cpp"), "-o", str(r354)], check=True)
        rows = []
        for arm in ARMS:
            L = INPUT / f"{arm}_L.bin"; U = INPUT / f"{arm}_U.bin"
            if not L.exists() or not U.exists(): raise FileNotFoundError(arm)
            stem = OUT / arm
            cmd = [str(r357), "256", "8", str(B0), str(B1), str(L), str(U), "2000000", "5", str(empty), "16"]
            before = child_cpu(); wall0 = time.perf_counter()
            p = subprocess.run(cmd, capture_output=True, text=True)
            cpu = child_cpu() - before; wall = time.perf_counter() - wall0
            (stem.with_suffix(".log")).write_text(p.stdout + p.stderr)
            status = next((x.split()[1] for x in p.stdout.splitlines() if x.startswith("status ")), "ERROR")
            created_line = next((x for x in p.stdout.splitlines() if x.startswith("pops ")), "")
            created = int(created_line.split()[5]) if created_line else -1
            paths = parse_paths(p.stdout)
            row = {"arm": arm, "status": status, "created_labels": created, "sink_candidates": len(paths),
                   "process_cpu_seconds": cpu, "wall_seconds": wall, "returncode": p.returncode}
            if status == "COMPLETE" and p.returncode == 0:
                support_file = stem.with_suffix(".support.txt")
                q = subprocess.run([str(r354), "256", "8", str(B0), str(B1), str(L), str(U), str(stem.with_suffix(".log")), str(support_file)], capture_output=True, text=True)
                if q.returncode: raise RuntimeError(f"R354 filter failed for {arm}: {q.stderr}")
                actual = parse_exact_support(p.stdout, support_file.read_text())
                refdir = R392 / arm
                reference = parse_exact_support((refdir / "certifier.log").read_text(), (refdir / "exact_support.txt").read_text())
                row["exact_possible_paths"] = len(actual)
                row["matches_R392_exact_support"] = actual == reference
                if actual != reference: raise RuntimeError(f"Exact support mismatch against R392 for {arm}")
            else:
                row["exact_possible_paths"] = None
                row["matches_R392_exact_support"] = None
            rows.append(row)
            (OUT / "summary.json").write_text(json.dumps({"protocol": {"h": 8, "beta": [B0, B1], "cells": 16,
                "max_created": 2_000_000, "time_limit_seconds": 5, "witnesses": "empty"}, "arms": rows}, indent=2))
            print(json.dumps(row), flush=True)
        with (OUT / "summary.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
        print(f"Completed {len(rows)} pairwise-only arms. Summary: {OUT / 'summary.json'}")


if __name__ == "__main__":
    main()
