#!/usr/bin/env python3
"""Recreate the exact X/y windows used by the R333 fitted-cost acquisition runs.

This deliberately calls the archived R333 full_pipeline.preprocess implementation,
so feature order, stable time sorting, imputation, scaling and trend construction
are inherited from the original experiment source rather than reimplemented here.
"""
from pathlib import Path
import importlib.util
import os
import sys
import numpy as np

ACQ = Path(__file__).resolve().parents[1]
REPO = ACQ.parent
COMMON = ACQ / "common"
R333 = ACQ / "r333"
OUT = ACQ / "data" / "processed" / "generated"
REF = ACQ / "data" / "processed" / "frozen_reference"

WINDOWS = [
    ("occupancy", 1280, 256, "occupancy_1280_n256_Xy.npz"),
    ("eeg", 1408, 256, "eeg_1408_n256_Xy.npz"),
    ("secom", 512, 256, "secom_512_n256_Xy.npz"),
    ("occupancy", 1536, 1024, "occupancy_1536_n1024_Xy.npz"),
    ("occupancy", 1536, 256, "occupancy_1536_n256_R372_Xy.npz"),
    ("occupancy", 768, 512, "occupancy_768_n512_Xy.npz"),
]


def load_pipeline():
    sys.path.insert(0, str(COMMON))
    sys.path.insert(0, str(R333))
    spec = importlib.util.spec_from_file_location("sac_r333_full_pipeline", R333 / "full_pipeline.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    os.chdir(ACQ)
    for dataset, start, n, filename in WINDOWS:
        expected = {"occupancy": "occupancy/uci_occupancy.zip", "eeg": "eeg/uci_eeg.zip", "secom": "secom/uci_secom.zip"}[dataset]
        if not (ACQ / expected).exists():
            raise SystemExit(f"Missing {expected}; first run: python scripts/download_uci_data.py")
    pipeline = load_pipeline()
    OUT.mkdir(parents=True, exist_ok=True)
    for dataset, start, n, filename in WINDOWS:
        X, y = pipeline.preprocess(start, dataset, n)
        output = OUT / filename
        np.savez_compressed(output, X=X, y=y)
        reference = REF / filename
        comparison = "no frozen reference"
        if reference.exists():
            with np.load(reference) as old:
                same_y = np.array_equal(y, old["y"])
                same_shape = X.shape == old["X"].shape
                max_abs = float(np.max(np.abs(X - old["X"]))) if same_shape else None
                exact_x = bool(same_shape and np.array_equal(X, old["X"]))
            comparison = f"reference_exact_X={exact_x}, labels_equal={same_y}, max_abs_X_diff={max_abs}"
        print(f"{dataset} start={start} n={n} X={X.shape} y={y.shape} -> {output}; {comparison}", flush=True)


if __name__ == "__main__":
    main()
