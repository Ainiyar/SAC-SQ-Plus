#!/usr/bin/env python3
"""Run one R372 arm in a fresh process using the frozen Occupancy X/y window."""
from pathlib import Path
import importlib.util, os, sys
import numpy as np

ACQ=Path(__file__).resolve().parents[1]
COMMON=ACQ/"common"; R372=ACQ/"r372"
INPUT=ACQ/"data"/"processed"/"frozen_reference"/"occupancy_1536_n256_R372_Xy.npz"


def main():
    import argparse
    p=argparse.ArgumentParser();p.add_argument("--arm",choices=["baseline","certificate_bank_pricing"],required=True)
    p.add_argument("--out",type=Path,required=True);p.add_argument("--cost-model-from",type=Path)
    a=p.parse_args()
    if not INPUT.exists(): raise SystemExit(f"Missing {INPUT}")
    for k in ["OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"]: os.environ[k]="1"
    sys.path.insert(0,str(COMMON));sys.path.insert(0,str(R372));os.chdir(ACQ)
    with np.load(INPUT) as z: X=z["X"].copy();y=z["y"].copy()
    spec=importlib.util.spec_from_file_location("full_pipeline",R372/"full_pipeline.py")
    fp=importlib.util.module_from_spec(spec);sys.modules["full_pipeline"]=fp;spec.loader.exec_module(fp)
    def frozen_preprocess(start,dataset="occupancy",n=128):
        if (dataset,start,n)!=("occupancy",1536,256): raise ValueError((dataset,start,n))
        return X.copy(),y.copy()
    fp.preprocess=frozen_preprocess
    args=[str(R372/"full_pipeline.py"),"--method","lookahead","--dataset","occupancy","--start","1536","--n","256",
          "--budget","400","--initial-policy","sparse4","--pricing","persistent","--source-priority","witness",
          "--target-priority","entropy","--entropy-temperature","debt","--source-pool","1","--obligation-weight","uniform",
          "--bank-cap","0","--closure","--closure-kernel","incremental","--closure-geometry","tiling",
          "--bank-geometry","full","--source-policy","target","--reuse-proposal-fits","--charge-round-overhead","--out",str(a.out.resolve())]
    if a.arm=="certificate_bank_pricing":
        if not a.cost_model_from or not a.cost_model_from.exists(): raise SystemExit("certificate-bank arm requires baseline --cost-model-from")
        args += ["--certificate-bank-pricing","--cost-model-from",str(a.cost_model_from.resolve())]
    sys.argv=args;fp.main()

if __name__=="__main__":main()
