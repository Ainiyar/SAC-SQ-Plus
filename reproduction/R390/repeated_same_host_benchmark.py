#!/usr/bin/env python3
"""Portable reproduction of the six-arm R390 15-repeat n=256 CPU comparison."""
from pathlib import Path
import json, resource, statistics, subprocess, tempfile, time

ROOT=Path(__file__).resolve().parents[2]
INPUT=ROOT/"reproduction"/"R392-R393"/"real_boxes"/"inputs"
SRC=ROOT/"reproduction"/"R399"/"src"
OUT=Path(__file__).resolve().parent/"results"/"n256_15rep"
B0,B1=100_000_000_000,1_200_000_000_000
CASES=[
"whole_round_holdout_occupancy_1280_lookahead",
"whole_round_holdout_occupancy_1280_crops_pelt",
"whole_round_holdout_eeg_1408_lookahead",
"whole_round_holdout_eeg_1408_crops_pelt",
"whole_round_holdout_secom_512_lookahead",
"whole_round_holdout_secom_512_crops_pelt",
]

def child_cpu(cmd):
    before=resource.getrusage(resource.RUSAGE_CHILDREN)
    p=subprocess.run(cmd,capture_output=True,text=True)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu=(after.ru_utime-before.ru_utime)+(after.ru_stime-before.ru_stime)
    if p.returncode: raise RuntimeError((cmd,p.stderr,p.stdout[-1000:]))
    return cpu,p.stdout

def main():
    OUT.mkdir(parents=True,exist_ok=True); rows=[]; start=time.time()
    with tempfile.TemporaryDirectory(prefix="sac_r390_") as td:
        td=Path(td); old=td/"grid16"; new=td/"exact_interval"; sup=td/"exact_support"
        subprocess.run(["g++","-O3","-std=c++17",str(SRC/"aug_tropical_multiwitness.cpp"),"-o",str(old)],check=True)
        subprocess.run(["g++","-O3","-std=c++17",str(ROOT/"reproduction/R392-R393/R392/source/exact_interval_cover_hull.cpp"),"-o",str(new)],check=True)
        subprocess.run(["g++","-O3","-std=c++17",str(SRC/"exact_best_scenario_support.cpp"),"-o",str(sup)],check=True)
        for name in CASES:
            L=INPUT/f"{name}_L.bin"; U=INPUT/f"{name}_U.bin"; W=INPUT/f"{name}_archive_witnesses.txt"
            old_search=[];new_search=[];old_support=[];new_support=[];old_total=[];new_total=[]
            for rep in range(15):
                methods=["old","new"] if rep%2==0 else ["new","old"]; local={}
                for method in methods:
                    exe=old if method=="old" else new
                    args=[str(exe),"256","8",str(B0),str(B1),str(L),str(U),"2000000","120",str(W)]
                    if method=="old": args.append("16")
                    cpu,stdout=child_cpu(args); local[method]=(cpu,stdout)
                for method in ["old","new"]:
                    scpu,stdout=local[method]
                    log=OUT/f"{name}_{method}_{rep:02d}.log"; sf=OUT/f"{name}_{method}_{rep:02d}_support.txt"
                    log.write_text(stdout)
                    cpu2,_=child_cpu([str(sup),"256","8",str(B0),str(B1),str(L),str(U),str(log),str(sf)])
                    if method=="old": old_search.append(scpu);old_support.append(cpu2);old_total.append(scpu+cpu2)
                    else: new_search.append(scpu);new_support.append(cpu2);new_total.append(scpu+cpu2)
            row={"run":name,"n":256,"repeats":15,
                 "old_search_median":statistics.median(old_search),"new_search_median":statistics.median(new_search),
                 "old_support_median":statistics.median(old_support),"new_support_median":statistics.median(new_support),
                 "old_total_median":statistics.median(old_total),"new_total_median":statistics.median(new_total),
                 "certificate_speedup":statistics.median(old_total)/statistics.median(new_total),
                 "old_total_samples":old_total,"new_total_samples":new_total,
                 "old_search_samples":old_search,"new_search_samples":new_search}
            rows.append(row); (OUT/"summary.json").write_text(json.dumps({"rows":rows,"elapsed_sec":time.time()-start},indent=2))
            print(json.dumps({k:v for k,v in row.items() if "samples" not in k}),flush=True)
    print(f"Completed {len(rows)} six-arm 15-repeat runs: {OUT}")

if __name__=="__main__": main()
