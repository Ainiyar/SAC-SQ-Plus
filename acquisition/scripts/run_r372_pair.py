#!/usr/bin/env python3
"""Run the R372 baseline and certificate-bank-pricing arms in separate processes."""
from pathlib import Path
import json, os, resource, subprocess, sys, time
ACQ=Path(__file__).resolve().parents[1]
SCRIPT=Path(__file__).resolve().with_name("run_r372_arm.py")
OUT=ACQ/"data"/"fresh_runs"/"r372_pair"

def child_cpu():
 r=resource.getrusage(resource.RUSAGE_CHILDREN);return r.ru_utime+r.ru_stime

def main():
 OUT.mkdir(parents=True,exist_ok=True); rows=[]
 env=os.environ.copy()
 for k in ["OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"]:env[k]="1"
 for arm in ["baseline","certificate_bank_pricing"]:
  result=OUT/f"{arm}.json";log=OUT/f"{arm}.log";before=child_cpu();t=time.perf_counter()
  cmd=[sys.executable,str(SCRIPT),"--arm",arm,"--out",str(result)]
  if arm=="certificate_bank_pricing":cmd += ["--cost-model-from",str((OUT/"baseline.json").resolve())]
  with log.open("w",encoding="utf-8") as f:p=subprocess.run(cmd,cwd=ACQ,env=env,stdout=f,stderr=subprocess.STDOUT)
  cpu=child_cpu()-before;wall=time.perf_counter()-t
  if p.returncode:raise RuntimeError(f"R372 {arm} failed; inspect {log}")
  data=json.loads(result.read_text())
  row={"arm":arm,"status":data.get("status"),"total_pipeline_cpu_s":data.get("total_cpu"),
       "actions":data.get("actions"),"archive_paths":len(data.get("archive_paths",[])),
       "process_cpu_seconds":cpu,"wall_seconds":wall,"result":result.name,"log":log.name}
  rows.append(row);(OUT/"summary.json").write_text(json.dumps({"runs":rows},indent=2));print(json.dumps(row),flush=True)
 print(f"R372 pair completed: {OUT}")

if __name__=="__main__":main()
