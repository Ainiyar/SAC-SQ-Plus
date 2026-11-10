"""An exact singleton-support path at an unsampled beta, with three segment counts.
Count-DP recovers k=2 singleton that 17-point penalty sweep does not discover.
"""
from pathlib import Path
import sys,json
from fractions import Fraction
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'R392'/'source'))
from r392_full_reproduction import run,path_to_line
from run_same_count_stress import parse_support
n=24;h=8;N=n+1;H=400_000_000_000_000
TAU=665_000_000_000; AP=1_000_000_000_000
L=np.zeros((N,N),dtype=np.int64);U=np.zeros((N,N),dtype=np.int64)
for s in range(n-h+1):L[s,s+h:]=H;U[s,s+h:]=H
Q1=((0,8),(8,16),(16,24))
P=((0,12),(12,24))
Q2=((0,24),)
for e,c in [(Q1[0],AP-TAU),(Q1[1],0),(Q1[2],0),(P[0],AP),(P[1],0),(Q2[0],AP+TAU)]:
   a,b=e;L[a,b]=U[a,b]=c
inp=ROOT/'inputs'/'singleton';out=ROOT/'results'/'singleton';inp.mkdir(parents=True,exist_ok=True);out.mkdir(parents=True,exist_ok=True)
name='threeway_singleton_tau665'
(inp/(name+'_L.bin')).write_bytes(L.tobytes());(inp/(name+'_U.bin')).write_bytes(U.tobytes())
(inp/(name+'_archive_witnesses.txt')).write_text('\n'.join(path_to_line(x) for x in (Q1,P,Q2))+'\n')
rec=run(name,n,h,inp,out,ROOT/'bin'/'exact_interval',ROOT/'bin'/'exact_support')
supp=parse_support(out/name/'certifier.log',out/name/'exact_support.txt')
assert rec['independent_seed_count']==2 and rec['count_midpoint_candidate_count']==3,rec
assert supp[P]==(2,str(TAU),str(TAU)),supp
assert supp[Q1]==(3,str(100_000_000_000),str(TAU)),supp
assert supp[Q2]==(1,str(TAU),str(1_200_000_000_000)),supp
assert rec['status']=='COMPLETE'
rec['analytical_supports_verified']=True
rec['singleton_path_not_in_initial_seeds']=True
(out/'singleton_summary.json').write_text(json.dumps(rec,indent=2))
print('SINGLETON SUCCESS',json.dumps(rec,indent=2))
