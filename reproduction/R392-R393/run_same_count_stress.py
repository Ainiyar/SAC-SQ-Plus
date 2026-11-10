#!/usr/bin/env python3
"""Independent equal-cardinality path explosion test of R392 exact interval certification.
All expected paths/support are derived analytically from the gadget construction,
not from the certifier. We retain the exact real-input algorithm, unchanged.
"""
from pathlib import Path
import sys, json, re, time, itertools, math, subprocess
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'R392'/'source'))
from r392_full_reproduction import run, path_to_line
H=400_000_000_000_000
DELTA=1_000_000_000_000
CUT=650_000_000_000
B0=100_000_000_000
B1=1_200_000_000_000

def make(m,with_direct=True,with_decoy=True):
    h=2;n=6*m;N=n+1
    L=np.zeros((N,N),dtype=np.int64);U=np.zeros((N,N),dtype=np.int64)
    for i in range(n-h+1):
        L[i,i+h:]=H;U[i,i+h:]=H
    paths=[]
    for bits in itertools.product((0,1),repeat=m):
        path=[]
        for i,choice in enumerate(bits):
            st=6*i;mid=st+2+choice;end=st+6
            path.extend(((st,mid),(mid,end)))
        paths.append(tuple(path))
    for i in range(m):
        st=6*i;end=st+6
        for mid in (st+2,st+3):
            for a,b in ((st,mid),(mid,end)):
                L[a,b]=0;U[a,b]=DELTA
        if with_decoy:
            for a,b in ((st,st+4),(st+4,end)):
                L[a,b]=3*DELTA;U[a,b]=4*DELTA
    if with_direct:
        L[0,n]=U[0,n]=(2*m-1)*CUT
        paths.append(((0,n),))
    return n,h,L,U,paths

def parse_support(log,out):
    lines=Path(log).read_text().splitlines()
    paths=[tuple(tuple(map(int,t.split(':'))) for t in s.split()[1:]) for s in lines if s.startswith('PATH ')]
    supports={}
    for s in Path(out).read_text().splitlines():
        if s.startswith('INDEX '):
            t=s.split();i=int(t[1]);poss=int(t[3]); k=int(t[5]);lo=t[7];hi=t[9]
            if poss: supports[paths[i]]=(k,lo,hi)
    return supports

def main():
    output=[]
    inp=ROOT/'inputs'/'synthetic'; out=ROOT/'results'/'synthetic';inp.mkdir(parents=True,exist_ok=True);out.mkdir(parents=True,exist_ok=True)
    for m in (2,3,4,5,6,7,8):
        name=f'equal_count_m{m}'
        n,h,L,U,truth=make(m)
        (inp/(name+'_L.bin')).write_bytes(L.tobytes());(inp/(name+'_U.bin')).write_bytes(U.tobytes())
        (inp/(name+'_archive_witnesses.txt')).write_text('\n'.join(path_to_line(p) for p in truth)+'\n')
        rec=run(name,n,h,inp,out,ROOT/'bin'/'exact_interval',ROOT/'bin'/'exact_support')
        if rec['status']!='COMPLETE':
            rec['analytical_truth']=len(truth);output.append(rec);break
        supported=parse_support(out/name/'certifier.log',out/name/'exact_support.txt')
        expected=set(truth); assert set(supported)==expected,(m,len(supported),len(expected))
        exp_lo=str(B0);exp_hi=str(CUT);exp_direct_hi=str(B1)
        for P,(k,lo,hi) in supported.items():
            assert k==len(P)
            assert lo==exp_lo, (m,P,lo)
            assert hi==(exp_direct_hi if len(P)==1 else exp_hi),(m,P,hi)
        rec.update(analytical_truth=len(truth),unique_supported_segment_counts=len(set(len(P) for P in supported)),max_multiplicity_same_k=2**m,support_endpoint_verified=True)
        output.append(rec)
    (out/'same_count_summary.json').write_text(json.dumps(output,indent=2))
    print('SUMMARY')
    for z in output:
        print(z['box'],z['independent_seed_count'],z['count_midpoint_candidate_count'],z.get('analytical_truth'),z['status'],z['created_labels'],round(z['total_cpu_all_stages'],4))
if __name__=='__main__':main()
