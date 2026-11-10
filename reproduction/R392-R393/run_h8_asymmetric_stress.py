#!/usr/bin/env python3
"""Asymmetric, h=8, multi-path-per-cardinality strict interval-box benchmark.
Analytical ground truth includes all 2^m alternatives, each with distinct rational
support upper endpoint determined by its own lower-cost scenario, plus direct path.
"""
from pathlib import Path
import sys, json, random, itertools
from fractions import Fraction
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'R392'/'source'))
from r392_full_reproduction import run,path_to_line
from run_same_count_stress import parse_support
H=400_000_000_000_000; DELTA=1_000_000_000_000; CUT=650_000_000_000
B0=100_000_000_000;B1=1_200_000_000_000

def build(m):
  h=8;n=24*m;N=n+1;rng=random.Random(20261110+17*m)
  L=np.zeros((N,N),dtype=np.int64);U=np.zeros((N,N),dtype=np.int64)
  for i in range(n-h+1):L[i,i+h:]=H;U[i,i+h:]=H
  for i in range(m):
    st=24*i;end=st+24
    for mid in (st+8,st+12):
      for a,b in ((st,mid),(mid,end)):
        low=rng.randint(1,200_000_000_000)
        L[a,b]=low;U[a,b]=low+DELTA
    for a,b in ((st,st+16),(st+16,end)):
      L[a,b]=3*DELTA;U[a,b]=4*DELTA
  D=(2*m-1)*CUT;L[0,n]=D;U[0,n]=D
  truths=[]
  for bits in itertools.product((0,1),repeat=m):
    P=[]
    for i,choice in enumerate(bits):
      st=24*i; mid=st+(8 if choice==0 else 12);end=st+24
      P.extend(((st,mid),(mid,end)))
    truths.append(tuple(P))
  truths.append(((0,n),))
  return n,h,L,U,truths

def main():
  inp=ROOT/'inputs'/'h8_asym';out=ROOT/'results'/'h8_asym'
  inp.mkdir(parents=True,exist_ok=True);out.mkdir(parents=True,exist_ok=True)
  all=[]
  for m in (3,4,5,6,7,8):
    name=f'h8_asymmetric_m{m}';n,h,L,U,truth=build(m)
    (inp/(name+'_L.bin')).write_bytes(L.tobytes());(inp/(name+'_U.bin')).write_bytes(U.tobytes())
    (inp/(name+'_archive_witnesses.txt')).write_text('\n'.join(path_to_line(P) for P in truth)+'\n')
    rec=run(name,n,h,inp,out,ROOT/'bin'/'exact_interval',ROOT/'bin'/'exact_support')
    if rec['status']=='COMPLETE':
      supp=parse_support(out/name/'certifier.log',out/name/'exact_support.txt')
      assert set(supp)==set(truth),(m,len(supp),len(truth))
      upper_set=set()
      for P,(k,lo,hi) in supp.items():
        assert k==len(P)
        assert Fraction(lo)==Fraction(B0)
        if k==1:expected=Fraction(B1)
        else:
          cost=sum(int(L[a,b]) for a,b in P)
          expected=Fraction((2*m-1)*CUT-cost,2*m-1)
          upper_set.add(expected)
        assert Fraction(hi)==expected,(m,P,hi,expected)
      rec.update(analytical_truth=len(truth),support_endpoint_verified=True,distinct_nontrivial_upper_endpoints=len(upper_set),max_multiplicity_same_k=2**m)
    all.append(rec)
    print('h8',m,rec['status'],rec['independent_seed_count'],rec['count_midpoint_candidate_count'],rec['created_labels'],rec.get('supported_paths'),rec.get('distinct_nontrivial_upper_endpoints'),flush=True)
  (out/'h8_asymmetric_summary.json').write_text(json.dumps(all,indent=2))
if __name__=='__main__':main()
