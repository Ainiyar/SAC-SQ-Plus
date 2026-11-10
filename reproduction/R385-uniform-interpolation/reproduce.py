#!/usr/bin/env python3
"""Portable R357/R354 reproduction for the m=8 uniform-interpolation gadget."""
from pathlib import Path
from fractions import Fraction as F
import itertools, struct, subprocess
W=Path(__file__).resolve().parent
n=24; N=25; h=1; m=8; Q=10**12; eps=1; beta0=0; beta1=12*Q//10; H=1000*Q
L=[0]*(N*N); U=[0]*(N*N)
for s in range(n+1):
    for t in range(s+h,n+1):
        L[s*N+t]=H; U[s*N+t]=H

def path_for(bits):
    out=[]
    for i,bit in enumerate(bits):
        b=3*i
        out += ([(b,b+1),(b+1,b+3)] if bit==0 else [(b,b+2),(b+2,b+3)])
    return out
A=[];B=[]
for i in range(m):
    b=3*i; ae=[(b,b+1),(b+1,b+3)]; be=[(b,b+2),(b+2,b+3)]
    for s,t in ae: L[s*N+t]=0; U[s*N+t]=eps*Q
    for s,t in be: L[s*N+t]=eps*Q//2; U[s*N+t]=3*eps*Q//2
    A+=ae; B+=be
for name,arr in [('L',L),('U',U)]:
    (W/f'{name}.bin').write_bytes(struct.pack('<'+'q'*len(arr),*arr))
wits=[path_for([0]*m),path_for([1]*m),path_for([i%2 for i in range(m)])]
(W/'witnesses.txt').write_text('\n'.join(' '.join(f'{s}:{t}' for s,t in p) for p in wits)+'\n')
allpaths=[path_for(bits) for bits in itertools.product([0,1],repeat=m)]
(W/'all_256_candidates.log').write_text('\n'.join('PATH '+' '.join(f'{s}:{t}' for s,t in p) for p in allpaths)+'\n')
R357=W/'bin/verifier_r357_reference'; R354=W/'bin/exact_support_reference'
R357.parent.mkdir(exist_ok=True)
subprocess.run(['g++','-O3','-std=c++17',str(W/'src/aug_tropical_multiwitness.cpp'),'-o',str(R357)],check=True)
subprocess.run(['g++','-O3','-std=c++17',str(W/'src/exact_best_scenario_support.cpp'),'-o',str(R354)],check=True)
with (W/'R357.log').open('w') as f:
    subprocess.run([str(R357),str(n),str(h),str(beta0),str(beta1),str(W/'L.bin'),str(W/'U.bin'),'2000000','10',str(W/'witnesses.txt'),'16'],stdout=f,stderr=subprocess.STDOUT,check=True)
for cand,out in [('R357.log','R354_R357.txt'),('all_256_candidates.log','R354_all256.txt')]:
    with (W/(out.replace('.txt','.log'))).open('w') as f:
        subprocess.run([str(R354),str(n),str(h),str(beta0),str(beta1),str(W/'L.bin'),str(W/'U.bin'),str(W/cand),str(W/out)],stdout=f,stderr=subprocess.STDOUT,check=True)
subprocess.run(['python3',str(W/'verify_results.py')],check=True)
print((W/'R385_M8_VERIFICATION.md').read_text())
