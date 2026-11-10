from pathlib import Path
from fractions import Fraction as F
import itertools, json, struct
W=Path(__file__).resolve().parent
n=24; N=25; h=1; m=8; Q=10**12; eps=1; H=1000*Q; beta0=0; beta1=12*Q//10
L=list(struct.unpack('<'+'q'*(N*N),(W/'L.bin').read_bytes()))
U=list(struct.unpack('<'+'q'*(N*N),(W/'U.bin').read_bytes()))

def path_for(bits):
    out=[]
    for i,bit in enumerate(bits):
        b=3*i
        out += ([(b,b+1),(b+1,b+3)] if bit==0 else [(b,b+2),(b+2,b+3)])
    return tuple(out)
expected={path_for(bits) for bits in itertools.product([0,1],repeat=m)}
# Parse R357 paths.
r357=(W/'R357.log').read_text().splitlines()
status=next(x.split()[1] for x in r357 if x.startswith('status '))
created=int(next(x.split()[5] for x in r357 if x.startswith('pops ')))
pruned=int(next(x.split()[3] for x in r357 if x.startswith('archive_witnesses ')))
actual=set()
for line in r357:
    if line.startswith('PATH '):
        actual.add(tuple(tuple(map(int,t.split(':'))) for t in line.split()[1:]))
assert status=='COMPLETE' and actual==expected and len(actual)==256
# Parse exact R354 results.
rows=[x for x in (W/'R354_R357.txt').read_text().splitlines() if x.startswith('INDEX ')]
assert len(rows)==256
assert all('POSSIBLE 1' in x and 'LO 0 HI 1200000000000' in x and 'SAME_COUNT 1' in x for x in rows)
# Independently evaluate fixed interpolation at selected t and beta values by exact DAG DP.
all_a=path_for([0]*m)
checks=[]
for t in [F(0),F(1,2),F(1)]:
    for beta in [F(0),F(3,5),F(6,5)]:
        d=[None]*(n+1); cnt=[0]*(n+1); pred=[None]*(n+1)
        d[0]=F(0);cnt[0]=1
        for v in range(1,n+1):
            best=None; ways=0; best_s=None
            for s in range(v):
                if d[s] is None: continue
                idx=s*N+v
                c=F(L[idx])+t*F(U[idx]-L[idx])+beta*Q
                val=d[s]+c
                if best is None or val<best:
                    best=val;ways=cnt[s];best_s=s
                elif val==best:
                    ways+=cnt[s]
            d[v]=best;cnt[v]=ways;pred[v]=best_s
        path=[];v=n
        while v:
            s=pred[v];path.append((s,v));v=s
        path=tuple(reversed(path))
        assert path==all_a and cnt[n]==1
        checks.append({'t':str(t),'beta':str(beta),'unique_all_A':True,'scaled_objective':str(d[n])})
report={
 'n':n,'h':h,'m':m,'epsilon':eps,'H_scaled':H,'beta_interval':[0,'1.2'],
 'high_edge_condition_pass':H>2*m*(beta1+eps*Q),
 'R357_status':status,'R357_created_labels':created,'R357_witness_pruned_prefixes':pruned,
 'R357_sink_paths':len(actual),'R357_equals_all_2^m_block_paths':True,
 'R354_candidates':len(rows),'R354_possible':sum('POSSIBLE 1' in x for x in rows),
 'R354_all_supports_exactly_[0,1.2]':True,'uniform_interpolation_checks':checks,
 'theorem_uniform_interpolation_envelope':'unique all-A path for all t in [0,1] and beta in [0,1.2]'
}
(W/'R385_M8_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
(W/'R385_M8_VERIFICATION.md').write_text('''# R385 m=8 uniform-interpolation separation: executable validation\n\n## Construction\n\nEight blocks of length `3h` (`h=1`) each contain two two-edge routes. Route A edges have intervals `[0,1]`; route B edges have `[0.5,1.5]`; all other legal edges have fixed cost `H=1000`. The shared penalty is `beta in [0,1.2]`, and costs/penalties are represented at scale `Q=10^12`. The sufficient high-edge condition `H > 2m(beta_1+epsilon)` holds.\n\n## R357/R354 result\n\nThe actual R357 executable returned `COMPLETE`, created 1,021 labels, pruned 212 prefixes using three feasible witnesses and returned exactly the 256 blockwise A/B paths. The R354 exact support executable processed all 256 candidates and marked all 256 possible, each with exact closed support `[0,1.2]`. The returned R357 path set was compared with an independently generated set of all `2^8` route choices.\n\n## Uniform interpolation check\n\nAn independent exact-rational DAG dynamic program was run at `t in {0, 1/2, 1}` and `beta in {0, 3/5, 6/5}`. All nine checks returned the unique all-A path. The proposition in the manuscript establishes the same result for every real `t in [0,1]` and every `beta in [0,1.2]`; the finite checks are implementation cross-checks, not a replacement for the proof.\n\nSee `R357.log`, `R354_R357.txt`, `R354_all256.txt`, `L.bin`, `U.bin`, `witnesses.txt`, and `all_256_candidates.log` for the complete raw output and inputs. `verify_results.py` reruns the output comparison and fixed-interpolation checks; `build_and_run.py` rebuilds R354 from the archived source and reruns both algorithms (its binary paths are local to the original workspace).\n''')
print(json.dumps(report,indent=2))
