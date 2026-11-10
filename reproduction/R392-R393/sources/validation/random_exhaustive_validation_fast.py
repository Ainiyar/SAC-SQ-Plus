from pathlib import Path
from fractions import Fraction
import numpy as np, random, subprocess, tempfile, json, time, math, os
RELEASE_ROOT = Path(__file__).resolve().parents[2]
BIN = os.environ.get("SAC_CERT_BIN", str(RELEASE_ROOT / "bin" / "exact_interval"))
OUT = Path(os.environ.get("SAC_VALIDATION_OUT", str(RELEASE_ROOT / "output" / "random_validation")))
OUT.mkdir(parents=True, exist_ok=True)

def all_paths(n,h):
    out=[]
    def rec(s,edges):
        if s==n: out.append(tuple(edges)); return
        for t in range(s+h,n+1): rec(t,edges+[(s,t)])
    rec(0,[])
    return out

def possible_paths(paths,L,U,b0,b1,n):
    # Exact best-scenario support for every path, independently by exhaustive path enumeration.
    possible=set(); supports={}
    for P in paths:
        Pset=set(P); kP=len(P); A=sum(int(L[s,t]) for s,t in P)
        m_by_k={}
        for Q in paths:
            k=len(Q); val=sum(int(L[s,t]) if (s,t) in Pset else int(U[s,t]) for s,t in Q)
            m_by_k[k]=min(m_by_k.get(k,10**100),val)
        if m_by_k.get(kP)!=A: continue
        lo,hi=Fraction(b0),Fraction(b1)
        for k,m in m_by_k.items():
            if k>kP: lo=max(lo,Fraction(A-m,k-kP))
            elif k<kP: hi=min(hi,Fraction(m-A,kP-k))
        if lo<=hi:
            possible.add(P); supports[P]=(lo,hi)
    return possible,supports

def pathstr(P): return ' '.join(f'{s}:{t}' for s,t in P)

def run_case(idx,negative):
    rng=random.Random(0x5AC000+idx+(100000 if negative else 0))
    n=rng.randint(10,16); h=rng.randint(3,5)
    paths=all_paths(n,h)
    if len(paths)>15000: raise RuntimeError('too many')
    N=n+1; L=np.zeros((N,N),dtype='<i8'); U=np.zeros((N,N),dtype='<i8')
    for s in range(n-h+1):
        for t in range(s+h,n+1):
            l=rng.randint(0,10000); u=l+rng.randint(0,5000)
            L[s,t]=l; U[s,t]=u
    if negative:
        b0=rng.randint(-30000,-1000); b1=rng.randint(-900,-1)
        if b0>=b1: b0,b1=-10000,-100
    else:
        b0=rng.randint(0,1000); b1=rng.randint(max(1001,b0+1),30000)
    Wn=min(len(paths),rng.randint(1, min(5,len(paths))))
    witnesses=rng.sample(paths,Wn)
    case=OUT/(f'case_{idx:04d}_'+('neg' if negative else 'pos'))
    L.tofile(str(case)+'.L.bin'); U.tofile(str(case)+'.U.bin')
    Path(str(case)+'.w.txt').write_text('\n'.join(pathstr(p) for p in witnesses)+'\n')
    p=subprocess.run([BIN,str(n),str(h),str(b0),str(b1),str(case)+'.L.bin',str(case)+'.U.bin','1000000','10',str(case)+'.w.txt'],capture_output=True,text=True)
    if p.returncode: return {'idx':idx,'negative':negative,'error':p.stderr,'stdout':p.stdout[-1000:]}
    status=next((x.split()[1] for x in p.stdout.splitlines() if x.startswith('status ')),'?')
    cand=[tuple((int(x.split(':')[0]),int(x.split(':')[1])) for x in line.split()[1:]) for line in p.stdout.splitlines() if line.startswith('PATH ')]
    candset=set(cand)
    true,supports=possible_paths(paths,L,U,b0,b1,n)
    filtered={P for P in candset if P in true}
    missing=true-candset
    extra=candset-true
    result={'idx':idx,'negative':negative,'n':n,'h':h,'all_paths':len(paths),'witnesses':Wn,'status':status,'candidate_paths':len(candset),'true_possible':len(true),'missing_true':len(missing),'extra_candidates':len(extra),'filtered_count':len(filtered),'ok':status=='COMPLETE' and not missing,
            'b0':b0,'b1':b1}
    if not result['ok']:
        result['missing_examples']=[pathstr(x) for x in list(missing)[:3]]
        result['extra_examples']=[pathstr(x) for x in list(extra)[:3]]
        result['solver_summary']='\n'.join(p.stdout.splitlines()[:8])
    return result

if __name__=='__main__':
    t=time.time(); results=[]
    # 500 positive-beta and 500 negative-beta exhaustive cases.
    for i in range(500):
        r=run_case(i,False); results.append(r)
        if not r.get('ok',False): print('FAIL',json.dumps(r),flush=True); break
        if (i+1)%100==0: print('positive',i+1,'elapsed',round(time.time()-t,2),flush=True)
    pos_ok=all(x.get('ok',False) for x in results)
    if pos_ok:
        for i in range(500):
            r=run_case(i,True); results.append(r)
            if not r.get('ok',False): print('FAIL',json.dumps(r),flush=True); break
            if (i+1)%100==0: print('negative',i+1,'elapsed',round(time.time()-t,2),flush=True)
    summary={'positive_cases':sum(1 for x in results if not x['negative']),'negative_cases':sum(1 for x in results if x['negative']),
             'all_ok':all(x.get('ok',False) for x in results),'failures':[x for x in results if not x.get('ok',False)],'elapsed_sec':time.time()-t}
    (OUT/'summary_1000.json').write_text(json.dumps({'summary':summary,'results':results},indent=2))
    print(json.dumps(summary),flush=True)
