#!/usr/bin/env python3
"""SAC R392: archive-independent cardinality-conditioned witness discovery + R390 certifier.
Reads final archive only after candidate generation and certification, exclusively for reporting.
The initial beta-node sweep is a proposal heuristic. COMPLETE is established by R390's
exhausted search and exact R354 best-scenario support, not by the seed-generation heuristic.
"""
import argparse, json, re, resource, subprocess, time
from pathlib import Path
import numpy as np
INF=np.int64(2**62)
B0=100_000_000_000; B1=1_200_000_000_000

def path_to_line(p):return ' '.join(f'{s}:{t}' for s,t in p)
def fixed_beta_best(C,n,h,beta):
    cost=np.full(n+1,INF,dtype=np.int64);nxt=np.full(n+1,-1,dtype=np.int32);cost[n]=0
    for v in range(n-h,-1,-1):
        ts=np.arange(v+h,n+1)
        vals=C[v,ts]+beta+cost[ts]
        j=int(np.argmin(vals));cost[v]=vals[j];nxt[v]=ts[j]
    P=[];v=0
    while v<n:
        t=int(nxt[v]);P.append((v,t));v=t
    return tuple(P)

def midpoint_exact_count(C,n,h,kmax):
    dp=np.full((kmax+1,n+1),INF,dtype=np.int64)
    nxt=np.full((kmax+1,n+1),-1,dtype=np.int32);dp[0,n]=0
    for k in range(1,kmax+1):
        for v in range(n-h,-1,-1):
            ts=np.arange(v+h,n+1);mask=dp[k-1,ts]<INF
            if not mask.any():continue
            ts=ts[mask];vals=C[v,ts]+dp[k-1,ts]
            j=int(np.argmin(vals));dp[k,v]=vals[j];nxt[k,v]=ts[j]
    return nxt

def recover_k(nxt,n,k):
    out=[];v=0
    while k>0:
        t=int(nxt[k,v])
        if t<0:return None
        out.append((v,t));v=t;k-=1
    return tuple(out) if v==n else None

def child_cpu():
    z=resource.getrusage(resource.RUSAGE_CHILDREN)
    return z.ru_utime+z.ru_stime

def run(box,n,h,inputdir,outdir,cert_bin,support_bin):
    Lp=inputdir/(box+'_L.bin');Up=inputdir/(box+'_U.bin')
    L=np.fromfile(Lp,dtype=np.int64).reshape(n+1,n+1)
    U=np.fromfile(Up,dtype=np.int64).reshape(n+1,n+1)
    M=L+(U-L)//2
    t=time.process_time()
    seeds=[]
    for C in (L,M,U):
        for i in range(17):
            p=fixed_beta_best(C,n,h,B0+((B1-B0)*i)//16)
            if p not in seeds:seeds.append(p)
    t_seed=time.process_time()-t
    kmin=min(map(len,seeds));kmax=max(map(len,seeds))
    t=time.process_time()
    nxt=midpoint_exact_count(M,n,h,kmax)
    proposal=list(seeds)
    for k in range(kmin,kmax+1):
        p=recover_k(nxt,n,k)
        if p and p not in proposal:proposal.append(p)
    t_count=time.process_time()-t
    out=outdir/box;out.mkdir(parents=True,exist_ok=True)
    (out/'independent_LMU_seeds.txt').write_text('\n'.join(map(path_to_line,seeds))+'\n')
    (out/'count_midpoint_witnesses.txt').write_text('\n'.join(map(path_to_line,proposal))+'\n')
    base=[str(n),str(h),str(B0),str(B1),str(Lp),str(Up)]
    t=child_cpu();p=subprocess.run([str(cert_bin),*base,'500000','40',str(out/'count_midpoint_witnesses.txt')],capture_output=True,text=True)
    t_cert=child_cpu()-t;(out/'certifier.log').write_text(p.stdout+'\n'+p.stderr)
    status=re.search(r'^status (\S+)',p.stdout,re.M)
    created=re.search(r'^pops \d+ extensions \d+ created (\d+)',p.stdout,re.M)
    paths=[x[5:].strip() for x in p.stdout.splitlines() if x.startswith('PATH ')]
    rec={'box':box,'n':n,'kmin':kmin,'kmax':kmax,'independent_seed_count':len(seeds),'count_midpoint_candidate_count':len(proposal),'new_count_candidates':len(proposal)-len(seeds),'seed_cpu':t_seed,'midpoint_count_dp_cpu':t_count,'cert_cpu':t_cert,'status':status.group(1) if status else 'ERROR','created_labels':int(created.group(1)) if created else -1,'sink_candidates':len(paths),'cert_returncode':p.returncode}
    if rec['status']=='COMPLETE':
        t=child_cpu();q=subprocess.run([str(support_bin),*base,str(out/'certifier.log'),str(out/'exact_support.txt')],capture_output=True,text=True)
        rec['support_cpu']=child_cpu()-t;rec['support_returncode']=q.returncode
        supported={paths[int(z.split()[1])] for z in (out/'exact_support.txt').read_text().splitlines() if z.startswith('INDEX ') and z.split()[3]=='1'}
        rec['supported_paths']=len(supported)
    # Archive comparison only after entire reconstruction has finished.
    archived={z.strip() for z in (inputdir/(box+'_archive_witnesses.txt')).read_text().splitlines() if z.strip()}
    rec['archive_count']=len(archived)
    rec['candidate_archive_intersection']=len(set(map(path_to_line,proposal))&archived)
    rec['certified_identity_equals_archive']=rec['status']=='COMPLETE' and supported==archived
    rec['total_cpu_all_stages']=sum(rec[k] for k in ('seed_cpu','midpoint_count_dp_cpu','cert_cpu','support_cpu') if k in rec)
    (out/'result.json').write_text(json.dumps(rec,indent=2))
    print(json.dumps(rec),flush=True)
    return rec

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--cert',type=Path,required=True);p.add_argument('--support',type=Path,required=True)
    p.add_argument('--boxes',nargs='*',default=[])
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    boxes=sorted(x.name[:-6] for x in a.inputs.glob('*_L.bin') if (a.inputs/(x.name[:-6]+'_archive_witnesses.txt')).exists())
    if a.boxes:boxes=[x for x in boxes if any(t in x for t in a.boxes)]
    results=[]
    for box in boxes:
        n=1024 if 'n1024' in box or 'R343' in box else 256
        results.append(run(box,n,8,a.inputs,a.out,a.cert,a.support))
    (a.out/'R392_SUMMARY.json').write_text(json.dumps(results,indent=2))
if __name__=='__main__':main()
