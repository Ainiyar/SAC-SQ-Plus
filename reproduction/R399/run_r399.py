#!/usr/bin/env python3
"""R399 same-box/same-witness old16 vs exact continuous comparison.
Archived output is never passed to witness initialization or C++ certifiers.
"""
from pathlib import Path
import importlib.util, json, os, re, resource, statistics, subprocess, sys, time
import numpy as np
R=Path(__file__).resolve().parent; IN=R.parent/'R392-R393'/'real_boxes'/'inputs'; REF=R/'r392/results'; OUT=R/'runs';OUT.mkdir(exist_ok=True)
B0=100_000_000_000;B1=1_200_000_000_000
spec=importlib.util.spec_from_file_location('r392',R/'r392/source/r392_full_reproduction.py'); gen=importlib.util.module_from_spec(spec);spec.loader.exec_module(gen)
def cput():
 x=resource.getrusage(resource.RUSAGE_CHILDREN);return x.ru_utime+x.ru_stime
def invoke(cmd):
 t=cput();w=time.perf_counter();z=subprocess.run(cmd,capture_output=True,text=True,timeout=170);return z,cput()-t,time.perf_counter()-w
def parse(s):
 st=re.search(r'^status (\S+)',s,re.M); cr=re.search(r'^pops \d+ extensions \d+ created (\d+)',s,re.M)
 return (st.group(1) if st else 'ERROR',int(cr.group(1)) if cr else -1,[x[5:].strip() for x in s.splitlines() if x.startswith('PATH ')])
def exact_paths(log,support):
 p=parse(log)[2];ans={}
 for line in support.splitlines():
  if not line.startswith('INDEX '):continue
  z=line.split()
  if z[3]=='1':ans[p[int(z[1])]]=(int(z[5]),int(z[7]),int(z[9]))
 return ans
def init(box,n,rounds=5):
 L=np.fromfile(IN/(box+'_L.bin'),np.int64).reshape(n+1,n+1);U=np.fromfile(IN/(box+'_U.bin'),np.int64).reshape(n+1,n+1);M=L+(U-L)//2
 trials=[];first=None;numseeds=0
 for rep in range(rounds):
  t=time.process_time();seeds=[]
  for C in (L,M,U):
   for i in range(17):
    p=gen.fixed_beta_best(C,n,8,B0+((B1-B0)*i)//16)
    if p not in seeds:seeds.append(p)
  ts=time.process_time()-t; numseeds=len(seeds)
  lo=min(map(len,seeds));hi=max(map(len,seeds));t=time.process_time(); nxt=gen.midpoint_exact_count(M,n,8,hi);proposals=list(seeds)
  for k in range(lo,hi+1):
   p=gen.recover_k(nxt,n,k)
   if p and p not in proposals:proposals.append(p)
  tc=time.process_time()-t;lines=[gen.path_to_line(p) for p in proposals]
  if first is None:first=lines
  assert lines==first,(box,'nondeterminism')
  trials.append({'seed_cpu':ts,'count_cpu':tc,'total_cpu':ts+tc})
 old=[s.strip() for s in (REF/box/'count_midpoint_witnesses.txt').read_text().splitlines() if s.strip()]
 assert first==old,(box,'R392 seed regeneration mismatch')
 (OUT/box).mkdir(exist_ok=True)
 (OUT/box/'independent_witnesses.txt').write_text('\n'.join(first)+'\n')
 return {'name':box,'n':n,'seed_count':numseeds,'pool_count':len(first),'trial_cpu':trials,'initialization_cpu_median':statistics.median(t['total_cpu'] for t in trials),'regenerated_matches_R392':True}
def checkpoint(rows,inits):
 (OUT/'R399_RAW_RUNS.json').write_text(json.dumps({'runs':rows,'inits':inits},indent=2))
def main():
 start=time.time();rows=[];inits=[];names=sorted(x.name[:-6] for x in IN.glob('*_L.bin'))
 for j,box in enumerate(names,1):
  n=1024 if ('1024' in box or 'R343' in box) else 256
  repeats=5 if n==1024 else 15
  z=init(box,n);inits.append(z)
  reference=exact_paths((REF/box/'certifier.log').read_text(),(REF/box/'exact_support.txt').read_text())
  L=IN/(box+'_L.bin');U=IN/(box+'_U.bin');W=OUT/box/'independent_witnesses.txt'
  args=[str(n),'8',str(B0),str(B1),str(L),str(U),'500000','120',str(W)]
  for rep in range(repeats):
   methods=('baseline16','exact_interval') if rep%2==0 else ('exact_interval','baseline16')
   for method in methods:
    cmd=[str(R/method),*args,*(['16'] if method=='baseline16' else [])]
    p,cpu,wall=invoke(cmd);status,created,paths=parse(p.stdout)
    f=OUT/box/f'{method}_rep{rep:02}.log';f.write_text(p.stdout+'\n'+p.stderr)
    row={'box':box,'n':n,'rep':rep,'method':method,'returncode':p.returncode,'status':status,'created_labels':created,'candidate_paths':len(paths),'cert_cpu':cpu,'cert_wall':wall}
    if status=='COMPLETE' and p.returncode==0:
     sf=OUT/box/f'{method}_rep{rep:02}.support.txt'
     q,scpu,swall=invoke([str(R/'exact_support'),str(n),'8',str(B0),str(B1),str(L),str(U),str(f),str(sf)])
     if q.returncode:raise RuntimeError(('support_failed',box,method,rep,q.stderr))
     actual=exact_paths(p.stdout,sf.read_text());eq=actual==reference
     row.update({'support_cpu':scpu,'support_wall':swall,'cert_plus_support_cpu':cpu+scpu,'supported_paths':len(actual),'identity_support_equal_R392':eq})
    rows.append(row)
    if status!='COMPLETE' or p.returncode or not row.get('identity_support_equal_R392',False):
     checkpoint(rows,inits);raise RuntimeError(('DIFFERENCE',row))
   checkpoint(rows,inits)
  a=[r for r in rows if r['box']==box and r['method']=='baseline16'];b=[r for r in rows if r['box']==box and r['method']=='exact_interval']
  print(f'BOX_DONE {j}/11 {box} old_labels={a[0]["created_labels"]} new_labels={b[0]["created_labels"]} old_med={statistics.median(r["cert_plus_support_cpu"] for r in a):.6f} new_med={statistics.median(r["cert_plus_support_cpu"] for r in b):.6f} elapsed={time.time()-start:.1f}s',flush=True)
 groups=[]
 for z in inits:
  box=z['name'];a=[r for r in rows if r['box']==box and r['method']=='baseline16'];b=[r for r in rows if r['box']==box and r['method']=='exact_interval']
  med=lambda rr,key:statistics.median(x[key] for x in rr)
  ia=z['initialization_cpu_median'];o=med(a,'cert_plus_support_cpu');e=med(b,'cert_plus_support_cpu')
  groups.append({'box':box,'n':z['n'],'repetitions':len(a),'seed_count':z['seed_count'],'witness_count':z['pool_count'],'initialization_cpu_median':ia,
   'old_labels':a[0]['created_labels'],'new_labels':b[0]['created_labels'],'old_candidates':a[0]['candidate_paths'],'new_candidates':b[0]['candidate_paths'],
   'final_possible_paths':a[0]['supported_paths'],'all_complete_and_exact_support_equal':all(r['identity_support_equal_R392'] for r in a+b),
   'old_search_cpu_median':med(a,'cert_cpu'),'new_search_cpu_median':med(b,'cert_cpu'),'old_support_cpu_median':med(a,'support_cpu'),'new_support_cpu_median':med(b,'support_cpu'),
   'old_cert_support_cpu_median':o,'new_cert_support_cpu_median':e,'old_full_cpu':o+ia,'new_full_cpu':e+ia,'cert_speedup':o/e,'full_speedup':(o+ia)/(e+ia)})
 aggregate={'cases':len(groups),'all_complete':all(x['all_complete_and_exact_support_equal'] for x in groups),'geo_cert_speedup':statistics.geometric_mean(g['cert_speedup'] for g in groups),'geo_full_speedup':statistics.geometric_mean(g['full_speedup'] for g in groups),'geo_label_reduction':statistics.geometric_mean(g['old_labels']/g['new_labels'] for g in groups)}
 result={'protocol':{'same_L_U':True,'archive_independent_seed_initialization':True,'beta_bounds':[B0,B1],'h':8,'baseline_cells':16,'max_labels':500000,'time_limit':120,'new_old_return_contract':'COMPLETE then R354 exact supports','rep_n256':15,'rep_n1024':5},'system':{'platform':os.uname().sysname,'cpu_arch':os.uname().machine},'cases':groups,'aggregate':aggregate,'elapsed_seconds':time.time()-start}
 (OUT/'R399_SUMMARY.json').write_text(json.dumps(result,indent=2));print('FINISHED '+json.dumps(aggregate)+' elapsed='+str(result['elapsed_seconds']),flush=True)
if __name__=='__main__': main()
