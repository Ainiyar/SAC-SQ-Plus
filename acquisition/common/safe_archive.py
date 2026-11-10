"""Safe bilateral monotone all-penalty archive, research implementation R333.

Permanent A stores every discovered path. Complement proofs exclude A, while
every A-current path needs a strict all-beta retirement certificate. Never
delete a path from A: that would change the domain of cached complement proofs.
"""
import argparse
import json
import sys
import time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from numba import njit

HERE = Path(__file__).resolve().parent
PACKAGE = HERE / 'input' / 'SAC_SOL_HANDOFF_20261107'
if not PACKAGE.exists():
    PACKAGE = HERE.parent / 'SAC_SOL_REVIEW' / 'SAC_SOL_HANDOFF_20261107'
sys.path.insert(0, str(PACKAGE / 'R326_CORE'))
sys.path.insert(0, str(PACKAGE / 'R325_SUPPORT' / 'support'))
from R326_SHARED_DP import all_k_with_upper_pred, path_from_pred
from R326_HULL_ALGORITHM import count_hull
import interval_arith as ia
from scipy.special import expit
Q = 10**12
INF = np.int64(1 << 60)


@njit(cache=True)
def complement(L, U, mask, edges, bn, bd, nv, cm, ps, pt, pl, h):
    n = L.shape[0]-1
    d = np.full(n+1, INF, np.int64)
    nx = np.full(n+1, -1, np.int32)
    d[n] = 0
    for s in range(n-h, -1, -1):
        for t in range(s+h, n+1):
            if d[t] == INF:
                continue
            w = U[s,t] if mask[s,t] else L[s,t]
            val = w*bd + bn + d[t]
            if val < d[s]:
                d[s] = val
                nx[s] = t
    best = INF
    node = -1
    target = -1
    for z in range(len(nv)):
        s = nv[z]
        pref = np.int64(0)
        for j in range(pl[z]):
            a,b = ps[z,j],pt[z,j]
            w = U[a,b] if mask[a,b] else L[a,b]
            pref += w*bd+bn
        for t in range(s+h,n+1):
            if cm[z,t] or d[t] == INF:
                continue
            w = U[s,t] if mask[s,t] else L[s,t]
            val = pref+w*bd+bn+d[t]
            if val < best:
                best,node,target = val,z,t
    cand = np.int64(0)
    for a,b in edges:
        cand += U[a,b]*bd+bn
    if node < 0:
        return INF,-1,-1,nx
    return best-cand,node,target,nx


def trie(paths,n):
    nodes=[{'v':0,'children':{},'prefix':()}]
    for path in paths:
        idx=0
        for edge in path:
            if edge[1] not in nodes[idx]['children']:
                nodes[idx]['children'][edge[1]]=len(nodes)
                nodes.append({'v':edge[1],'children':{},
                              'prefix':nodes[idx]['prefix']+(edge,)})
            idx=nodes[idx]['children'][edge[1]]
    depth=max(len(x['prefix']) for x in nodes)
    nv=np.array([x['v'] for x in nodes],np.int32)
    cm=np.zeros((len(nodes),n+1),np.uint8)
    ps=np.zeros((len(nodes),depth),np.int32)
    pt=ps.copy()
    pl=np.zeros(len(nodes),np.int32)
    for i,x in enumerate(nodes):
        for t in x['children']:
            cm[i,t]=1
        pl[i]=len(x['prefix'])
        for j,(s,t) in enumerate(x['prefix']):
            ps[i,j],pt[i,j]=s,t
    return nodes,(nv,cm,ps,pt,pl)


def hull(L,U,h,b0,b1):
    _,uk,pred=all_k_with_upper_pred(L,U,h)
    cells=count_hull(uk,b0,b1)
    for c in cells:
        c['path']=tuple(path_from_pred(pred,c['k'],L.shape[0]-1))
    return cells


def diff_lower(R,P,L,U):
    r,p=set(R),set(P)
    return (sum(int(L[e]) for e in r-p)-sum(int(U[e]) for e in p-r),
            len(R)-len(P))


def max_lines_min(lines,lo,hi):
    """Exact min of a convex max of affine lines on a closed interval."""
    points={F(lo),F(hi)}
    for i,(a,b) in enumerate(lines):
        for c,d in lines[i+1:]:
            if b != d:
                x=F(c-a,b-d)
                if lo<x<hi:
                    points.add(x)
    value,x=min((max(a+b*x for a,b in lines),x) for x in points)
    active=max(range(len(lines)),key=lambda i:lines[i][0]+lines[i][1]*x)
    return value,x,active


def retirement(R,C,L,U,b0,b1):
    lines=[diff_lower(R,c['path'],L,U) for c in C]
    gap,x,j=max_lines_min(lines,b0,b1)
    return gap,x,C[j]['path']


def actions(R,P,L,U):
    return [{'edge':e,'width':(int(U[e])-int(L[e]))/Q}
            for e in sorted(set(R)^set(P))
            if int(U[e])-int(L[e])>100]


class SafeArchive:
    def __init__(self,n,h,b0,b1,mode='archive',incremental=False):
        self.n,self.h,self.b0,self.b1=n,h,F(b0),F(b1)
        self.mode=mode
        self.A={}
        self.cover={}
        self.retired={}
        self.nodes=None
        self.previous=set()
        self.incremental=incremental
        self.upper=None

    def refresh(self,L,U):
        if self.incremental:
            if self.upper is None:
                from incremental_upper import UpperLayers
                self.upper=UpperLayers(U,self.h)
            C=self.upper.hull(self.b0,self.b1)
        else:
            C=hull(L,U,self.h,self.b0,self.b1)
        current={c['path'] for c in C}
        added=[]
        for c in C:
            p=c['path']
            if p not in self.A:
                self.A[p]=len(self.A)
                self.cover[p]=None
                added.append(self.A[p])
        migrated=[self.A[p] for p in self.previous-current]
        if self.nodes is None or added or self.mode=='fresh':
            excluded=self.A if self.mode=='archive' else current
            self.nodes,self.arrays=trie(excluded,self.n)
        self.previous=current
        return C,added,migrated

    def check(self,L,U,C):
        failures=[]
        calls=0
        passed=0
        for c in C:
            p=c['path']
            cov=self.cover[p] if self.mode=='archive' else None
            points=[]
            if cov is None:
                points=[c['lo'],c['hi']]
            else:
                if c['lo']<cov[0]: points.append(c['lo'])
                if c['hi']>cov[1]: points.append(c['hi'])
            mask=np.zeros(L.shape,np.uint8)
            for e in p: mask[e]=1
            for x in dict.fromkeys(points):
                calls+=1
                v,z,t,nxt=complement(L,U,mask,np.asarray(p,np.int64),
                    x.numerator,x.denominator,*self.arrays,self.h)
                if v>0:
                    passed+=1
                    if self.mode=='archive':
                        cov=(x,x) if cov is None else (min(cov[0],x),max(cov[1],x))
                        self.cover[p]=cov
                    continue
                R=list(self.nodes[z]['prefix'])+[(int(self.arrays[0][z]),int(t))]
                cur=int(t)
                while cur<self.n:
                    t=int(nxt[cur]); R.append((cur,t)); cur=t
                R=tuple(R)
                failures.append({'type':'complement','candidate':p,'rival':R,
                    'beta':str(x),'margin':float(F(int(v),x.denominator*Q)),
                    'actions':actions(R,p,L,U)})
        retired_now=[]
        old_active=[]
        if self.mode=='archive':
            current={c['path'] for c in C}
            for R,path_id in self.A.items():
                if R in current or R in self.retired: continue
                gap,x,P=retirement(R,C,L,U,self.b0,self.b1)
                if gap>0:
                    rec={'id':path_id,'gap':str(gap),'beta':str(x),'against_id':self.A[P]}
                    self.retired[R]=rec
                    retired_now.append(rec)
                else:
                    old_active.append(path_id)
                    failures.append({'type':'retirement','candidate':P,'rival':R,
                        'beta':str(x),'margin':float(gap/Q),'actions':actions(R,P,L,U)})
        return {'failures':failures,'oracle_calls':calls,'positive_checks':passed,
                'retired_now':retired_now,'old_active':old_active}


def fit(X,y,s,t,tol=2e-11):
    xx,yy=X[s:t],y[s:t]
    m=t-s; th=np.zeros(X.shape[1])
    def fgh(z):
        eta=xx@z; pr=expit(eta)
        f=np.logaddexp(0,eta).sum()-yy@eta+.05*m*np.dot(z,z)
        g=xx.T@(pr-yy)+.1*m*z
        H=xx.T@(xx*(pr*(1-pr))[:,None])+.1*m*np.eye(len(z))
        return f,g,H
    for it in range(100):
        f,g,H=fgh(th)
        if max(abs(g))<tol: break
        direction=np.linalg.solve(H,g); step=1.
        for ls in range(60):
            nxt=th-step*direction
            if fgh(nxt)[0]<=f-1e-4*step*(g@direction):
                th=nxt; break
            step*=.5
        else: break
    return th,{'iterations':it+1,'grad_inf':float(max(abs(fgh(th)[1])))}


def direct(L,U,X,y,Xd,yd,e,fit_tol=2e-11):
    s,t=e; start=time.process_time()
    theta,rec=fit(X,y,s,t,tol=fit_tol)
    rec['fit_cpu']=time.process_time()-start
    start=time.process_time()
    bank=ia.PrefixBank.build(Xd[s:t],yd[s:t],theta)
    build_cpu=time.process_time()-start;query_start=time.process_time()
    low,high=bank.query(0,t-s,.1,theta)
    query_cpu=time.process_time()-query_start
    l=int(ia.mul_dn(ia.dfloat(Q),low).to_integral_value(rounding=ia.ROUND_FLOOR))
    u=int(ia.mul_up(ia.dfloat(Q),high).to_integral_value(rounding=ia.ROUND_CEILING))
    ol,ou=int(L[e]),int(U[e])
    nl,nu=max(ol,l),min(ou,u)
    if nl>nu: raise ValueError(('inconsistent certificates',e))
    L[e],U[e]=nl,nu
    rec.update(edge=e,strict_cpu=time.process_time()-start,
               strict_build_cpu=build_cpu,strict_target_query_cpu=query_cpu,
               lower_gain=(int(L[e])-ol)/Q,upper_gain=(ou-int(U[e]))/Q,
               theta=list(map(float,theta)))
    return rec


def phase_audit(L,U,C,b0,b1):
    """Exact candidate-internal phase topology; complement cover is a premise."""
    transitions=[]
    for i in range(len(C)-1):
        A,B=C[i]['path'],C[i+1]['path']
        intercept,slope=diff_lower(B,A,L,U)
        upper=-diff_lower(A,B,L,U)[0]
        lo,hi=F(intercept,-slope),F(upper,-slope)
        third=[]
        for j,c in enumerate(C):
            if j in (i,i+1): continue
            val,x,_=max_lines_min([diff_lower(c['path'],A,L,U),
                                   diff_lower(c['path'],B,L,U)],lo,hi)
            third.append({'k':c['k'],'gap':str(val),'beta':str(x),'pass':val>0})
        ok=b0<=lo<=hi<=b1 and all(t['pass'] for t in third)
        transitions.append({'lo':str(lo),'hi':str(hi),
            'width_beta':str((hi-lo)/Q),'pass':ok,'third':third})
    stable=[]
    for i,c in enumerate(C):
        lo=F(b0) if i==0 else F(transitions[i-1]['hi'])
        hi=F(b1) if i==len(C)-1 else F(transitions[i]['lo'])
        ok=lo<hi
        for j,d in enumerate(C):
            if i==j: continue
            a,b=diff_lower(d['path'],c['path'],L,U)
            v,w=a+b*lo,a+b*hi
            ok=ok and min(v,w)>=0 and max(v,w)>0
        stable.append({'k':c['k'],'lo':str(lo),'hi':str(hi),'pass':bool(ok)})
    max_width=max(sum(int(U[e])-int(L[e]) for e in c['path']) for c in C)/Q
    return {'counts':[c['k'] for c in C],'transitions':transitions,'stable':stable,
            'max_retained_path_value_width':max_width,'value_epsilon':.01,
            'value_pass':max_width<=.01,
            'pass':all(x['pass'] for x in transitions+stable) and max_width<=.01}


def run(case,mode,max_actions,out,incremental=False):
    if case=='sunspots':
        raw=PACKAGE/'workspace/SAC_R328_UNPACKED/R328_SUNSPOTS_INPUT.npz'
        bounds=PACKAGE/'workspace/SAC_R329_ADAPTIVE_PROBE_20261107/R329_BROADCAST_BOUNDS.npz'
        b0,b1=Q//10,6*Q//5
    else:
        raw=PACKAGE/'workspace/SAC_R328_UNPACKED/inputs/SKAB_VALVE1_1.npz'
        bounds=PACKAGE/'workspace/SAC_R331_METHOD/R331_SKAB_INITIAL_BOUNDS.npz'
        b0,b1=2*Q,30*Q
    raw=np.load(raw); z=np.load(bounds)
    X,y=raw['X'],raw['y']; L,U=z['low_int'].copy(),z['high_int'].copy()
    # Precompile reusable kernels outside the suffix runtime boundary for both arms.
    arch=SafeArchive(len(y),8,b0,b1,mode,incremental)
    C,_,_=arch.refresh(L,U); arch.check(L,U,C)
    if incremental:
        from incremental_upper import decrease
        decrease(arch.upper.D.copy(),arch.upper.P.copy(),U,0,8,8)
    start=time.process_time()
    arch=SafeArchive(len(y),8,b0,b1,mode,incremental)
    Xd=[[ia.dfloat(v) for v in row] for row in X]; yd=list(map(int,y))
    hist=[]; snapshots=[]; selected=set(); graph_cpu=0.; selection_cpu=0.
    for step in range(max_actions+1):
        g=time.process_time()
        C,new,migration=arch.refresh(L,U)
        state=arch.check(L,U,C)
        graph_cpu+=time.process_time()-g
        snap={k:v for k,v in state.items() if k!='failures'}
        snap.update(step=step,counts=[c['k'] for c in C],new=new,migration=migration,
                    failed=len(state['failures']),archive=len(arch.A))
        snapshots.append(snap)
        print(case,mode,step,'failed',snap['failed'],'archive',len(arch.A),
              'active_old',snap['old_active'],'retired',len(arch.retired),flush=True)
        if not state['failures'] or step==max_actions: break
        g=time.process_time();scores={}
        for f in state['failures']:
            for a in f['actions']:
                e=tuple(a['edge'])
                if e in selected: continue
                scores[e]=scores.get(e,0)+min(1.,a['width']/max(-f['margin'],1e-9))
        if not scores: break
        e=max(scores,key=lambda e:(scores[e],1/(1+e[1]-e[0])))
        selected.add(e); selection_cpu+=time.process_time()-g
        rec=direct(L,U,X,y,Xd,yd,e); rec['score']=scores[e];hist.append(rec)
        if incremental and rec['upper_gain']>0:
            g=time.process_time();arch.upper.update(*e);graph_cpu+=time.process_time()-g
    cover=not state['failures']
    g=time.process_time();phase=phase_audit(L,U,C,b0,b1)
    from path_support import all_supports
    supports=all_supports(arch.A,L,U,b0,b1)
    final_cpu=time.process_time()-g
    total=time.process_time()-start
    result={'case':case,'mode':mode,'status':'COMPLETE_PHASE_PASS' if cover and phase['pass']
        else 'UNRESOLVED','cover_pass':cover,'actions':len(hist),'total_suffix_cpu':total,
        'graph_cpu':graph_cpu,'selection_cpu':selection_cpu,'final_phase_cpu':final_cpu,
        'fit_cpu':sum(r['fit_cpu'] for r in hist),'strict_cpu':sum(r['strict_cpu'] for r in hist),
        'history':hist,'snapshots':snapshots,'phase':phase,
        'archive_paths':[list(p) for p in arch.A],
        'retirement_certificates':list(arch.retired.values()),
        'incremental_discovery':incremental,
        'upper_updates':None if arch.upper is None else {'updates':arch.upper.updates,
            'changed_states':arch.upper.states,'edge_relaxations':arch.upper.relaxations},
        'final_failed':state['failures'],
        'scope':'online suffix from frozen certified-cost snapshot; startup evidence excluded; final graph phase audit included; JIT precompiled'}
    result['complete_possible_path_supports']=supports
    Path(out).write_text(json.dumps(result,indent=2),encoding='utf8')
    np.savez_compressed(Path(out).with_suffix('.npz'),low_int=L,high_int=U)
    print('RESULT',case,mode,result['status'],len(hist),total,flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--case',choices=['sunspots','skab'],required=True)
    p.add_argument('--mode',choices=['fresh','archive'],required=True)
    p.add_argument('--max-actions',type=int,default=160)
    p.add_argument('--out',required=True)
    p.add_argument('--incremental',action='store_true')
    a=p.parse_args();run(a.case,a.mode,a.max_actions,a.out,a.incremental)
