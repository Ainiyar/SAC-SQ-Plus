"""Certified source-pool adaptation of Reliever's deterministic relief intervals.
The statistical proxy is turned into target-valid directed L/U evidence; the
same full proof and refinement engine then determines completion.
"""
import time
import numpy as np
from scipy.special import expit
from safe_archive import ia,Q,fit

def relief_intervals(n,h=8,r=.5):
    b=r**(-.5);w=b-1;out=set()
    max_layer=int(np.floor(np.log((1+w)*n/h)/np.log(b)))
    for k in range(max_layer+1):
        length=b**k*h/(1+w);shift=w*length
        count=int(np.floor((n-length)/shift));offset=n/2-(length+count*shift)/2
        m=max(1,min(n,int(np.floor(length))))
        for q in range(count+1):
            s=max(0,min(n-m,int(round(q*shift+offset))))
            out.add((s,s+m))
        out.add((0,m));out.add((n-m,n))
    return sorted(out,key=lambda e:(e[1]-e[0],e[0]))

def initialize(L,U,X,y,Xd,yd,r=.5,sources=None):
    n=len(y);custom=sources is not None
    if sources is None:sources=relief_intervals(n,r=r)
    ss,tt=np.triu_indices(n+1,8);m=tt-ss
    bestL=np.zeros(len(ss));bestU=U[ss,tt].astype(float)/Q
    li=np.full(len(ss),-1,np.int32);ui=li.copy();records=[];thetas=[]
    guide_start=time.process_time();fit_total=0.
    for idx,e in enumerate(sources):
        t0=time.process_time()
        if e is None:theta=np.zeros(X.shape[1]);info={}
        else:theta,info=fit(X,y,*e,tol=1e-7)
        fitcpu=time.process_time()-t0;fit_total+=fitcpu;thetas.append(theta)
        z=X@theta;p=expit(z)
        fp=np.r_[0,np.cumsum(np.logaddexp(0,z)-y*z)]
        gp=np.vstack([np.zeros((1,X.shape[1])),np.cumsum(X*(p-y)[:,None],axis=0)])
        f=fp[tt]-fp[ss]+.05*m*float(theta@theta)
        g=gp[tt]-gp[ss]+.1*m[:,None]*theta
        lower=np.maximum(0,f-np.sum(g*g,axis=1)/(.2*m))
        better=lower>bestL;bestL[better]=lower[better];li[better]=idx
        better=f<bestU;bestU[better]=f[better];ui[better]=idx
        records.append({'source':e,'fit_cpu':fitcpu,'bank_cpu':0.,'strict_query_cpu':0.,
                        'queries':0,'theta':list(map(float,theta)),**info})
    guidecpu=time.process_time()-guide_start-fit_total
    used=sorted(set(int(x) for x in np.r_[li,ui] if x>=0));banks={}
    for idx in used:
        t0=time.process_time();banks[idx]=ia.PrefixBank.build(Xd,yd,thetas[idx])
        records[idx]['bank_cpu']=time.process_time()-t0
    for j,(s,t) in enumerate(zip(ss,tt)):
        s,t=int(s),int(t)
        for idx in set([int(li[j]),int(ui[j])])-{ -1 }:
            t0=time.process_time();lo,up=banks[idx].query(s,t,.1,thetas[idx])
            l=int(ia.mul_dn(ia.Decimal(Q),lo).to_integral_value(rounding=ia.ROUND_FLOOR))
            u=int(ia.mul_up(ia.Decimal(Q),up).to_integral_value(rounding=ia.ROUND_CEILING))
            nl,nu=max(int(L[s,t]),l),min(int(U[s,t]),u)
            if nl>nu:raise ValueError('inconsistent certified relief evidence')
            L[s,t],U[s,t]=nl,nu
            records[idx]['strict_query_cpu']+=time.process_time()-t0
            records[idx]['queries']+=1
    return {'banks':records,'source_policy':('four seed anchors; best lower/upper target-valid strict queries' if custom else 'Reliever-inspired deterministic relief intervals; best lower/upper target-valid strict queries'),
            'relief_r':r,'candidate_relief_sources':len(sources),'materialized_banks':len(used),
            'relief_guide_cpu':guidecpu}
