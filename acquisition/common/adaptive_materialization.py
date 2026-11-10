"""Choose local or shared strict materialization after fitting a live blocker.
Float predictions only schedule work; every published endpoint is directed.
"""
import time
import numpy as np
from scipy.special import expit
from safe_archive import Q,ia,fit

def predicted_gains(X,y,theta,L,U,targets):
    z=X@theta;p=expit(z)
    fl=np.r_[0,np.cumsum(np.logaddexp(0,z)-y*z)]
    gg=np.vstack([np.zeros((1,X.shape[1])),np.cumsum(X*(p-y)[:,None],axis=0)])
    norm=float(theta@theta);out={}
    for s,t in targets:
        m=t-s;f=fl[t]-fl[s]+.05*m*norm
        g=gg[t]-gg[s]+.1*m*theta
        lower=max(0.,f-float(g@g)/(.2*m))
        out[s,t]=(max(0.,lower-int(L[s,t])/Q),max(0.,int(U[s,t])/Q-f))
    return out

def utility(failures,gains,targets):
    selected=set(targets);value=0.
    for f in failures:
        R,P=set(f['rival']),set(f['candidate'])
        progress=sum(gains[e][0] for e in (R-P)&selected)
        progress+=sum(gains[e][1] for e in (P-R)&selected)
        value+=min(1.,progress/max(-f['margin'],1e-9))
    return value

def acquire(L,U,X,y,Xd,yd,e,failures,row_cost,query_cost,bank_allowed):
    t0=time.process_time();theta,conv=fit(X,y,*e,tol=1e-7)
    fitcpu=time.process_time()-t0;t0=time.process_time()
    targets={tuple(a['edge']) for f in failures for a in f['actions']}
    targets.add(e)
    gains=predicted_gains(X,y,theta,L,U,targets)
    score={edge:utility(failures,gains,[edge]) for edge in targets}
    ranked=sorted(targets,key=lambda x:(score[x],-x[1]+x[0]),reverse=True)
    direct_utility=utility(failures,gains,[e])
    directcost=(e[1]-e[0])*row_cost+query_cost
    candidates=[(direct_utility/directcost,'direct',[e],direct_utility,directcost)]
    if bank_allowed:
        for count in [8,32,128,256]:
            chosen=list(dict.fromkeys([e]+ranked[:count]))
            chosen=[edge for edge in chosen if edge==e or score[edge]>1e-8]
            val=utility(failures,gains,chosen)
            cost=len(y)*row_cost+len(chosen)*query_cost
            candidates.append((val/cost,'bank',chosen,val,cost))
    _,kind,chosen,value,cost=max(candidates,key=lambda x:x[0])
    guidecpu=time.process_time()-t0;t0=time.process_time()
    if kind=='bank':
        bank=ia.PrefixBank.build(Xd,yd,theta);offset=0
    else:
        bank=ia.PrefixBank.build(Xd[e[0]:e[1]],yd[e[0]:e[1]],theta);offset=e[0]
    buildcpu=time.process_time()-t0;t0=time.process_time();rows=[];upper_changed=[]
    for s,t in chosen:
        low,high=bank.query(s-offset,t-offset,.1,theta)
        l=int(ia.mul_dn(ia.Decimal(Q),low).to_integral_value(rounding=ia.ROUND_FLOOR))
        u=int(ia.mul_up(ia.Decimal(Q),high).to_integral_value(rounding=ia.ROUND_CEILING))
        oldl,oldu=int(L[s,t]),int(U[s,t])
        nl,nu=max(oldl,l),min(oldu,u)
        if nl>nu:raise ValueError('invalid strict intersection')
        L[s,t],U[s,t]=nl,nu
        dl,du=(int(L[s,t])-oldl)/Q,(oldu-int(U[s,t]))/Q
        if du>0:upper_changed.append((s,t))
        if dl>0 or du>0:rows.append({'edge':[s,t],'lower_gain':dl,'upper_gain':du})
    querycpu=time.process_time()-t0
    return {'source_edge':e,'kind':kind,'strict_queries':len(chosen),'fit_cpu':fitcpu,
            'guide_cpu':guidecpu,'bank_build_cpu':buildcpu,'strict_query_cpu':querycpu,
            'predicted_utility':value,'predicted_materialization_cpu':cost,
            'direct_predicted_utility':direct_utility,'direct_predicted_cpu':directcost,
            'improved':rows,'upper_changed':upper_changed,'theta':list(map(float,theta)),**conv}
