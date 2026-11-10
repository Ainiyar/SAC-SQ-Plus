"""Stable all-nonarchive-path exponential margin sum on the segment DAG.

Advisory floating arithmetic. Strict certification uses the integer oracle.
Every outside path has one unique first departure from the archive trie.
"""
import numpy as np
from numba import njit
from safe_archive import Q

@njit(cache=True)
def log_partition(L,U,mask,edges,beta,tau,nv,cm,ps,pt,pl,h):
    n=len(L)-1;suffix=np.full(n+1,-np.inf);suffix[n]=0.
    for s in range(n-h,-1,-1):
        for t in range(s+h,n+1):
            w=U[s,t] if mask[s,t] else L[s,t]
            term=-(w/Q+beta)/tau+suffix[t]
            suffix[s]=np.logaddexp(suffix[s],term)
    total=-np.inf
    for z in range(len(nv)):
        s=nv[z];prefix=0.
        for j in range(pl[z]):
            a,b=ps[z,j],pt[z,j]
            w=U[a,b] if mask[a,b] else L[a,b]
            prefix-=(w/Q+beta)/tau
        for t in range(s+h,n+1):
            if cm[z,t]:continue
            w=U[s,t] if mask[s,t] else L[s,t]
            term=prefix-(w/Q+beta)/tau+suffix[t]
            total=np.logaddexp(total,term)
    candidate=0.
    for s,t in edges:candidate+=(U[s,t]/Q+beta)/tau
    return total+candidate


@njit(cache=True)
def log_partition_fast(L,U,mask,edges,beta,tau,nv,cm,ps,pt,pl,h):
    # Reweight with exact-in-real-arithmetic shortest suffix potentials. All
    # exponential arguments are nonpositive; scaled mass is at most path count.
    n=len(L)-1
    if n*np.log(2.)>650.:
        return log_partition(L,U,mask,edges,beta,tau,nv,cm,ps,pt,pl,h)
    distance=np.full(n+1,np.inf);distance[n]=0.
    mass=np.zeros(n+1);mass[n]=1.
    for s in range(n-h,-1,-1):
        for t in range(s+h,n+1):
            w=U[s,t] if mask[s,t] else L[s,t]
            distance[s]=min(distance[s],w/Q+beta+distance[t])
        if not np.isfinite(distance[s]):continue
        for t in range(s+h,n+1):
            if mass[t]==0.:continue
            w=U[s,t] if mask[s,t] else L[s,t]
            reduced=w/Q+beta+distance[t]-distance[s]
            mass[s]+=np.exp(-max(0.,reduced)/tau)*mass[t]
    prefixes=np.zeros(len(nv));best=np.inf
    for z in range(len(nv)):
        for j in range(pl[z]):
            a,b=ps[z,j],pt[z,j]
            w=U[a,b] if mask[a,b] else L[a,b]
            prefixes[z]+=w/Q+beta
        s=nv[z]
        for t in range(s+h,n+1):
            if cm[z,t] or mass[t]==0.:continue
            w=U[s,t] if mask[s,t] else L[s,t]
            best=min(best,prefixes[z]+w/Q+beta+distance[t])
    if not np.isfinite(best):return -np.inf
    total=0.
    for z in range(len(nv)):
        s=nv[z]
        for t in range(s+h,n+1):
            if cm[z,t] or mass[t]==0.:continue
            w=U[s,t] if mask[s,t] else L[s,t]
            reduced=prefixes[z]+w/Q+beta+distance[t]-best
            total+=np.exp(-max(0.,reduced)/tau)*mass[t]
    candidate=0.
    for s,t in edges:candidate+=U[s,t]/Q+beta
    return (candidate-best)/tau+np.log(total)
