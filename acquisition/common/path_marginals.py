"""All nonarchive path edge marginals via unique-departure forward/backward."""
import numpy as np
from numba import njit
from safe_archive import Q

@njit(cache=True)
def partition_marginals(L,U,mask,edges,beta,tau,nv,cm,ps,pt,pl,h):
    n=len(L)-1
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
            mass[s]+=np.exp(-max(0.,w/Q+beta+distance[t]-distance[s])/tau)*mass[t]
    prefixes=np.zeros(len(nv));best=np.inf
    for z in range(len(nv)):
        for j in range(pl[z]):
            a,b=ps[z,j],pt[z,j];w=U[a,b] if mask[a,b] else L[a,b]
            prefixes[z]+=w/Q+beta
        for t in range(nv[z]+h,n+1):
            if cm[z,t] or mass[t]==0.:continue
            w=U[nv[z],t] if mask[nv[z],t] else L[nv[z],t]
            best=min(best,prefixes[z]+w/Q+beta+distance[t])
    edge_mass=np.zeros(L.shape);arrivals=np.zeros(n+1);total=0.
    if not np.isfinite(best):return -np.inf,edge_mass
    for z in range(len(nv)):
        s=nv[z];prefix_mass=0.
        for t in range(s+h,n+1):
            if cm[z,t] or mass[t]==0.:continue
            w=U[s,t] if mask[s,t] else L[s,t]
            weight=np.exp(-max(0.,prefixes[z]+w/Q+beta+distance[t]-best)/tau)*mass[t]
            total+=weight;prefix_mass+=weight;edge_mass[s,t]+=weight;arrivals[t]+=weight
        for j in range(pl[z]):edge_mass[ps[z,j],pt[z,j]]+=prefix_mass
    edge_mass/=total;arrivals/=total
    for s in range(n-h+1):
        if arrivals[s]==0. or mass[s]==0.:continue
        for t in range(s+h,n+1):
            if mass[t]==0.:continue
            w=U[s,t] if mask[s,t] else L[s,t]
            probability=np.exp(-max(0.,w/Q+beta+distance[t]-distance[s])/tau)*mass[t]/mass[s]
            flow=arrivals[s]*probability
            edge_mass[s,t]+=flow;arrivals[t]+=flow
    candidate=0.
    for s,t in edges:candidate+=U[s,t]/Q+beta
    return (candidate-best)/tau+np.log(total),edge_mass
