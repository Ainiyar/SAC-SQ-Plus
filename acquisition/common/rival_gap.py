"""Two distinct non-archive rival values for certified action-value screening."""
import numpy as np
from numba import njit
from safe_archive import INF

@njit(cache=True)
def two_rivals(L,U,mask,edges,bn,bd,nv,cm,ps,pt,pl,h):
    n=len(L)-1;first=np.full(n+1,INF,np.int64);second=first.copy();first[n]=0
    for s in range(n-h,-1,-1):
        for t in range(s+h,n+1):
            w=U[s,t] if mask[s,t] else L[s,t]
            for rank in range(2):
                suffix=first[t] if rank==0 else second[t]
                if suffix==INF:continue
                val=w*bd+bn+suffix
                if val<first[s]:second[s]=first[s];first[s]=val
                elif val<second[s]:second[s]=val
    one=INF;two=INF
    for z in range(len(nv)):
        s=nv[z];prefix=np.int64(0)
        for j in range(pl[z]):
            a,b=ps[z,j],pt[z,j]
            w=U[a,b] if mask[a,b] else L[a,b]
            prefix+=w*bd+bn
        for t in range(s+h,n+1):
            if cm[z,t]:continue
            w=U[s,t] if mask[s,t] else L[s,t]
            for rank in range(2):
                suffix=first[t] if rank==0 else second[t]
                if suffix==INF:continue
                val=prefix+w*bd+bn+suffix
                if val<one:two=one;one=val
                elif val<two:two=val
    cand=np.int64(0)
    for s,t in edges:cand+=U[s,t]*bd+bn
    return (INF if one==INF else one-cand),(INF if two==INF else two-cand)
