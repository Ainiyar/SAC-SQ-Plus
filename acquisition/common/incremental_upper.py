"""Decrease-only count-layer shortest paths on a segmentation DAG.
An edge decrease seeds affected states and propagates only changed distances
or canonical predecessors. Canonical ties use the smallest predecessor.
"""
import numpy as np
from numba import njit
from safe_archive import INF, count_hull, path_from_pred

@njit(cache=True)
def initialize(U,h):
    n=len(U)-1;K=n//h
    D=np.full((K+1,n+1),INF,np.int64)
    P=np.full((K+1,n+1),-1,np.int32)
    D[0,0]=0
    for k in range(1,K+1):
        for t in range(k*h,n+1):
            for s in range((k-1)*h,t-h+1):
                if D[k-1,s]==INF:continue
                val=D[k-1,s]+U[s,t]
                if val<D[k,t]: D[k,t],P[k,t]=val,s
    return D,P

@njit(cache=True)
def decrease(D,P,U,s,t,h):
    n=len(U)-1;K=n//h
    changed=np.zeros(D.shape,np.uint8)
    states=0;relaxations=0
    # The source node s is upstream of this edge; all prefix costs at s persist.
    for k in range(1,K+1):
        if D[k-1,s]==INF:continue
        val=D[k-1,s]+U[s,t]
        if val<D[k,t]:
            D[k,t],P[k,t]=val,s;changed[k,t]=1
        elif val==D[k,t] and (P[k,t]<0 or s<P[k,t]):
            P[k,t]=s
    for k in range(1,K+1):
        for v in range(k*h,n+1):
            if not changed[k,v]: continue
            states+=1
            if k==K:continue
            for w in range(v+h,n+1):
                relaxations+=1
                val=D[k,v]+U[v,w]
                if val<D[k+1,w]:
                    D[k+1,w],P[k+1,w]=val,v;changed[k+1,w]=1
                elif val==D[k+1,w] and (P[k+1,w]<0 or v<P[k+1,w]):
                    P[k+1,w]=v
    return states,relaxations


class UpperLayers:
    def __init__(self,U,h):
        self.U=U;self.h=h
        self.D,self.P=initialize(U,h)
        self.updates=0;self.states=0;self.relaxations=0
    def update(self,s,t):
        a,b=decrease(self.D,self.P,self.U,s,t,self.h)
        self.updates+=1;self.states+=a;self.relaxations+=b
    def hull(self,b0,b1):
        C=count_hull(self.D[:,-1],b0,b1)
        for c in C:c['path']=tuple(path_from_pred(self.P,c['k'],len(self.U)-1))
        return C
