"""Exact monotone max-plus closure updates for evidence and action forecasts.

A closed lower matrix stores the longest legal subpartition value. Raising one
edge(s,t) changes only intervals containing(s,t), via one max-plus outer update.
Dense batches use the existing full closure; sparse batches avoid cubic work.
"""
import numpy as np
from numba import njit
from lower_closure import close_lower

@njit(cache=True)
def sparse_raise(D,starts,ends,values,h):
    n=len(D)-1;relaxations=0;propagated=0
    for z in range(len(starts)):
        s,t,w=starts[z],ends[z],values[z]
        if w<=D[s,t]:continue
        D[s,t]=w;propagated+=1
        for i in range(s+1):
            if i!=s and s-i<h:continue
            left=0 if i==s else D[i,s]
            for j in range(t,n+1):
                if j!=t and j-t<h:continue
                right=0 if j==t else D[t,j]
                v=left+w+right;relaxations+=1
                if v>D[i,j]:D[i,j]=v
    return relaxations,propagated

def close_updates(base,starts,ends,values,h=8,force_sparse=False):
    """Return exact closure after monotone updates to an already closed base."""
    starts=np.asarray(starts,np.int64);ends=np.asarray(ends,np.int64)
    values=np.asarray(values,np.int64)
    changed=values>base[starts,ends]
    starts,ends,values=starts[changed],ends[changed],values[changed]
    D=base.copy();n=len(D)-1;m=max(0,n-2*h+1)
    full_work=m*(m+1)*(m+2)//6
    sparse_work=int(np.sum((np.maximum(0,starts-h+1)+1)*
                           (np.maximum(0,n-ends-h+1)+1)))
    if not len(starts):
        return D,{'method':'unchanged','updates':0,'estimated_work':0,'relaxations':0}
    if force_sparse or sparse_work<full_work:
        order=np.argsort(ends-starts,kind='stable')
        relax,prop=sparse_raise(D,starts[order],ends[order],values[order],h)
        return D,{'method':'sparse_outer','updates':len(starts),'estimated_work':sparse_work,
                  'relaxations':int(relax),'propagated_updates':int(prop)}
    D[starts,ends]=values;close_lower(D,h)
    return D,{'method':'full','updates':len(starts),'estimated_work':full_work,'relaxations':full_work}
