"""First-order salience of every exchange edge under all-rival proof debt."""
from fractions import Fraction as F
import numpy as np
from safe_archive import Q
from path_partition import log_partition_fast
from path_marginals import partition_marginals

def salience(arch,L,U,failures,temperature='debt',h=8):
    wl=np.zeros(L.shape);wu=np.zeros(L.shape);calls=0
    logcount=1.
    if temperature=='entropy_gap':
        zero=np.zeros_like(L)
        logcount=log_partition_fast(zero,zero,np.zeros(L.shape,np.uint8),
            np.empty((0,2),np.int64),0.,1.,*arch.arrays,h)
    for f in failures:
        P=tuple(map(tuple,f['candidate']));R=tuple(map(tuple,f['rival']))
        tau=max(-f['margin'],1e-3)/(max(logcount,1.) if temperature=='entropy_gap' else 1.)
        if f['type']=='complement':
            mask=np.zeros(L.shape,np.uint8)
            for e in P:mask[e]=1
            _,prob=partition_marginals(L,U,mask,np.asarray(P,np.int64),
                float(F(f['beta']))/Q,tau,*arch.arrays,h)
            prob=np.clip(prob,0.,1.)
            wl+=(1-mask)*prob/tau;wu+=mask*(1-prob)/tau;calls+=1
        else:
            for e in set(R)-set(P):wl[e]+=1/tau
            for e in set(P)-set(R):wu[e]+=1/tau
    return (wl,wu),{'marginal_distribution_calls':calls,'target_priority':'entropy'}
