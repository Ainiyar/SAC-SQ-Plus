"""Exact superadditive segment lower closure, existing SAC mathematical asset."""
from numba import njit

@njit(cache=True)
def close_lower(L,h):
    n=len(L)-1;changed=0
    for length in range(2*h,n+1):
        for s in range(n-length+1):
            t=s+length;best=L[s,t]
            for r in range(s+h,t-h+1):
                val=L[s,r]+L[r,t]
                if val>best:best=val
            if best>L[s,t]:L[s,t]=best;changed+=1
    return changed
