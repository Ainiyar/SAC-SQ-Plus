"""One pass over count-layers, retaining upper predecessor to avoid re-solving every active count."""
import numpy as np
from numba import njit
INF=np.int64(1<<60)
@njit(cache=True)
def all_k_with_upper_pred(lo,hi,h):
 n=lo.shape[0]-1;K=n//h
 L=np.full((K+1,n+1),INF,np.int64)
 U=np.full((K+1,n+1),INF,np.int64)
 P=np.full((K+1,n+1),-1,np.int32)
 L[0,0]=0;U[0,0]=0
 for k in range(1,K+1):
  for t in range(k*h,n+1):
   bl=INF;bu=INF;ps=-1
   for s in range((k-1)*h,t-h+1):
    if L[k-1,s]!=INF:
     v=L[k-1,s]+lo[s,t]
     if v<bl:bl=v
    if U[k-1,s]!=INF:
     v=U[k-1,s]+hi[s,t]
     if v<bu:bu=v;ps=s
   L[k,t]=bl;U[k,t]=bu;P[k,t]=ps
 return L[:,n],U[:,n],P

def path_from_pred(pred,k,n):
 out=[];t=n
 for j in range(k,0,-1):
  s=int(pred[j,t]);out.append((s,t));t=s
 if t!=0:raise ValueError('path does not end at 0')
 return out[::-1]
