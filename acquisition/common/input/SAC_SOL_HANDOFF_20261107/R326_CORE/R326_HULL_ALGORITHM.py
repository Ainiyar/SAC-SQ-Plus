"""O(K) exact-rational lower envelope of affine count lines Uk+k*beta.
Line slopes added in decreasing order. Return segments of positive length within [b0,b1].
This algorithm discovers the same open phase cells as exhaustive pair-intersections.
"""
from fractions import Fraction as F

def count_hull(uk,low,high):
 S=[]
 for k in range(len(uk)-1,0,-1):
  b=int(uk[k]);start=None
  while S:
   prev=S[-1]
   x=F(b-prev['b'],prev['k']-k)
   if prev['start'] is not None and x<=prev['start']:
    S.pop()
   else:
    start=x;break
  S.append({'k':k,'b':b,'start':start})
 out=[]
 for i,c in enumerate(S):
  left=F(low) if c['start'] is None else max(F(low),c['start'])
  right=F(high) if i==len(S)-1 else min(F(high),S[i+1]['start'])
  if left<right:out.append({'k':c['k'],'lo':left,'hi':right})
 return out
