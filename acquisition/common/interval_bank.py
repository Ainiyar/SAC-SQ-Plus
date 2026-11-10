"""Directed source evidence restricted to the union of requested intervals."""
from bisect import bisect_right
from safe_archive import ia

def interval_union(edges):
    out=[]
    for s,t in sorted(set((int(s),int(t)) for s,t in edges)):
        if out and s<=out[-1][1]:out[-1]=(out[-1][0],max(t,out[-1][1]))
        else:out.append((s,t))
    return out

class UnionBank:
    def __init__(self,Xd,yd,theta,edges):
        self.components=interval_union(edges)
        self.starts=[s for s,t in self.components]
        self.banks=[ia.PrefixBank.build(Xd[s:t],yd[s:t],theta) for s,t in self.components]
        self.rows_built=sum(t-s for s,t in self.components)

    def query(self,s,t,lam,theta):
        j=bisect_right(self.starts,s)-1
        if j<0 or t>self.components[j][1]:raise ValueError('query outside source evidence footprint')
        offset=self.components[j][0]
        return self.banks[j].query(s-offset,t-offset,lam,theta)
