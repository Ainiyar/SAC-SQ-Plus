"""Lazy, cached, target-valid Reliever source policy with exact recovery.

Same proof-debt target and materialization controller as SAC. The source is the
longest deterministic relief interval contained in the requested target. A
second request for a target uses its own fitted parameter. Both numerical fits
and full directed PrefixBanks persist and are charged only when constructed.
"""
import time
import numpy as np
from relief_initial import relief_intervals
from proof_lookahead import acquire

def geometry(n,h=8,r=.9):
    sources=relief_intervals(n,h,r);best=np.zeros((n+1,n+1),np.int32)
    index=np.full(best.shape,-1,np.int32)
    for i,(s,t) in enumerate(sources):
        length=t-s;view=best[:s+1,t:];iv=index[:s+1,t:]
        better=view<length;view[better]=length;iv[better]=i
    ss,tt=np.triu_indices(n+1,h)
    return sources,index,float(np.min(best[ss,tt]/(tt-ss)))

class LazyRelief:
    def __init__(self,X,y,Xd,yd,costs,r=.9,pricing='entropy',temperature='debt',
                 target_priority='entropy',closure_kernel='incremental',theta_cache=None):
        t=time.process_time();self.sources,self.index,self.coverage=geometry(len(y),8,r)
        self.geometry_cpu=time.process_time()-t
        self.X,self.y,self.Xd,self.yd=X,y,Xd,yd;self.costs=costs;self.r=r
        self.pricing,self.temperature=pricing,temperature
        self.target_priority,self.closure_kernel=target_priority,closure_kernel
        self.theta_cache={} if theta_cache is None else theta_cache
        self.bank_cache={};self.tried=set();self.proxy_actions=0;self.exact_actions=0
    def __call__(self,arch,L,U,target,failures):
        target=tuple(target);idx=int(self.index[target])
        if target not in self.tried and idx>=0:
            source=self.sources[idx];self.tried.add(target);self.proxy_actions+=1
            policy='MAX_CONTAINED_RELIEF'
        else:
            source=target;self.exact_actions+=1;policy='EXACT_TARGET_RECOVERY'
        rec=acquire(arch,L,U,self.X,self.y,self.Xd,self.yd,target,failures,*self.costs,True,
                    pricing=self.pricing,temperature=self.temperature,source_edges=[source],
                    theta_cache=self.theta_cache,target_priority=self.target_priority,
                    closure_kernel=self.closure_kernel,source_targets={source:target},bank_cache=self.bank_cache)
        rec.update(source_policy=policy,target_complete=source==target,
                   proxy_coverage=(source[1]-source[0])/(target[1]-target[0]))
        return rec
    def summary(self):
        return {'policy':'lazy max-contained relief, cached finite parameters and directed PrefixBanks, exact recovery',
                'coverage_r':self.r,'candidate_sources':len(self.sources),
                'actual_minimum_coverage':self.coverage,'geometry_cpu':self.geometry_cpu,
                'materialized_full_banks':len(self.bank_cache),'cached_parameter_count':len(self.theta_cache),
                'proxy_actions':self.proxy_actions,'exact_recovery_actions':self.exact_actions}
