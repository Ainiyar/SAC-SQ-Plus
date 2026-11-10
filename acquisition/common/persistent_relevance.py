"""Exact repeated relevance pruning on the surviving automaton only.

For fixed P,beta and outside-archive language, every nonpositive path survives
all prior masks. Therefore pruning a prior mask gives exactly the full-graph
mask at the new monotone bounds. Archive growth starts a new language epoch.
"""
from fractions import Fraction as F
import numpy as np
from safe_archive import Q,INF
from proof_relevance import topology,exact_prune

def begin_epoch(arch,failures,h):
    version=(len(arch.A),arch.n if hasattr(arch,'n') else None,h)
    if getattr(arch,'_relevance_survivor_version',None)!=version:
        arch._relevance_survivors={};arch._relevance_survivor_version=version
    active={(tuple(map(tuple,f['candidate'])),F(f['beta'])) for f in failures if f['type']=='complement'}
    arch._relevance_survivors={k:v for k,v in arch._relevance_survivors.items() if k in active}

def pruned_graph_cached(arch,L,U,P,beta,h=8):
    P=tuple(P);beta=F(beta);key=(P,beta);full=topology(arch,len(L)-1,h)
    saved=arch._relevance_survivors.get(key);graph=full if saved is None else saved[0]
    mask=np.zeros(L.shape,np.uint8)
    for edge in P:mask[edge]=1
    pp,dd,gap,cutgap=exact_prune(L,U,mask,np.asarray(P,np.int64),beta.numerator,beta.denominator,*graph)
    if saved is not None:cutgap=min(cutgap,saved[1])
    reduced=(pp,dd,*full[2:]);arch._relevance_survivors[key]=(reduced,cutgap)
    stats={'total_arcs':len(full[1]),'retained_arcs':len(dd),'pruning_input_arcs':len(graph[1]),
           'relevance_cache_hit':saved is not None,
           # The retained pre-pruning minimum is the global minimum whenever
           # it is nonpositive. Context preparation is only for failed debts.
           'exact_min_margin':None if gap==INF or gap>0 else int(gap)/(beta.denominator*Q),
           'minimum_removed_arc_margin':None if cutgap==INF else int(cutgap)/(beta.denominator*Q)}
    return reduced,mask,stats
