"""Forecast the certified relevance projection after a monotone evidence menu.

The pre-action reduced graph contains every possible future nonpositive path.
Re-pruning this graph yields the exact post-action relevance operator. Its mass
is zero exactly at external proof completion, and dominates the progress of the
frozen-graph submodular surrogate. Experimental pricing module, separate from
the frozen whole-round main controller.
"""
import numpy as np
from safe_archive import Q,diff_lower
from proof_relevance import exact_prune,partition_core

def projected_logmass(item,L,U):
    beta=item['beta'];graph=item['graph']
    if graph is None:
        a,b=diff_lower(item['R'],item['P'],L,U);value=a+b*beta
        return (-np.inf if value>0 else -float(value/Q)/item['tau']),value>0,0
    ptr,dest,pos,order,accept=graph
    pp,dd,gap,cut=exact_prune(L,U,item['mask'],item['edges'],beta.numerator,beta.denominator,*graph)
    if gap>0:return -np.inf,True,len(dd)
    value,*_=partition_core(L,U,item['mask'],item['edges'],float(beta)/Q,item['tau'],pp,dd,pos,order,accept)
    return value,False,len(dd)

def price(context,L,U,alternatives):
    records=[];calls=0
    for alt in alternatives:
        gain=0.;complete=arcs=0
        for item in context['items']:
            current,done,remaining=projected_logmass(item,alt['L'],alt['U'])
            gain+=item['weight']*(-np.expm1(min(0.,current-item['initial'])))
            complete+=done;arcs+=remaining;calls+=item['graph'] is not None
        records.append({'gain':float(gain),'score':float(gain/alt['cost']),
                        'predicted_complete_obligations':int(complete),'projected_arcs':arcs,'evaluated':True})
    best=max(range(len(records)),key=lambda i:(records[i]['score'],-alternatives[i]['cost']))
    return best,records,{'pricing':'projected','forecast_relevance_projections':calls,**context['stats']}
