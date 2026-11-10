"""Lazy whole-envelope action pricing using two-rival gap bounds.

The bounds control advisory action selection. Directed evidence publication is
independent of this module. Full and screened modes maximize the same forecast.
"""
from fractions import Fraction as F
import numpy as np
from safe_archive import Q, INF, complement, diff_lower
from rival_gap import two_rivals


def price(arch,L,U,failures,alternatives,mode='full',h=8,temperature='entropy_gap',relevance_context=None):
    # alternatives contain monotone forecast matrices, after lower closure.
    if mode in ['relevant','automaton','persistent','relevant_hard','projected']:
        from proof_relevance import prepare_context,price_context,price_hard_context
        if relevance_context is None:
            relevance_context=prepare_context(arch,L,U,failures,temperature,h,False,'persistent' if mode in ['persistent','relevant_hard','projected'] else mode=='relevant')
        if mode=='projected':
            from projected_proof_debt import price as projected_price
            return projected_price(relevance_context,L,U,alternatives)
        return (price_hard_context if mode=='relevant_hard' else price_context)(relevance_context,L,U,alternatives)
    prepared=[];gap_queries=0;entropy_calls=0;entropy=[]
    if mode=='entropy':
        from path_partition import log_partition_fast as log_partition
        zero=np.zeros_like(L);mask0=np.zeros(L.shape,np.uint8)
        logcount=log_partition(zero,zero,mask0,np.empty((0,2),np.int64),0.,1.,*arch.arrays,h)
    for f in failures:
        P=tuple(map(tuple,f['candidate']));R=tuple(map(tuple,f['rival']))
        beta=F(f['beta']);mask=np.zeros(L.shape,np.uint8)
        for edge in P:mask[edge]=1
        gap=0.
        if mode=='gap' and f['type']=='complement':
            m1,m2=two_rivals(L,U,mask,np.asarray(P,np.int64),
                            beta.numerator,beta.denominator,*arch.arrays,h)
            gap_queries+=1
            gap=float('inf') if m2==INF else (int(m2)-int(m1))/(beta.denominator*Q)
        prepared.append((f,P,R,beta,mask,gap))
        if mode=='entropy':
            tau=max(-f['margin'],1e-3)/(max(logcount,1.) if temperature=='entropy_gap' else 1.)
            if f['type']=='complement':
                initial=log_partition(L,U,mask,np.asarray(P,np.int64),float(beta)/Q,tau,*arch.arrays,h)
                entropy_calls+=1
            else:initial=-f['margin']/tau
            entropy.append((tau,initial))
    records=[];calls=0
    for alt in alternatives:
        LL,UU=alt['L'],alt['U'];lower=upper=0.
        for f,P,R,beta,mask,gap in prepared:
            a,b=diff_lower(R,P,LL,UU);new=float((a+b*beta)/Q)
            gain=max(0.,new-f['margin']);den=max(-f['margin'],1e-9)
            upper+=min(1.,gain/den)
            lower+=min(1.,min(gain,gap)/den) if f['type']=='complement' else min(1.,gain/den)
        records.append({'lower_gain':lower,'upper_gain':upper,'evaluated':False})

    def evaluate(i):
        nonlocal calls
        alt=alternatives[i];gain=0.
        for f,P,R,beta,mask,gap in prepared:
            if f['type']=='complement':
                v,*_=complement(alt['L'],alt['U'],mask,np.asarray(P,np.int64),
                               beta.numerator,beta.denominator,*arch.arrays,h)
                calls+=1;new=int(v)/(beta.denominator*Q)
            else:
                a,b=diff_lower(R,P,alt['L'],alt['U']);new=float((a+b*beta)/Q)
            gain+=min(1.,max(0.,new-f['margin'])/max(-f['margin'],1e-9))
        records[i].update(evaluated=True,gain=gain,score=gain/alt['cost'])

    if mode=='entropy':
        for r,a in zip(records,alternatives):
            gain=0.
            for (f,P,R,beta,mask,gap),(tau,initial) in zip(prepared,entropy):
                if f['type']=='complement':
                    current=log_partition(a['L'],a['U'],mask,np.asarray(P,np.int64),float(beta)/Q,tau,*arch.arrays,h)
                    entropy_calls+=1
                else:
                    intercept,slope=diff_lower(R,P,a['L'],a['U'])
                    current=-float((intercept+slope*beta)/Q)/tau
                gain+=-np.expm1(min(0.,current-initial))
            r.update(gain=gain,score=gain/a['cost'])
    elif mode=='witness':
        for r,a in zip(records,alternatives):
            r.update(gain=r['upper_gain'],score=r['upper_gain']/a['cost'])
    elif mode=='full':
        for i in range(len(records)):evaluate(i)
    elif mode=='gap':
        first=max(range(len(records)),key=lambda i:(records[i]['lower_gain']/alternatives[i]['cost'],-alternatives[i]['cost']))
        evaluate(first)
        while True:
            done=[i for i,r in enumerate(records) if r['evaluated']]
            best=max(done,key=lambda i:(records[i]['score'],-alternatives[i]['cost']))
            # A small scheduling tolerance keeps roundoff away from pruning ties.
            remaining=[i for i,r in enumerate(records) if not r['evaluated'] and
                       r['upper_gain']/alternatives[i]['cost']>=records[best]['score']-1e-10]
            if not remaining:break
            evaluate(max(remaining,key=lambda i:records[i]['upper_gain']/alternatives[i]['cost']))
    else:raise ValueError(mode)
    eligible=[i for i,r in enumerate(records) if 'score' in r]
    best=max(eligible,key=lambda i:(records[i]['score'],-alternatives[i]['cost']))
    return best,records,{'pricing':mode,'two_rival_calls':gap_queries,
        'full_complement_calls':calls,'partition_calls':entropy_calls,
        'full_alternatives':sum(r['evaluated'] for r in records),
        'screened_alternatives':sum('score' not in r for r in records)}
