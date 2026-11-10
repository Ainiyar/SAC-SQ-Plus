"""Whole-complement marginal prediction for strict evidence materialization.

Predictions choose work only. Actual updates always use directed Decimal.
The complete-complement gap or relevant-path partition credits replacement
rivals when estimating an action. Union geometry prices and constructs only
the source rows covered by the requested target intervals.
"""
import time
import numpy as np
from scipy.special import expit
from safe_archive import Q,ia,fit,complement
from lower_closure import close_lower

def forecasts(L,U,X,y,theta,e,failures,row_cost,query_cost,bank_allowed,pricing,edge_weights=None,closure_kernel='full',bank_cached=False,closure_piece_h=8,bank_geometry='full'):
    n=len(y)
    ss,tt=np.triu_indices(n+1,8)
    z=X@theta;p=expit(z)
    prefix=np.r_[0,np.cumsum(np.logaddexp(0,z)-y*z)]
    gp=np.vstack([np.zeros((1,X.shape[1])),np.cumsum(X*(p-y)[:,None],axis=0)])
    m=tt-ss
    f=prefix[tt]-prefix[ss]+.05*m*float(theta@theta)
    g=gp[tt]-gp[ss]+.1*m[:,None]*theta
    low=np.maximum(0.,f-np.sum(g*g,axis=1)/(.2*m))
    pl=np.maximum(L[ss,tt],np.minimum(U[ss,tt],np.rint(low*Q).astype(np.int64)))
    pu=np.minimum(U[ss,tt],np.maximum(pl,np.rint(f*Q).astype(np.int64)))
    dl=(pl-L[ss,tt])/Q;du=(U[ss,tt]-pu)/Q
    weights={}
    for f0 in failures:
        R,P=set(f0['rival']),set(f0['candidate']);den=max(-f0['margin'],1e-9)
        for edge in R-P:weights[edge]=weights.get(edge,(0.,0.))[0]+1/den,weights.get(edge,(0.,0.))[1]
        for edge in P-R:weights[edge]=weights.get(edge,(0.,0.))[0],weights.get(edge,(0.,0.))[1]+1/den
    if edge_weights is None:
        importance=np.array([weights.get((int(s),int(t)),(0.,0.))[0]*a+
                             weights.get((int(s),int(t)),(0.,0.))[1]*b
                             for s,t,a,b in zip(ss,tt,dl,du)])
    else:importance=edge_weights[0][ss,tt]*dl+edge_weights[1][ss,tt]*du
    ranks=np.lexsort((dl+du,importance))[::-1]
    informative=ranks[(dl[ranks]+du[ranks])>1e-10]
    ei=int(np.flatnonzero((ss==e[0])&(tt==e[1]))[0])
    candidates=[('direct',np.array([ei],np.int64),(0 if bank_cached else (e[1]-e[0])*row_cost)+query_cost)]
    if bank_allowed:
        seen=set()
        for size in [64,256,1024,len(informative)]:
            indices=np.unique(np.r_[ei,informative[:size]])
            key=tuple(indices)
            if key in seen:continue
            seen.add(key)
            if bank_geometry=='union':
                from interval_bank import interval_union
                rows=sum(t-s for s,t in interval_union(zip(ss[indices],tt[indices])))
            else:rows=n
            candidates.append(('bank',indices,(0 if bank_cached else rows*row_cost)+len(indices)*query_cost))
    alternatives=[];closure_base=L
    for kind,indices,cost in candidates:
        closure_start=time.process_time();UU=U.copy()
        UU[ss[indices],tt[indices]]=pu[indices]
        # Each standalone closed forecast is a fixed action outcome. Their
        # edgewise union retains diminishing returns; cross-action closure may
        # supply additional progress beyond that lower-bound surrogate.
        if closure_kernel=='incremental' and pricing!='entropy_raw':
            from incremental_lower_closure import close_updates
            LL,closure_stats=close_updates(closure_base,ss[indices],tt[indices],pl[indices],closure_piece_h)
            closure_base=LL
            closure_stats['nested_menu_reuse']=True
        else:
            LL=L.copy();LL[ss[indices],tt[indices]]=pl[indices]
            if pricing!='entropy_raw':close_lower(LL,closure_piece_h)
            closure_stats={'method':'raw' if pricing=='entropy_raw' else 'full'}
        alternatives.append({'kind':kind,'indices':indices,'cost':cost,'L':LL,'U':UU,
                             'closure':closure_stats,'bound_closure_cpu':time.process_time()-closure_start})
    return ss,tt,alternatives


def acquire(arch,L,U,X,y,Xd,yd,e,failures,row_cost,query_cost,bank_allowed,
            pricing='full',temperature='entropy_gap',source_edges=None,theta_cache=None,
            target_priority='witness',closure_kernel='full',source_targets=None,bank_cache=None,
            source_priority='witness',source_candidates=None,charge_round_overhead=False,round_prework_cpu=0.,closure_piece_h=8,bank_geometry='full'):
    start=time.process_time();fitcpu=0.;alternatives=[];sources=[]
    edge_weights=None;salience_stats={'target_priority':'witness','marginal_distribution_calls':0}
    salience_start=time.process_time()
    relevance_context=None
    if pricing in ['relevant','automaton','persistent','relevant_hard','projected']:
        from proof_relevance import prepare_context
        relevance_context=prepare_context(arch,L,U,failures,temperature,8,target_priority=='entropy','persistent' if pricing in ['persistent','relevant_hard','projected'] else pricing=='relevant',
                                          source_candidates if source_priority=='relevant_gain' else None)
        if target_priority=='entropy':
            edge_weights=relevance_context['weights']
            salience_stats={'target_priority':pricing,'marginal_distribution_calls':relevance_context['stats']['relevance_graphs']}
    elif target_priority=='entropy':
        from entropy_salience import salience
        edge_weights,salience_stats=salience(arch,L,U,failures,temperature)
    salience_cpu=time.process_time()-salience_start
    source_selection={'source_priority':source_priority}
    if source_priority=='relevant_gain':
        values=relevance_context['source_values'];edges=relevance_context['source_candidates']
        j=max(range(len(edges)),key=lambda j:(values[j]/((edges[j][1]-edges[j][0])*row_cost+query_cost),-(edges[j][1]-edges[j][0])))
        source_edges=[edges[j]]
        source_selection.update(source_priority_value=float(values[j]),source_priority_candidates=len(edges))
    for source in (source_edges if source_edges is not None else [e]):
        source=tuple(source);cached=theta_cache is not None and source in theta_cache
        if cached:theta,conv=theta_cache[source];fcpu=0.
        else:
            t0=time.process_time();theta,conv=fit(X,y,*source,tol=1e-7)
            fcpu=time.process_time()-t0
            if theta_cache is not None:theta_cache[source]=(theta,conv)
        fitcpu+=fcpu
        target=source if source_targets is None else source_targets.get(source,source)
        ss,tt,options=forecasts(L,U,X,y,theta,target,failures,row_cost,query_cost,bank_allowed,pricing,edge_weights,closure_kernel,
                               bank_cache is not None and source in bank_cache,closure_piece_h,bank_geometry)
        for option in options:option['source_index']=len(sources)
        alternatives.extend(options)
        sources.append((source,theta,conv,fcpu,cached))
    from action_pricing import price
    pricing_start=time.process_time()
    chosen,evaluations,pricing_stats=price(arch,L,U,failures,alternatives,
                                         'entropy' if pricing=='entropy_raw' else pricing,
                                         temperature=temperature,relevance_context=relevance_context)
    pricing_cpu=time.process_time()-pricing_start
    charged_overhead=0.
    if charge_round_overhead:
        charged_overhead=round_prework_cpu+time.process_time()-start
        for a,r in zip(alternatives,evaluations):
            r['materialization_only_score']=r['score']
            r['score']=r['gain']/(a['cost']+charged_overhead)
        chosen=max(range(len(alternatives)),key=lambda i:(evaluations[i]['score'],-alternatives[i]['cost']))
    selected=alternatives[chosen]
    e,theta,conv,_,_=sources[selected['source_index']]
    target=e if source_targets is None else source_targets.get(e,e)
    kind,indices,cost=selected['kind'],selected['indices'],selected['cost']
    value=evaluations[chosen]['gain']
    guidecpu=time.process_time()-start-fitcpu;t0=time.process_time()
    bank_hit=bank_cache is not None and e in bank_cache
    if bank_hit:bank=bank_cache[e];offset=0;rows_built=0
    elif kind=='bank' and bank_geometry=='union':
        from interval_bank import UnionBank
        bank=UnionBank(Xd,yd,theta,zip(ss[indices],tt[indices]));offset=0;rows_built=bank.rows_built
    elif kind=='bank':
        bank=ia.PrefixBank.build(Xd,yd,theta);offset=0;rows_built=len(y)
        if bank_cache is not None:bank_cache[e]=bank
    else:
        bank=ia.PrefixBank.build(Xd[target[0]:target[1]],yd[target[0]:target[1]],theta);offset=target[0]
        rows_built=target[1]-target[0]
    buildcpu=time.process_time()-t0;t0=time.process_time();improved=[];upper=[]
    for idx in indices:
        s,t=int(ss[idx]),int(tt[idx]);low,high=bank.query(s-offset,t-offset,.1,theta)
        l=int(ia.mul_dn(ia.Decimal(Q),low).to_integral_value(rounding=ia.ROUND_FLOOR))
        u=int(ia.mul_up(ia.Decimal(Q),high).to_integral_value(rounding=ia.ROUND_CEILING))
        ol,ou=int(L[s,t]),int(U[s,t]);nl,nu=max(ol,l),min(ou,u)
        if nl>nu:raise ValueError('invalid directed bounds')
        L[s,t],U[s,t]=nl,nu
        a,b=(int(L[s,t])-ol)/Q,(ou-int(U[s,t]))/Q
        if b>0:upper.append((s,t))
        if a>0 or b>0:improved.append({'edge':[s,t],'lower_gain':a,'upper_gain':b})
    return {'source_edge':e,'kind':kind,'strict_queries':len(indices),
        'direct_target_edge':target,'bank_cache_hit':bank_hit,'strict_rows_built':rows_built,
        'bank_geometry':bank_geometry,
        'bank_components':getattr(bank,'components',None),
        'fit_cpu':fitcpu,'guide_cpu':guidecpu,'bank_build_cpu':buildcpu,
        'strict_query_cpu':time.process_time()-t0,'predicted_full_proof_gain':value,
        'predicted_materialization_cpu':cost,'improved':improved,'upper_changed':upper,
        'charged_round_overhead_cpu':charged_overhead,
        'predicted_round_cpu':cost+charged_overhead,
        'theta':list(map(float,theta)),
        'advisory_lower_closure':pricing!='entropy_raw',
        'entropy_temperature_policy':temperature,
        'target_salience_cpu':salience_cpu,'pricing_cpu':pricing_cpu,
        'forecast_bound_closure_cpu':sum(a['bound_closure_cpu'] for a in alternatives),
        'closure_kernel':closure_kernel,
        'numerical_source_fits':sum(not s[4] for s in sources),
        'source_menu':[{'edge':s[0],'theta':list(map(float,s[1])),
                        'fit_cpu':s[3],'cached':s[4]} for s in sources],
        'alternatives':[{'kind':a['kind'],'source_index':a['source_index'],'queries':len(a['indices']),
                         'lower_closure':a['closure'],
                         'predicted_cpu':a['cost'],**r} for a,r in zip(alternatives,evaluations)],
        **pricing_stats,**salience_stats,**source_selection,**conv}
