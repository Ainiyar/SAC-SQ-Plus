"""Exact proof-relevant automaton pruning for exponential acquisition.

The automaton follows archived prefixes until the first departure, then uses
the ordinary segment DAG. It accepts exactly paths outside the archive.
Delete an arc only when its exact shortest full accepted path has positive
exchange margin. Every nonpositive rival survives; surviving paths are fixed
throughout a materialization menu, preserving its submodular exponential gain.
"""
from fractions import Fraction as F
import time
import numpy as np
from numba import njit
from safe_archive import Q,INF,diff_lower

def topology(arch,n,h):
    key=(len(arch.A),n,h)
    cached=getattr(arch,'_proof_relevance_topology',None)
    if cached is not None and cached[0]==key:return cached[1]
    nodes=arch.nodes;nt=len(nodes)
    positions=np.r_[[x['v'] for x in nodes],np.arange(n+1)].astype(np.int32)
    ptr=[0];dest=[]
    for u,s in enumerate(positions):
        children=nodes[u]['children'] if u<nt else {}
        for t in range(int(s)+h,n+1):dest.append(children.get(t,nt+t))
        ptr.append(len(dest))
    graph=(np.asarray(ptr,np.int32),np.asarray(dest,np.int32),positions,
           np.argsort(positions,kind='stable').astype(np.int32),nt+n)
    arch._proof_relevance_topology=(key,graph)
    return graph

@njit(cache=True)
def exact_prune(L,U,mask,edges,bn,bd,ptr,dest,pos,order,accept):
    N=len(pos);fw=np.full(N,INF,np.int64);bw=fw.copy();fw[0]=0;bw[accept]=0
    for u in order:
        if fw[u]==INF:continue
        s=pos[u]
        for k in range(ptr[u],ptr[u+1]):
            v=dest[k];t=pos[v];w=(U[s,t] if mask[s,t] else L[s,t])*bd+bn
            fw[v]=min(fw[v],fw[u]+w)
    for j in range(len(order)-1,-1,-1):
        u=order[j];s=pos[u]
        for k in range(ptr[u],ptr[u+1]):
            v=dest[k]
            if bw[v]==INF:continue
            t=pos[v];w=(U[s,t] if mask[s,t] else L[s,t])*bd+bn
            bw[u]=min(bw[u],w+bw[v])
    base=np.int64(0)
    for s,t in edges:base+=U[s,t]*bd+bn
    active=np.zeros(len(dest),np.uint8);counts=np.zeros(N,np.int32);removed_gap=INF
    for u in range(N):
        if fw[u]==INF:continue
        s=pos[u]
        for k in range(ptr[u],ptr[u+1]):
            v=dest[k]
            if bw[v]==INF:continue
            t=pos[v];w=(U[s,t] if mask[s,t] else L[s,t])*bd+bn
            gap=fw[u]+w+bw[v]-base
            if gap<=0:active[k]=1;counts[u]+=1
            else:removed_gap=min(removed_gap,gap)
    newptr=np.zeros(N+1,np.int32)
    for u in range(N):newptr[u+1]=newptr[u]+counts[u]
    newdest=np.empty(newptr[-1],np.int32);offset=0
    for k in range(len(dest)):
        if active[k]:newdest[offset]=dest[k];offset+=1
    margin=INF if fw[accept]==INF else fw[accept]-base
    return newptr,newdest,margin,removed_gap

@njit(cache=True)
def minimum_margin(L,U,mask,edges,bn,bd,ptr,dest,pos,order,accept):
    distance=np.full(len(pos),INF,np.int64);distance[accept]=0
    for j in range(len(order)-1,-1,-1):
        u=order[j];s=pos[u]
        for k in range(ptr[u],ptr[u+1]):
            v=dest[k]
            if distance[v]==INF:continue
            t=pos[v];w=(U[s,t] if mask[s,t] else L[s,t])*bd+bn
            distance[u]=min(distance[u],w+distance[v])
    base=np.int64(0)
    for s,t in edges:base+=U[s,t]*bd+bn
    return INF if distance[0]==INF else distance[0]-base

@njit(cache=True)
def partition_core(L,U,mask,edges,beta,tau,ptr,dest,pos,order,accept):
    N=len(pos);distance=np.full(N,np.inf);distance[accept]=0.
    mass=np.zeros(N);mass[accept]=1.
    for j in range(len(order)-1,-1,-1):
        u=order[j];s=pos[u]
        for k in range(ptr[u],ptr[u+1]):
            v=dest[k];t=pos[v];w=(U[s,t] if mask[s,t] else L[s,t])/Q+beta
            distance[u]=min(distance[u],w+distance[v])
        if not np.isfinite(distance[u]):continue
        for k in range(ptr[u],ptr[u+1]):
            v=dest[k]
            if mass[v]==0.:continue
            t=pos[v];w=(U[s,t] if mask[s,t] else L[s,t])/Q+beta
            mass[u]+=np.exp(-max(0.,w+distance[v]-distance[u])/tau)*mass[v]
    base=0.
    for s,t in edges:base+=U[s,t]/Q+beta
    value=-np.inf if mass[0]==0. else (base-distance[0])/tau+np.log(mass[0])
    return value,distance,mass

@njit(cache=True)
def edge_marginals(L,U,mask,beta,tau,ptr,dest,pos,order,accept,distance,mass):
    probability=np.zeros(L.shape);arrivals=np.zeros(len(pos));arrivals[0]=1.
    for u in order:
        if arrivals[u]==0. or mass[u]==0.:continue
        s=pos[u]
        for k in range(ptr[u],ptr[u+1]):
            v=dest[k]
            if mass[v]==0.:continue
            t=pos[v];w=(U[s,t] if mask[s,t] else L[s,t])/Q+beta
            conditional=np.exp(-max(0.,w+distance[v]-distance[u])/tau)*mass[v]/mass[u]
            flow=arrivals[u]*conditional
            probability[s,t]+=flow;arrivals[v]+=flow
    return probability

@njit(cache=True)
def log_path_count(ptr,dest,pos,order,accept):
    value=np.full(len(pos),-np.inf);value[accept]=0.
    for j in range(len(order)-1,-1,-1):
        u=order[j]
        for k in range(ptr[u],ptr[u+1]):value[u]=np.logaddexp(value[u],value[dest[k]])
    return value[0]

def pruned_graph(arch,L,U,P,beta,h=8):
    ptr,dest,pos,order,accept=topology(arch,len(L)-1,h)
    mask=np.zeros(L.shape,np.uint8)
    for e in P:mask[e]=1
    edges=np.asarray(P,np.int64);beta=F(beta)
    pp,dd,gap,cutgap=exact_prune(L,U,mask,edges,beta.numerator,beta.denominator,ptr,dest,pos,order,accept)
    stats={'total_arcs':len(dest),'retained_arcs':len(dd),
           'exact_min_margin':None if gap==INF else int(gap)/(beta.denominator*Q),
           'minimum_removed_arc_margin':None if cutgap==INF else int(cutgap)/(beta.denominator*Q)}
    return (pp,dd,pos,order,accept),mask,stats

def prepare_context(arch,L,U,failures,temperature='debt',h=8,need_salience=True,prune=True,source_candidates=None):
    start=time.process_time();items=[];wl=np.zeros(L.shape);wu=np.zeros(L.shape)
    source_values=None
    if source_candidates is not None:
        source_candidates=list(map(tuple,source_candidates));source_index={e:i for i,e in enumerate(source_candidates)}
        ss=np.array([e[0] for e in source_candidates],np.int64);tt=np.array([e[1] for e in source_candidates],np.int64)
        widths=(U[ss,tt]-L[ss,tt])/Q;source_values=np.zeros(len(ss))
    original_arcs=0;retained_arcs=0;comp=0;input_arcs=0;cache_hits=0
    if prune=='persistent':
        from persistent_relevance import begin_epoch
        begin_epoch(arch,failures,h)
    for failure in failures:
        weight=failure.get('weight',1.)
        P=tuple(map(tuple,failure['candidate']));R=tuple(map(tuple,failure['rival']))
        beta=F(failure['beta']);tau=max(-failure['margin'],.001)
        item={'failure':failure,'P':P,'R':R,'beta':beta,'tau':tau,'graph':None,'weight':weight}
        if failure['type']=='complement':
            if prune:
                if prune=='persistent':
                    from persistent_relevance import pruned_graph_cached
                    graph,mask,stats=pruned_graph_cached(arch,L,U,P,beta,h)
                else:graph,mask,stats=pruned_graph(arch,L,U,P,beta,h)
                if stats['exact_min_margin'] is None or stats['exact_min_margin']>0:
                    raise ValueError('failed obligation has no nonpositive external rival')
            else:
                graph=topology(arch,len(L)-1,h);mask=np.zeros(L.shape,np.uint8)
                for edge in P:mask[edge]=1
                stats={'total_arcs':len(graph[1]),'retained_arcs':len(graph[1])}
            if temperature=='entropy_gap':tau/=max(1.,log_path_count(*graph))
            edges=np.asarray(P,np.int64)
            initial,distance,mass=partition_core(L,U,mask,edges,float(beta)/Q,tau,*graph)
            if need_salience or source_values is not None:
                prob=edge_marginals(L,U,mask,float(beta)/Q,tau,*graph,distance,mass)
                prob=np.clip(prob,0.,1.)
                if need_salience:wl+=weight*(1-mask)*prob/tau;wu+=weight*mask*(1-prob)/tau
                if source_values is not None:
                    q=prob[ss,tt];pflag=mask[ss,tt]
                    source_values+=weight*((1-pflag)*q+pflag*(1-q))*(-np.expm1(-widths/tau))
            item.update(graph=graph,mask=mask,edges=edges,initial=initial,tau=tau,pruning=stats)
            original_arcs+=stats['total_arcs'];retained_arcs+=stats['retained_arcs'];comp+=1
            input_arcs+=stats.get('pruning_input_arcs',stats['total_arcs']);cache_hits+=stats.get('relevance_cache_hit',False)
        else:
            item['initial']=-failure['margin']/tau
            if need_salience:
                for edge in set(R)-set(P):wl[edge]+=weight/tau
                for edge in set(P)-set(R):wu[edge]+=weight/tau
            if source_values is not None:
                for edge in set(R)^set(P):
                    if edge in source_index:
                        j=source_index[edge];source_values[j]+=weight*(-np.expm1(-widths[j]/tau))
        items.append(item)
    return {'items':items,'weights':(wl,wu),'source_values':source_values,'source_candidates':source_candidates,
        'mode':'persistent' if prune=='persistent' else 'relevant' if prune else 'automaton','stats':{
        'relevance_graphs':comp,'relevance_original_arcs':original_arcs,
        'relevance_retained_arcs':retained_arcs,
        'relevance_pruning_input_arcs':input_arcs,'relevance_cache_hits':cache_hits,
        'relevance_arc_fraction':None if not original_arcs else retained_arcs/original_arcs,
        'proof_relevance_preparation_cpu':time.process_time()-start}}

def price_context(context,L,U,alternatives):
    records=[];calls=context['stats']['relevance_graphs']
    for alt in alternatives:
        gain=0.
        for item in context['items']:
            P,R,beta,tau=item['P'],item['R'],item['beta'],item['tau']
            if item['graph'] is not None:
                current,*_=partition_core(alt['L'],alt['U'],item['mask'],item['edges'],
                                           float(beta)/Q,tau,*item['graph']);calls+=1
            else:
                a,b=diff_lower(R,P,alt['L'],alt['U']);current=-float((a+b*beta)/Q)/tau
            gain+=item['weight']*(-np.expm1(min(0.,current-item['initial'])))
        records.append({'gain':float(gain),'score':float(gain/alt['cost']),'evaluated':False})
    best=max(range(len(records)),key=lambda i:(records[i]['score'],-alternatives[i]['cost']))
    return best,records,{'pricing':context['mode'],'partition_calls':calls,
        'full_complement_calls':0,'two_rival_calls':0,'full_alternatives':0,'screened_alternatives':0,
        **context['stats']}

def price_hard_context(context,L,U,alternatives):
    records=[];calls=0
    for alt in alternatives:
        gain=0.
        for item in context['items']:
            f=item['failure'];P,R,beta=item['P'],item['R'],item['beta']
            if item['graph'] is not None:
                v=minimum_margin(alt['L'],alt['U'],item['mask'],item['edges'],beta.numerator,beta.denominator,*item['graph'])
                new=float('inf') if v==INF else int(v)/(beta.denominator*Q);calls+=1
            else:
                a,b=diff_lower(R,P,alt['L'],alt['U']);new=float((a+b*beta)/Q)
            # Every margin at this beta is an integer multiple of this quantum.
            # Deleted paths already exceed zero, hence are at least quantum.
            # Capping here is exactly the full-graph strict-proof utility,
            # including zero/tiny initial deficits.
            target=1/(beta.denominator*Q)
            gain+=item['weight']*min(1.,max(0.,new-f['margin'])/(target-f['margin']))
        records.append({'gain':gain,'score':gain/alt['cost'],'evaluated':True})
    best=max(range(len(records)),key=lambda i:(records[i]['score'],-alternatives[i]['cost']))
    return best,records,{'pricing':'relevant_hard','relevant_minimum_calls':calls,
        'partition_calls':context['stats']['relevance_graphs'],'full_complement_calls':0,
        'two_rival_calls':0,'full_alternatives':len(records),'screened_alternatives':0,**context['stats']}
