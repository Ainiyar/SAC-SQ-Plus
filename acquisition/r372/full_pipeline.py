"""Full-cost prospective Occupancy pilot and official R CROPS bridge."""
import argparse,gzip,json,resource,socket,subprocess,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from safe_archive import (SafeArchive,Q,ia,fit,direct as strict_direct,phase_audit,complement)
from incremental_upper import initialize,decrease
from lower_closure import close_lower

def direct(*args):
    return strict_direct(*args,fit_tol=1e-7)

def init_bounds(n,h):
    L=np.zeros((n+1,n+1),np.int64)
    U=np.zeros_like(L)
    _,log2=ia.ln_i_positive((ia.TWO,ia.TWO))
    for s in range(n-h+1):
        for t in range(s+h,n+1):
            u=ia.mul_up(ia.Decimal(t-s),log2)
            U[s,t]=int(ia.mul_up(ia.Decimal(Q),u).to_integral_value(rounding=ia.ROUND_CEILING))
    return L,U

def preprocess(start,dataset='occupancy',n=128):
    # Frozen-input paired R337 comparison: load the exact X/y saved by R333.
    if dataset=='occupancy' and start==768 and n==512:
        z=np.load('data/processed/frozen_reference/occupancy_768_n512_Xy.npz')
        return z['X'],z['y']
    # Recreate preprocessing from raw source so its cost is in each arm.
    import io,zipfile,pandas as pd
    if dataset=='secom':
        from secom_data import preprocess as secom_preprocess
        return secom_preprocess(start,n)
    if dataset=='eeg':
        from scipy.io import arff
        with zipfile.ZipFile('eeg/uci_eeg.zip') as z:
            member=next(x for x in z.namelist() if x.lower().endswith('.arff'))
            data,meta=arff.loadarff(io.StringIO(z.read(member).decode('utf8')))
        a=np.column_stack([data[x][start:start+n] for x in data.dtype.names[:-1]]).astype(float)
        yy=data[data.dtype.names[-1]][start:start+n].astype(np.int64)
        sd=a.std(axis=0);sd[sd<1e-12]=1.
        X=np.column_stack([np.ones(len(a)),(a-a.mean(axis=0))/sd,np.linspace(-1,1,len(a))])
        return X,yy
    with zipfile.ZipFile('occupancy/uci_occupancy.zip') as z:
        member=next(x for x in z.namelist() if Path(x).name=='datatraining.txt')
        df=pd.read_csv(io.BytesIO(z.read(member)))
    df=df.sort_values('date',kind='stable').reset_index(drop=True)
    block=df.iloc[start:start+n]
    a=block[['Temperature','Humidity','Light','CO2','HumidityRatio']].to_numpy(float)
    sd=a.std(axis=0);sd[sd<1e-12]=1.
    X=np.column_stack([np.ones(len(a)),(a-a.mean(axis=0))/sd,np.linspace(-1,1,len(a))])
    return X,block.Occupancy.to_numpy(np.int64)

def shared_start(L,U,X,y,Xd,yd):
    n=len(y);records=[]
    for e in [(0,n),(0,n//2),(n//2,n),None]:
        start=time.process_time()
        theta=np.zeros(X.shape[1]) if e is None else fit(X,y,*e)[0]
        fit_cpu=time.process_time()-start
        start=time.process_time();bank=ia.PrefixBank.build(Xd,yd,theta)
        build_cpu=time.process_time()-start;start=time.process_time();queries=0
        for s in range(n-8+1):
            for t in range(s+8,n+1):
                low,high=bank.query(s,t,.1,theta)
                l=int(ia.mul_dn(ia.Decimal(Q),low).to_integral_value(rounding=ia.ROUND_FLOOR))
                u=int(ia.mul_up(ia.Decimal(Q),high).to_integral_value(rounding=ia.ROUND_CEILING))
                L[s,t]=max(int(L[s,t]),l);U[s,t]=min(int(U[s,t]),u);queries+=1
        records.append({'source':e,'fit_cpu':fit_cpu,'bank_cpu':build_cpu,
                        'strict_query_cpu':time.process_time()-start,'queries':queries,
                        'theta':list(map(float,theta))})
    return records

class StrictPELT:
    """Independent delayed-pruning recurrence; exact interval pruning tests."""
    def __init__(self,L,U,X,y,Xd,yd):
        self.L,self.U,self.X,self.y,self.Xd,self.yd=L,U,X,y,Xd,yd
        self.cache={};self.calls=[];self.n=len(y)
    def query(self,s,t):
        if (s,t) not in self.cache:
            self.cache[s,t]=direct(self.L,self.U,self.X,self.y,self.Xd,self.yd,(s,t))
            if self.n>=1024 and len(self.cache)%10000==0:
                print('STRICT_PROGRESS','fits',len(self.cache),'process_cpu',time.process_time(),flush=True)
    def solve(self,beta):
        n=self.n;h=8;b=F.from_float(beta)*Q
        DL=[None]*(n+1);DU=[None]*(n+1);pred=[-1]*(n+1)
        DL[0]=DU[0]=F(0)
        expiry={0:n+1};active=[0];prunes=0;evals=0;start=time.process_time()
        for t in range(h,n+1):
            active=[s for s in active if expiry[s]>t]
            eligible=[s for s in active if s<=t-h]
            for s in eligible:self.query(s,t)
            valsL=[DL[s]+int(self.L[s,t])+b for s in eligible]
            valsU=[DU[s]+int(self.U[s,t])+b for s in eligible]
            j=min(range(len(eligible)),key=lambda i:valsU[i])
            DL[t]=min(valsL);DU[t]=valsU[j];pred[t]=eligible[j];evals+=len(eligible)
            for s in eligible:
                if DL[s]+int(self.L[s,t])>DU[t]:
                    expiry[s]=min(expiry[s],t+h);prunes+=1
            if t<n:active.append(t);expiry[t]=n+1
        P=[];t=n
        while t:
            s=pred[t];P.append((s,t));t=s
        P=P[::-1]
        unpen=sum(int(self.U[e]) for e in P)/Q
        rec={'beta':beta,'path':P,'ncp':len(P)-1,'unpenalized':unpen,
             'penalized_interval':[str(DL[n]),str(DU[n])],
             'prune_events':prunes,'edge_evals':evals,'cpu':time.process_time()-start}
        self.calls.append(rec)
        print('PELT',beta,'counts',len(P),'fits',len(self.cache),flush=True)
        return rec

def official_crops(oracle,outbase):
    before=resource.getrusage(resource.RUSAGE_CHILDREN)
    with socket.socket() as server:
        server.bind(('127.0.0.1',0));server.listen(1);server.settimeout(60)
        port=server.getsockname()[1]
        with open(str(outbase)+'_R.log','w') as log:
            proc=subprocess.Popen(['Rscript','crops_bridge.R',str(port),'.1','1.2',
                                   str(outbase)+'_R.json'],stdout=log,stderr=log)
            conn,_=server.accept()
            with conn,conn.makefile('rwb') as file:
                while True:
                    line=file.readline().decode().strip()
                    if line=='DONE':break
                    if not line:raise RuntimeError('CROPS connection closed early')
                    rec=oracle.solve(float(line))
                    file.write((json.dumps({'unpenalized':rec['unpenalized'],
                                             'ncp':rec['ncp']})+'\n').encode());file.flush()
            code=proc.wait(timeout=60)
            if code:raise RuntimeError('official R CROPS failed: '+str(code))
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    return after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime

def finish(L,U,X,y,Xd,yd,budget,adaptive=False,costs=None,closure=False,lookahead=False,bank_cap=12,pricing='full',temperature='entropy_gap',source_pool=1,theta_cache=None,target_priority='witness',resume=None,closure_kernel='full',acquirer=None,source_priority='witness',obligation_weight='uniform',charge_round_overhead=False,closure_piece_h=8,bank_geometry='full',acquisition_objective='all_rivals',joint_bank_pricing=False,certificate_bank_pricing=False):
    arch=SafeArchive(len(y),8,Q//10,6*Q//5,'archive',True)
    hist=[];snapshots=[];selected=set();state=None;banks=0;stop_reason=None
    if resume is not None:
        paths=[tuple(map(tuple,p)) for p in resume['archive_paths']]
        arch.A={p:i for i,p in enumerate(paths)};arch.cover={p:None for p in paths}
        for c in resume.get('proven_cover_intervals',[]):
            arch.cover[paths[c['path_id']]]=(F(c['lo_beta'])*Q,F(c['hi_beta'])*Q)
        for c in resume.get('retirement_certificates',[]):arch.retired[paths[c['id']]]=c
        selected={tuple(r.get('direct_target_edge',r.get('source_edge',r.get('edge'))))
                  for r in resume['history'] if r.get('target_complete',True)}
        banks=sum(r.get('kind')=='bank' for r in resume['history'])
    lower_closed=False
    for step in range(budget+1):
        round_start=time.process_time()
        if closure and (closure_kernel=='full' or not lower_closed):close_lower(L,closure_piece_h)
        lower_closed=True
        C,new,migrations=arch.refresh(L,U);state=arch.check(L,U,C)
        if obligation_weight=='unique_beta':
            from collections import Counter
            multiplicity=Counter(F(f['beta']) for f in state['failures'] if f['type']=='complement')
            for f in state['failures']:
                f['weight']=1/multiplicity[F(f['beta'])] if f['type']=='complement' else 1.
        phase=phase_audit(L,U,C,Q//10,6*Q//5)
        snap={k:v for k,v in state.items() if k!='failures'}
        snap.update(step=step,counts=[c['k'] for c in C],failed=len(state['failures']),
                    phase_pass=phase['pass'],new=new,migrations=migrations)
        snapshots.append(snap)
        if not state['failures'] and phase['pass']:break
        if step==budget:stop_reason='ACTION_BUDGET';break
        # Strong single-obligation Partition/Gibbs baseline: retain the
        # identical evidence oracle, bank menu, closure and full certificate.
        # Only the proof obligation guiding the next acquisition is restricted.
        active_failures=state['failures']
        if acquisition_objective=='single_gibbs' and active_failures:
            eligible=[f for f in active_failures if any(
                tuple(a['edge']) not in selected for a in f['actions'])]
            if eligible:
                # Pick the hardest currently actionable certificate obligation.
                # Ties respect the verifier's deterministic obligation order.
                active_failures=[max(eligible,key=lambda f:-float(f['margin']))]
            else:
                active_failures=[]
        scores={}
        for f in active_failures:
            for a in f['actions']:
                e=tuple(a['edge'])
                if e not in selected:
                    scores[e]=scores.get(e,0)+f.get('weight',1.)*min(1.,a['width']/max(-f['margin'],1e-9))
        if not state['failures']:
            # Topology debt: refine only edges of the retained path set.
            for c in C:
                for e in c['path']:
                    width=int(U[e])-int(L[e])
                    if width>100 and e not in selected:scores[e]=width/Q
        if not scores:stop_reason='PRECISION_OR_NO_INFORMATIVE_TARGET';break
        e=max(scores,key=lambda e:(scores[e],1/(1+e[1]-e[0])))
        prior_lower=L.copy() if closure and closure_kernel=='incremental' else None
        try:
            if adaptive and state['failures']:
                if acquirer is not None:
                    rec=acquirer(arch,L,U,e,state['failures'])
                elif lookahead:
                    from proof_lookahead import acquire
                    source_edges=sorted(scores,key=lambda e:(scores[e],1/(1+e[1]-e[0])),reverse=True)[:source_pool]
                    rec=acquire(arch,L,U,X,y,Xd,yd,e,active_failures,costs[0],costs[1],bank_cap==0 or banks<bank_cap,pricing,temperature,
                                source_edges,theta_cache,target_priority,closure_kernel,
                                source_priority=source_priority,source_candidates=list(scores),
                                charge_round_overhead=charge_round_overhead,round_prework_cpu=time.process_time()-round_start,
                                closure_piece_h=closure_piece_h,bank_geometry=bank_geometry,joint_bank=joint_bank_pricing,certificate_bank=certificate_bank_pricing)
                else:
                    from adaptive_materialization import acquire
                    rec=acquire(L,U,X,y,Xd,yd,e,state['failures'],costs[0],costs[1],bank_cap==0 or banks<bank_cap)
                hist.append(rec);banks+=int(rec['kind']=='bank')
                if rec.get('target_complete',True):
                    selected.add(tuple(rec.get('direct_target_edge',rec['source_edge'])))
                for changed in rec['upper_changed']:arch.upper.update(*changed)
            else:
                rec=direct(L,U,X,y,Xd,yd,e);hist.append(rec)
                selected.add(e)
                if rec['upper_gain']>0:arch.upper.update(*e)
            lower_closed=False
            if prior_lower is not None:
                from incremental_lower_closure import close_updates
                closure_start=time.process_time();ss,tt=np.where(L>prior_lower)
                closed,stats=close_updates(prior_lower,ss,tt,L[ss,tt],closure_piece_h)
                L[:]=closed;lower_closed=True
                rec['published_lower_closure']=stats
                rec['published_lower_closure_cpu']=time.process_time()-closure_start
        except (ArithmeticError,ValueError,np.linalg.LinAlgError) as exc:
            stop_reason={'kind':'SOURCE_OR_NUMERICAL_FAILURE','edge':e,
                         'error_type':type(exc).__name__,'detail':str(exc)}
            break
        if step%10==0:print('REFINE',step,'debt',snap['failed'],'paths',len(arch.A),flush=True)
    status='COMPLETE_PHASE_PASS' if stop_reason is None and not state['failures'] and phase['pass'] else 'UNRESOLVED'
    from path_support import all_supports
    supports=all_supports(arch.A,L,U,Q//10,6*Q//5)
    external_complete=not any(f['type']=='complement' for f in state['failures'])
    return {'status':status,'actions':len(hist),'history':hist,'snapshots':snapshots,
            'stop_reason':stop_reason,
            'phase':phase,'archive_paths':[list(p) for p in arch.A],
            'archive_supports':supports,'support_penalty_integer_scale':Q,
            'complete_box_support_certificate':external_complete,
            'proven_cover_intervals':[{'path_id':arch.A[p],
                'lo_beta':str(interval[0]/Q),'hi_beta':str(interval[1]/Q)}
                for p,interval in arch.cover.items() if interval is not None],
            'retirement_certificates':list(arch.retired.values()),
            'last_failed':state['failures'],'upper_updates':{'updates':arch.upper.updates,
                'changed_states':arch.upper.states,'relaxations':arch.upper.relaxations}}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--start',type=int,required=True)
    ap.add_argument('--method',choices=['shared','adaptive','hybrid','lookahead','reliever','crops_pelt'],required=True)
    ap.add_argument('--budget',type=int,default=200);ap.add_argument('--out',required=True)
    ap.add_argument('--closure',action='store_true')
    ap.add_argument('--dataset',choices=['occupancy','eeg','secom'],default='occupancy')
    ap.add_argument('--bank-cap',type=int,default=12)
    ap.add_argument('--n',type=int,default=128)
    ap.add_argument('--relief-r',type=float,default=.5)
    ap.add_argument('--initial-policy',choices=['full4','sparse4','proposal_only'],default='full4')
    ap.add_argument('--bank-geometry',choices=['full','union'],default='full')
    ap.add_argument('--joint-bank-pricing',action='store_true',help='R371 exact joint full-beta menu pricing')
    ap.add_argument('--certificate-bank-pricing',action='store_true',help='R372 exact projected archive debt pricing')
    ap.add_argument('--pricing',choices=['full','gap','witness','entropy','entropy_raw','relevant','automaton','persistent','relevant_hard','projected'],default='full')
    ap.add_argument('--source-priority',choices=['witness','relevant_gain'],default='witness')
    ap.add_argument('--acquisition-objective',choices=['all_rivals','single_gibbs'],default='all_rivals')
    ap.add_argument('--obligation-weight',choices=['uniform','unique_beta'],default='uniform')
    ap.add_argument('--charge-round-overhead',action='store_true')
    ap.add_argument('--closure-geometry',choices=['tiling','packing'],default='tiling')
    ap.add_argument('--entropy-temperature',choices=['entropy_gap','debt'],default='entropy_gap')
    ap.add_argument('--source-pool',type=int,default=1)
    ap.add_argument('--target-priority',choices=['witness','entropy'],default='witness')
    ap.add_argument('--closure-kernel',choices=['full','incremental'],default='full')
    ap.add_argument('--cost-model-from')
    ap.add_argument('--source-policy',choices=['target','lazy_relief'],default='target')
    ap.add_argument('--reuse-proposal-fits',action='store_true')
    args=ap.parse_args()
    if args.joint_bank_pricing and args.certificate_bank_pricing:
        ap.error('choose only one optional bank pricer')
    closure_piece_h=1 if args.closure_geometry=='packing' else 8
    if args.bank_geometry=='union' and (args.source_policy!='target' or args.method!='lookahead'):
        ap.error('union evidence uses the lookahead target-source controller')
    if args.closure_geometry=='packing' and args.source_policy!='target':
        ap.error('packing evidence geometry currently uses the target-source controller')
    if args.charge_round_overhead and (args.source_pool!=1 or args.source_policy!='target' or args.pricing not in ['relevant','persistent','projected']):
        ap.error('round overhead pricing uses one target source and relevant/persistent pricing')
    if args.obligation_weight=='unique_beta' and args.pricing not in ['relevant','persistent','relevant_hard']:
        ap.error('unique_beta weights use relevant/persistent/relevant_hard pricing')
    if args.source_priority=='relevant_gain' and (args.source_policy!='target' or args.source_pool!=1 or args.pricing not in ['relevant','persistent','relevant_hard']):
        ap.error('relevant_gain uses one target source with relevant/persistent/relevant_hard pricing')
    # Kernels compile before both runtime boundaries. pandas import also warmed.
    X,y=preprocess(args.start,args.dataset,args.n);L,U=init_bounds(len(y),8)
    warm=SafeArchive(len(y),8,Q//10,6*Q//5,'archive',True)
    C,_,_=warm.refresh(L,U);warm_state=warm.check(L,U,C)
    decrease(warm.upper.D.copy(),warm.upper.P.copy(),U,0,8,8)
    close_lower(L.copy(),8)
    if args.closure_kernel=='incremental':
        from incremental_lower_closure import close_updates
        close_updates(L,[0],[8],[int(L[0,8])+1],closure_piece_h,force_sparse=True)
    if args.method in ['lookahead','reliever']:
        from proof_lookahead import forecasts
        from action_pricing import price
        price(warm,L,U,warm_state['failures'],[{'L':L,'U':U,'cost':1.}],
              'entropy' if args.pricing=='entropy_raw' else args.pricing,
              temperature=args.entropy_temperature)
        if args.target_priority=='entropy':
            if args.pricing in ['relevant','automaton','persistent','relevant_hard','projected']:
                from proof_relevance import prepare_context
                prepare_context(warm,L,U,warm_state['failures'],args.entropy_temperature,need_salience=True,prune='persistent' if args.pricing in ['persistent','relevant_hard','projected'] else args.pricing=='relevant')
            else:
                from entropy_salience import salience
                salience(warm,L,U,warm_state['failures'],args.entropy_temperature)
    # Shared pre-timing warmup for both arms; selection remains conditional in finish().
    from certificate_bank_pricing import price_by_debt
    price_by_debt(warm,warm_state['failures'],[{'L':L,'U':U,'cost':1.}])
    start=time.process_time();wall=time.perf_counter()
    X,y=preprocess(args.start,args.dataset,args.n);L,U=init_bounds(len(y),8)
    Xd=[[ia.dfloat(v) for v in row] for row in X];yd=list(map(int,y))
    prep=time.process_time()-start
    childcpu=0.;initial={};t=time.process_time()
    costs=None;theta_cache={} if args.source_pool>1 or args.reuse_proposal_fits or args.source_policy=='lazy_relief' else None
    if args.method in ['shared','adaptive','hybrid','lookahead','reliever']:
        if args.method=='reliever':
            from relief_initial import initialize as relief_start
            initial=relief_start(L,U,X,y,Xd,yd,r=args.relief_r)
        elif args.initial_policy=='proposal_only':
            initial={'banks':[],'source_policy':'numerical CROPS candidate edges only; adaptive proof-directed materialization'}
        elif args.initial_policy=='sparse4':
            from relief_initial import initialize as seed_start
            n=len(y)
            initial=seed_start(L,U,X,y,Xd,yd,sources=[(0,n),(0,n//2),(n//2,n),None])
        else:initial['banks']=shared_start(L,U,X,y,Xd,yd)
        if initial['banks']:
            nb=sum(b['bank_cpu']>0 for b in initial['banks'])
            rows=sum(b['bank_cpu'] for b in initial['banks'])/(nb*len(y))
            query=sum(b['strict_query_cpu'] for b in initial['banks'])/sum(b['queries'] for b in initial['banks'])
            costs=(rows,query)
        if args.method in ['hybrid','lookahead','reliever']:
            from float_proposal import FloatPELT
            oracle=FloatPELT(X,y)
            childcpu=official_crops(oracle,Path(args.out).with_suffix(''))
            edges=sorted({tuple(e) for c in oracle.calls for e in c['path']})
            initial['proposal_calls']=oracle.calls
            initial['numerical_proposal_fits']=len(oracle.cache)
            initial['candidate_strict_updates']=[direct(L,U,X,y,Xd,yd,e) for e in edges]
            if costs is None:
                updates=initial['candidate_strict_updates']
                costs=(sum(r['strict_build_cpu'] for r in updates)/sum(t-s for s,t in edges),
                       sum(r['strict_target_query_cpu'] for r in updates)/len(updates))
                initial['materialization_cost_model']={'row_cpu':costs[0],'query_cpu':costs[1],
                    'calibration':'actual strict candidate edge build/query CPU'}
            if theta_cache is not None:
                theta_cache.update({e:(theta,{'proposal_cached':True}) for e,(value,theta) in oracle.cache.items()})
    else:
        oracle=StrictPELT(L,U,X,y,Xd,yd)
        childcpu=official_crops(oracle,Path(args.out).with_suffix(''))
        initial={'calls':oracle.calls,'direct_cost_fits':len(oracle.cache),
                 'direct_history':list(oracle.cache.values())}
    if args.cost_model_from and costs is not None:
        prior=json.loads(gzip.open(args.cost_model_from,'rt').read() if args.cost_model_from.endswith('.gz') else Path(args.cost_model_from).read_text());bb=prior['initial']['banks']
        count=sum(b['bank_cpu']>0 for b in bb)
        costs=(sum(b['bank_cpu'] for b in bb)/(count*prior['input_shape'][0]),
               sum(b['strict_query_cpu'] for b in bb)/sum(b['queries'] for b in bb))
    acquirer=None
    if args.source_policy=='lazy_relief':
        from lazy_relief_policy import LazyRelief
        acquirer=LazyRelief(X,y,Xd,yd,costs,args.relief_r,args.pricing,args.entropy_temperature,
                            args.target_priority,args.closure_kernel,theta_cache)
    initialcpu=time.process_time()-t
    t=time.process_time()
    result=finish(L,U,X,y,Xd,yd,args.budget,args.method in ['adaptive','hybrid','lookahead','reliever'],costs,args.closure,args.method in ['lookahead','reliever'],args.bank_cap,args.pricing,args.entropy_temperature,args.source_pool,theta_cache,args.target_priority,closure_kernel=args.closure_kernel,acquirer=acquirer,source_priority=args.source_priority,obligation_weight=args.obligation_weight,charge_round_overhead=args.charge_round_overhead,closure_piece_h=closure_piece_h,bank_geometry=args.bank_geometry,acquisition_objective=args.acquisition_objective,joint_bank_pricing=args.joint_bank_pricing,certificate_bank_pricing=args.certificate_bank_pricing);finishcpu=time.process_time()-t
    joint_cpp_cpu=sum(h.get('joint_cpp_cpu',0.) for h in result['history'])
    total=time.process_time()-start+childcpu+joint_cpp_cpu
    result.update(method=args.method,dataset=args.dataset,row_start=args.start,initial=initial,
        total_cpu=total,python_cpu=total-childcpu,R_process_cpu=childcpu,joint_cpp_cpu=joint_cpp_cpu,
        preparation_cpu=prep,initial_evidence_and_search_cpu=initialcpu,
        refinement_and_final_proof_cpu=finishcpu,wall_seconds=time.perf_counter()-wall,
        input_shape=list(X.shape),beta_range=[.1,1.2],h=8,lambda_ridge=.1,
        superadditive_lower_closure=args.closure,
        fit_gradient_tolerance=1e-7,
        bank_cap=args.bank_cap,
        initial_policy=args.initial_policy,
        action_pricing=args.pricing,
        entropy_temperature_policy=args.entropy_temperature,
        source_pool_size=args.source_pool,
        source_policy=args.source_policy,reuse_proposal_fits=args.reuse_proposal_fits,source_priority=args.source_priority,obligation_weight=args.obligation_weight,charge_round_overhead=args.charge_round_overhead,acquisition_objective=args.acquisition_objective,certificate_bank_pricing=args.certificate_bank_pricing,
        lazy_relief_policy=(None if acquirer is None else acquirer.summary()),
        target_priority=args.target_priority,
        closure_kernel=args.closure_kernel,
        closure_geometry=args.closure_geometry,closure_piece_h=closure_piece_h,bank_geometry=args.bank_geometry,
        materialization_cost_model=(None if costs is None else {'row_cpu':costs[0],'query_cpu':costs[1],
            'source':args.cost_model_from or 'current initial materialization'}),
        initial_seed_fit_gradient_tolerance=(2e-11 if args.initial_policy=='full4' and args.method!='reliever' else 1e-7),
        arithmetic='70-digit directed Decimal, Q=1e12',
        runtime_scope='from frozen X/y load through initial evidence/search and complete final graph proof; includes official R process CPU; shared data/graph warmup and method pricing kernels warmed before this boundary; fresh-process metrics include all startup/imports/warmup')
    Path(args.out).write_text(json.dumps(result,indent=2))
    np.savez_compressed(Path(args.out).with_suffix('.npz'),low_int=L,high_int=U,X=X,y=y)
    print('FULL_RESULT',args.start,args.method,result['status'],'actions',result['actions'],
          'CPU',total,'RCPU',childcpu,flush=True)

if __name__=='__main__':main()
