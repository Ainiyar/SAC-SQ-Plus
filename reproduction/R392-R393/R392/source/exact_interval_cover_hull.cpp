#include <bits/stdc++.h>
using namespace std;
using i64 = long long;
using i128 = __int128_t;

struct Edge { int to; int id; i64 l,u; };
struct Label { int node; vector<int> es; i64 sumL=0,sumU=0; bool active=true, accepted=false; };
struct Item { i64 pri; long long seq; int id; };
struct Cmp { bool operator()(const Item&a,const Item&b)const { if(a.pri!=b.pri) return a.pri>b.pri; return a.seq>b.seq; } };

struct Rat { i128 num=0, den=1; bool neg_inf=false, pos_inf=false; };
static inline int cmp_rat(const Rat&a,const Rat&b){
    if(a.neg_inf) return b.neg_inf?0:-1;
    if(a.pos_inf) return b.pos_inf?0:1;
    if(b.neg_inf) return 1;
    if(b.pos_inf) return -1;
    i128 x=a.num*b.den, y=b.num*a.den;
    return (x>y)-(x<y);
}
static inline int cmp_rat_i64(const Rat&a,i64 b){
    if(a.neg_inf) return -1;
    if(a.pos_inf) return 1;
    i128 x=a.num, y=(i128)b*a.den;
    return (x>y)-(x<y);
}
static inline Rat make_rat(i128 n,i128 d){ if(d<0){n=-n;d=-d;} return Rat{n,d,false,false}; }
static inline string i128str(i128 x){
    if(x==0) return "0";
    bool neg=x<0; if(neg)x=-x;
    string s; while(x){s.push_back('0'+(int)(x%10));x/=10;}
    if(neg)s.push_back('-'); reverse(s.begin(),s.end()); return s;
}
struct Line { int m; i128 b; }; // b + m*beta
struct HullLine { int m; i128 b; };
struct Witness {
    vector<unsigned char> has;
    vector<vector<HullLine>> hull; // exact lower envelope over [beta0,beta1]
    i128 sumU=0;
    int k=0;
};
struct Interval { Rat lo,hi; };

static vector<HullLine> make_hull_sorted(const vector<Line>& uniq, i64 beta0, i64 beta1){
    if(uniq.empty()) return {};
    vector<Line> hs; vector<Rat> starts; hs.reserve(uniq.size()); starts.reserve(uniq.size());
    for(const auto &x:uniq){
        if(hs.empty()){
            hs.push_back(x); starts.push_back(Rat{0,1,true,false}); continue;
        }
        Rat cross;
        while(true){
            const auto &last=hs.back();
            cross=make_rat(x.b-last.b,(i128)last.m-x.m);
            if(hs.size()==1 || cmp_rat(cross,starts.back())>0) break;
            hs.pop_back(); starts.pop_back();
        }
        hs.push_back(x); starts.push_back(cross);
    }
    vector<HullLine> out; out.reserve(hs.size());
    for(size_t i=0;i<hs.size();++i){
        bool starts_before_end = cmp_rat_i64(starts[i],beta1)<=0;
        bool ends_after_start = (i+1==hs.size()) || cmp_rat_i64(starts[i+1],beta0)>=0;
        if(starts_before_end && ends_after_start) out.push_back({hs[i].m,hs[i].b});
    }
    return out;
}

static inline bool witness_interval(const i64 prefix_mix, const Witness&w, int prefix_k,
                                    int node, i64 beta0,i64 beta1, Interval &out){
    const auto &h=w.hull[node];
    if(h.empty()) return false;
    Rat lo; lo.neg_inf=true;
    Rat hi; hi.pos_inf=true;
    for(const auto &ln:h){
        i128 A=(i128)prefix_mix + ln.b - w.sumU;
        int d=prefix_k + ln.m - w.k;
        if(d>0){
            Rat r=make_rat(-A,d);
            if(cmp_rat(r,lo)>0) lo=r;
        }else if(d<0){
            Rat r=make_rat(A,-d);
            if(cmp_rat(r,hi)<0) hi=r;
        }else if(A<=0){
            return false;
        }
    }
    if(cmp_rat(lo,hi)>=0) return false;
    // The interval must intersect the beta domain in a way that can cover a point.
    if(cmp_rat_i64(lo,beta1)>=0 || cmp_rat_i64(hi,beta0)<=0) return false;
    out={lo,hi};
    return true;
}

static inline bool union_covers_closed_domain(vector<Interval>& iv, i64 beta0,i64 beta1){
    if(iv.empty()) return false;
    sort(iv.begin(),iv.end(),[](const Interval&a,const Interval&b){
        int c=cmp_rat(a.lo,b.lo);
        if(c) return c<0;
        return cmp_rat(a.hi,b.hi)>0;
    });
    bool started=false; Rat reach;
    for(const auto &x:iv){
        if(!started){
            if(cmp_rat_i64(x.lo,beta0)<0 && cmp_rat_i64(x.hi,beta0)>0){
                reach=x.hi; started=true;
                if(cmp_rat_i64(reach,beta1)>0) return true;
            }
        }else{
            if(cmp_rat(x.lo,reach)<0 && cmp_rat(x.hi,reach)>0){
                if(cmp_rat(x.hi,reach)>0) reach=x.hi;
                if(cmp_rat_i64(reach,beta1)>0) return true;
            }
            // Sorted starts: once the next interval starts at or beyond reach,
            // no later interval can bridge the uncovered point at reach.
            else if(cmp_rat(x.lo,reach)>=0 && cmp_rat_i64(reach,beta1)<=0){
                break;
            }
        }
    }
    return started && cmp_rat_i64(reach,beta1)>0;
}

static inline i64 mixed_edge_cost(const Witness&w,const vector<i64>&L,const vector<i64>&U,int id){
    return w.has[id]?U[id]:L[id];
}

static bool archive_prunes(const Label&p,int v,const vector<Witness>&ws,
                           const i64 *mix,int W,
                           i64 beta0,i64 beta1,long long&checked,long long&intervals_created){
    if(ws.empty()) return false;
    // No legal suffix means the prefix has no complete path at all.
    bool any_suffix=false;
    for(const auto&w:ws) if(!w.hull[v].empty()){any_suffix=true;break;}
    if(!any_suffix) return true;
    vector<Interval> iv; iv.reserve(W);
    for(int j=0;j<W;++j){
        ++checked;
        Interval x;
        if(witness_interval(mix[j],ws[j],(int)p.es.size(),v,beta0,beta1,x)){
            iv.push_back(x); ++intervals_created;
        }
    }
    return union_covers_closed_domain(iv,beta0,beta1);
}

int main(int argc,char**argv){
    if(argc<9){cerr<<"usage n h beta0 beta1 L.bin U.bin max_created time_limit_sec [witnesses.txt]\n";return 2;}
    int n=atoi(argv[1]),h=atoi(argv[2]); i64 beta0=atoll(argv[3]),beta1=atoll(argv[4]);
    long long max_created=atoll(argv[7]); double time_limit=atof(argv[8]);
    int N=n+1,M=N*N;
    vector<i64>L(M),U(M),width(M);
    ifstream lf(argv[5],ios::binary),uf(argv[6],ios::binary);
    lf.read((char*)L.data(),M*sizeof(i64)); uf.read((char*)U.data(),M*sizeof(i64));
    if(!lf||!uf){cerr<<"matrix read failed\n";return 3;}
    vector<vector<Edge>> adj(N); long long edge_count=0;
    for(int s=0;s<=n-h;s++) for(int t=s+h;t<=n;t++){
        int id=s*N+t; adj[s].push_back({t,id,L[id],U[id]}); width[id]=U[id]-L[id]; ++edge_count;
    }
    vector<Witness> ws;
    if(argc>9){
        ifstream wf(argv[9]); string ln;
        while(getline(wf,ln)){
            if(ln.empty()) continue;
            Witness w; w.has.assign(M,0);
            istringstream is(ln); string tok;
            while(is>>tok){auto z=tok.find(':');if(z==string::npos)continue;int a=stoi(tok.substr(0,z)),b=stoi(tok.substr(z+1));int id=a*N+b;w.has[id]=1;w.sumU+=(i128)U[id];++w.k;}
            if(w.k==0) continue;
            ws.push_back(std::move(w));
        }
    }
    if(ws.empty()){cerr<<"no witnesses\n";return 4;}

    double prep_cpu0=clock(); auto prep_wall0=chrono::steady_clock::now();
    long long hull_candidates=0,hull_vertices=0,max_hull=0;
    const int K=n/h;
    const i128 INF_HULL=((i128)1<<120);
    for(auto &w:ws){
        w.hull.resize(N); w.hull[n].push_back({0,0});
        vector<i128> best(K+1,INF_HULL);
        vector<Line> cand; cand.reserve(K+1);
        for(int v=n-1;v>=0;--v){
            fill(best.begin(),best.end(),INF_HULL);
            for(const auto&e:adj[v]){
                i64 c=w.has[e.id]?U[e.id]:L[e.id];
                for(const auto&ln:w.hull[e.to]){
                    int m=ln.m+1; i128 b=ln.b+(i128)c;
                    if(b<best[m]) best[m]=b;
                    ++hull_candidates;
                }
            }
            cand.clear();
            for(int m=K;m>=0;--m) if(best[m]!=INF_HULL) cand.push_back({m,best[m]});
            w.hull[v]=make_hull_sorted(cand,beta0,beta1);
            hull_vertices+=(long long)w.hull[v].size();
            max_hull=max<long long>(max_hull,w.hull[v].size());
        }
    }
    double prep_cpu=(clock()-prep_cpu0)/CLOCKS_PER_SEC;
    double prep_wall=chrono::duration<double>(chrono::steady_clock::now()-prep_wall0).count();

    vector<Label> labs; labs.reserve((size_t)max_created+1);
    vector<i64> mixflat; mixflat.reserve(((size_t)max_created+1)*ws.size());
    vector<vector<int>> at(N);
    priority_queue<Item,vector<Item>,Cmp> pq;
    labs.push_back({0,{},0,0,true,false}); mixflat.resize(ws.size(),0);
    pq.push({0,0,0});
    long long seq=0,pops=0,extensions=0,skipped=0,accepted=0,max_at=1;
    long long archived_pruned=0,archived_checked=0,intervals_created=0;
    double search_cpu0=clock(); auto search_wall0=chrono::steady_clock::now();
    string status="COMPLETE";
    while(!pq.empty()){
        double elapsed=(clock()-search_cpu0)/CLOCKS_PER_SEC;
        if(elapsed>time_limit){status="TIME_LIMIT";break;}
        Item it=pq.top();pq.pop(); ++pops;
        Label &p=labs[it.id];
        if(!p.active) continue;
        bool dom=false;
        for(int qid:at[p.node]){
            const Label&q=labs[qid];
            if(q.active&&q.accepted){
                i64 shared=0; size_t i=0,j=0;
                while(i<q.es.size()&&j<p.es.size()){
                    if(q.es[i]==p.es[j]){shared+=width[q.es[i]];++i;++j;}
                    else if(q.es[i]<p.es[j])++i;else ++j;
                }
                i128 base=(i128)q.sumU-p.sumL-shared;
                i128 dk=(i128)q.es.size()-(i128)p.es.size();
                i128 d0=base+dk*beta0,d1=base+dk*beta1;
                if(max(d0,d1)<0){dom=true;break;}
            }
        }
        if(dom){p.active=false;++skipped;continue;}
        for(int qid:at[p.node]){
            Label&q=labs[qid];
            if(q.active&&q.accepted){
                i64 shared=0;size_t i=0,j=0;
                while(i<p.es.size()&&j<q.es.size()){
                    if(p.es[i]==q.es[j]){shared+=width[p.es[i]];++i;++j;}
                    else if(p.es[i]<q.es[j])++i;else ++j;
                }
                i128 base=(i128)p.sumU-q.sumL-shared;
                i128 dk=(i128)p.es.size()-(i128)q.es.size();
                i128 d0=base+dk*beta0,d1=base+dk*beta1;
                if(max(d0,d1)<0) q.active=false;
            }
        }
        p.accepted=true;at[p.node].push_back(it.id);++accepted;
        if((long long)at[p.node].size()>max_at)max_at=at[p.node].size();
        int v=p.node;if(v==n)continue;
        const i64 *pmix=&mixflat[(size_t)it.id*ws.size()];
        vector<i64> childmix(ws.size());
        for(const Edge&e:adj[v]){
            Label c;c.node=e.to;c.es=p.es;c.es.push_back(e.id);c.sumL=p.sumL+e.l;c.sumU=p.sumU+e.u;
            bool cd=false;
            for(int qid:at[e.to]){
                const Label&q=labs[qid];if(!q.active||!q.accepted)continue;
                i64 shared=0;size_t i=0,j=0;
                while(i<q.es.size()&&j<c.es.size()){
                    if(q.es[i]==c.es[j]){shared+=width[q.es[i]];++i;++j;}
                    else if(q.es[i]<c.es[j])++i;else ++j;
                }
                i128 base=(i128)q.sumU-c.sumL-shared;
                i128 dk=(i128)q.es.size()-(i128)c.es.size();
                i128 d0=base+dk*beta0,d1=base+dk*beta1;
                if(max(d0,d1)<0){cd=true;break;}
            }
            if(cd){++skipped;continue;}
            for(size_t j=0;j<ws.size();++j) childmix[j]=pmix[j]+(ws[j].has[e.id]?e.u:e.l);
            bool prune=archive_prunes(c,e.to,ws,childmix.data(),(int)ws.size(),beta0,beta1,archived_checked,intervals_created);
            if(prune){++archived_pruned;++skipped;continue;}
            if((long long)labs.size()>=max_created){status="LABEL_LIMIT";goto done;}
            int id=labs.size();labs.push_back(std::move(c));
            mixflat.insert(mixflat.end(),childmix.begin(),childmix.end());
            ++extensions;++seq;pq.push({labs[id].sumL+(i64)labs[id].es.size()*beta0,seq,id});
        }
    }

done:{
        double search_cpu=(clock()-search_cpu0)/CLOCKS_PER_SEC;
        double search_wall=chrono::duration<double>(chrono::steady_clock::now()-search_wall0).count();
        vector<int> finalids;for(int id:at[n])if(labs[id].active&&labs[id].accepted)finalids.push_back(id);
        cout<<"status "<<status<<"\n";
        cout<<"n "<<n<<" h "<<h<<" edges "<<edge_count<<"\n";
        cout<<"pops "<<pops<<" extensions "<<extensions<<" created "<<labs.size()<<" skipped "<<skipped<<" accepted "<<accepted<<"\n";
        cout<<"max_labels_at_node "<<max_at<<" active_sink_paths "<<finalids.size()<<"\n";
        cout<<"archive_witnesses "<<ws.size()<<" archive_pruned "<<archived_pruned<<" archive_checked "<<archived_checked<<" intervals_created "<<intervals_created<<"\n";
        cout<<"hull_candidates "<<hull_candidates<<" hull_vertices_total "<<hull_vertices<<" max_hull_vertices_at_node "<<max_hull<<"\n";
        cout<<fixed<<setprecision(6)<<"preprocess_cpu_seconds "<<prep_cpu<<" preprocess_wall_seconds "<<prep_wall<<"\n";
        cout<<"search_cpu_seconds "<<search_cpu<<" search_wall_seconds "<<search_wall<<" total_cpu_seconds "<<(prep_cpu+search_cpu)<<"\n";
        for(int id:finalids){cout<<"PATH";for(int e:labs[id].es)cout<<" "<<(e/N)<<":"<<(e%N);cout<<"\n";}
    }
    return 0;
}
