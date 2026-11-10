
#include <bits/stdc++.h>
using namespace std;
using i64 = long long;
using i128 = __int128_t;
struct Edge { int to; int id; i64 l,u; };
struct Label { int node; vector<int> es; i64 sumL=0,sumU=0; bool active=true, accepted=false; };
struct Item { i64 pri; long long seq; int id; };
struct Cmp { bool operator()(const Item&a,const Item&b)const{
    if(a.pri!=b.pri) return a.pri>b.pri;
    return a.seq>b.seq;
}};
static inline i64 shared_width(const vector<int>&a,const vector<int>&b,
                               const vector<i64>&width){
    size_t i=0,j=0; i64 s=0;
    while(i<a.size() && j<b.size()){
        if(a[i]==b[j]){s+=width[a[i]]; ++i;++j;}
        else if(a[i]<b[j]) ++i;
        else ++j;
    }
    return s;
}
static inline bool dominates(const Label&q,const Label&p,i64 b0,i64 b1,
                             const vector<i64>&width){
    // max_{x,beta} J_q-J_p = sum U(q)-sum L(p)-sum_{shared}(U-L)+(kq-kp)beta
    i64 base=q.sumU-p.sumL-shared_width(q.es,p.es,width);
    i64 dk=(i64)q.es.size()-(i64)p.es.size();
    i64 d0=base+dk*b0, d1=base+dk*b1;
    return max(d0,d1)<0;
}

struct Witness{
    vector<unsigned char> has;
    vector<vector<i128>> tail; // [beta grid point][node]
    i128 sumU=0;
    int k=0;
};
static inline bool archive_prunes(const Label& p,int v,const vector<Witness>& witnesses,
        const vector<i64>&width,const vector<i64>&betas,long long&checked,
        long long&covered_cells){
    const int G=(int)betas.size()-1;
    vector<unsigned char> covered(G,0);
    int remaining=G;
    for(const Witness &w:witnesses){
        ++checked;
        i128 prefix=(i128)p.sumL;
        for(int e:p.es) if(w.has[e]) prefix+=(i128)width[e];
        vector<i128> val(G+1);
        for(int j=0;j<=G;++j){
            i128 b=betas[j];
            val[j]=prefix + (i128)p.es.size()*b + w.tail[j][v]
                   - w.sumU - (i128)w.k*b;
        }
        // For each cell, one witness positive at both endpoints certifies
        // strict exclusion throughout that entire cell by concavity.
        for(int j=0;j<G;++j){
            if(!covered[j] && val[j]>0 && val[j+1]>0){
                covered[j]=1;
                --remaining;
                ++covered_cells;
            }
        }
        if(remaining==0) return true;
    }
    return false;
}

int main(int argc,char**argv){
    if(argc<9){cerr<<"usage n h beta0 beta1 L.bin U.bin max_created time_limit_sec\n";return 2;}
    int n=atoi(argv[1]),h=atoi(argv[2]); i64 b0=atoll(argv[3]),b1=atoll(argv[4]);
    long long max_created=atoll(argv[7]); double time_limit=atof(argv[8]);
    int grid_cells=(argc>10?max(1,atoi(argv[10])):8);
    if(grid_cells>128) grid_cells=128;
    vector<i64> betas(grid_cells+1);
    for(int j=0;j<=grid_cells;++j)
        betas[j]=(i64)((i128)b0 + ((i128)(b1-b0)*j)/grid_cells);
    int N=n+1, M=N*N;
    vector<i64>L(M),U(M),width(M);
    ifstream lf(argv[5],ios::binary),uf(argv[6],ios::binary);
    lf.read((char*)L.data(),M*sizeof(i64)); uf.read((char*)U.data(),M*sizeof(i64));
    if(!lf||!uf){cerr<<"matrix read failed\n";return 3;}
    vector<vector<Edge>> adj(N);
    long long m=0;
    for(int s=0;s<=n-h;s++) for(int t=s+h;t<=n;t++){
        int id=s*N+t;
        adj[s].push_back({t,id,L[id],U[id]}); width[id]=U[id]-L[id]; ++m;
    }

    vector<Witness> witnesses;
    if(argc>9){
       ifstream wf(argv[9]); string ln;
       while(getline(wf,ln)){
          if(ln.empty()) continue;
          Witness w; w.has.assign(M,0);
          istringstream is(ln);string tok;
          while(is>>tok){
             auto z=tok.find(':');if(z==string::npos)continue;
             int a=stoi(tok.substr(0,z)),b=stoi(tok.substr(z+1));
             int e=a*N+b;w.has[e]=1;w.sumU+=(i128)U[e];++w.k;
          }
          if(w.k==0)continue;
          const i128 INF=((i128)1<<110);
          w.tail.assign(grid_cells+1,vector<i128>(N,INF));
          for(int j=0;j<=grid_cells;++j){
             w.tail[j][n]=0;
             i128 beta=betas[j];
             for(int v=n-1;v>=0;--v){
                for(const Edge &e:adj[v]){
                   i128 c=w.has[e.id]? (i128)U[e.id]:(i128)L[e.id];
                   w.tail[j][v]=min(w.tail[j][v],c+beta+w.tail[j][e.to]);
                }
             }
          }
          witnesses.push_back(std::move(w));
       }
    }
    long long archived_pruned=0, archived_checked=0, archived_cells_covered=0;
    vector<Label> labs; labs.reserve((size_t)max_created+1);
    vector<vector<int>> at(N);
    priority_queue<Item,vector<Item>,Cmp> pq;
    labs.push_back({0,{},0,0,true,false});
    pq.push({0,0,0});
    long long seq=0,pops=0,extensions=0,skipped=0,accepted=0,max_at=1;
    auto start=chrono::steady_clock::now();
    string status="COMPLETE";
    while(!pq.empty()){
        auto now=chrono::steady_clock::now();
        double elapsed=chrono::duration<double>(now-start).count();
        if(elapsed>time_limit){status="TIME_LIMIT";break;}
        Item it=pq.top();pq.pop(); ++pops;
        Label &p=labs[it.id];
        if(!p.active) continue;
        bool dom=false;
        for(int qid:at[p.node]){
            const Label&q=labs[qid];
            if(q.active && q.accepted && dominates(q,p,b0,b1,width)){dom=true;break;}
        }
        if(dom){p.active=false;++skipped;continue;}
        // Remove any older labels dominated by p (should be rare under label-setting order).
        for(int qid:at[p.node]){
            Label &q=labs[qid];
            if(q.active && q.accepted && dominates(p,q,b0,b1,width)) q.active=false;
        }
        p.accepted=true; at[p.node].push_back(it.id); ++accepted;
        if((long long)at[p.node].size()>max_at) max_at=at[p.node].size();
        int v=p.node;
        if(v==n) continue;
        for(const Edge&e:adj[v]){
            Label c;
            c.node=e.to; c.es=p.es; c.es.push_back(e.id);
            c.sumL=p.sumL+e.l; c.sumU=p.sumU+e.u;
            // Cheaply prune against already accepted labels at target.
            bool cd=false;
            for(int qid:at[e.to]){
                const Label&q=labs[qid];
                if(q.active && q.accepted && dominates(q,c,b0,b1,width)){cd=true;break;}
            }
            if(cd){++skipped;continue;}
            if(!witnesses.empty() && archive_prunes(c,e.to,witnesses,width,betas,archived_checked,archived_cells_covered)){
                ++archived_pruned; ++skipped; continue;
            }
            if((long long)labs.size()>=max_created){status="LABEL_LIMIT";goto done;}
            int id=labs.size(); labs.push_back(std::move(c));
            ++extensions; ++seq;
            pq.push({labs[id].sumL+(i64)labs[id].es.size()*b0,seq,id});
        }
    }
done:
    {
      auto now=chrono::steady_clock::now();
      double elapsed=chrono::duration<double>(now-start).count();
      vector<int> finalids;
      for(int id:at[n]) if(labs[id].active && labs[id].accepted) finalids.push_back(id);
      cout<<"status "<<status<<"\n";
      cout<<"n "<<n<<" h "<<h<<" edges "<<m<<"\n";
      cout<<"pops "<<pops<<" extensions "<<extensions<<" created "<<labs.size()<<" skipped "<<skipped<<" accepted "<<accepted<<"\n";
      cout<<"max_labels_at_node "<<max_at<<" active_sink_paths "<<finalids.size()<<" seconds "<<elapsed<<"\n";
      cout<<"archive_witnesses "<<witnesses.size()<<" archive_pruned "<<archived_pruned<<" archive_checked "<<archived_checked<<" beta_cells "<<grid_cells<<" archive_cells_covered "<<archived_cells_covered<<"\n";
      for(int id:finalids){
        cout<<"PATH";
        for(int e:labs[id].es) cout<<" "<<(e/N)<<":"<<(e%N);
        cout<<"\n";
      }
    }
}
