#include <bits/stdc++.h>
using namespace std; using i64=long long; using i128=__int128_t;
struct Rat {i128 x; int d;};
bool lt(Rat a,Rat b){return a.x*b.d<b.x*a.d;}
string s128(i128 v){if(!v)return "0";bool neg=v<0;if(neg)v=-v;string s;while(v){s+=char('0'+v%10);v/=10;}if(neg)s+='-';reverse(s.begin(),s.end());return s;}
string fmt(Rat r){i128 a=r.x;int b=r.d;if(a==0)return "0";i64 aa=(i64)(a<0?-a:a);i64 g=std::gcd(aa,(i64)b); a/=g; b/=g; return s128(a)+(b==1?"":"/"+to_string(b));}
int main(int argc,char**argv){
 if(argc<9){cerr<<"usage n h b0 b1 L.bin U.bin candidates.log results.txt\n";return 2;}
 int n=atoi(argv[1]),h=atoi(argv[2]),N=n+1, K=n/h; i64 b0=atoll(argv[3]),b1=atoll(argv[4]);
 vector<i64>L((size_t)N*N),U((size_t)N*N);
 ifstream lf(argv[5],ios::binary),uf(argv[6],ios::binary);lf.read((char*)L.data(),L.size()*8);uf.read((char*)U.data(),U.size()*8);if(!lf || !uf)return 3;
 vector<vector<pair<int,int>>> paths; ifstream in(argv[7]);string line;
 while(getline(in,line))if(line.rfind("PATH ",0)==0){istringstream ss(line.substr(5));string p;vector<pair<int,int>> edges;while(ss>>p){auto i=p.find(':');if(i!=string::npos)edges.push_back({stoi(p.substr(0,i)),stoi(p.substr(i+1))});}if(!edges.empty())paths.push_back(move(edges));}
 ofstream out(argv[8]);int pass=0;const i64 INF=LLONG_MAX/4;vector<unsigned char>inP((size_t)N*N,0);
 auto start=chrono::steady_clock::now();
 for(int index=0;index<(int)paths.size();++index){
  const auto &P=paths[index];fill(inP.begin(),inP.end(),0);i64 A=0;int kp=P.size();
  for(auto [s,t]:P){inP[(size_t)s*N+t]=1;A+=L[(size_t)s*N+t];}
  vector<i64>prev(n+1,INF),cur(n+1,INF),m(K+1,INF);prev[0]=0;
  for(int k=1;k<=K;++k){fill(cur.begin(),cur.end(),INF);
   for(int t=k*h;t<=n;++t){
    i64 best=INF;int s0=(k-1)*h;for(int s=s0;s<=t-h;++s){i64 z=prev[s];if(z==INF)continue;
     size_t e=(size_t)s*N+t;i64 v=z+(inP[e]?L[e]:U[e]);if(v<best)best=v;}
    cur[t]=best;
   }
   m[k]=cur[n];prev.swap(cur);
  }
  Rat lo{b0,1},hi{b1,1};bool same_count_ok=m[kp]==A;
  for(int k=1;k<=K;++k){if(k==kp||m[k]==INF)continue;
    if(k>kp){Rat r{(i128)A-m[k],k-kp};if(lt(lo,r))lo=r;}
    else{Rat r{(i128)m[k]-A,kp-k};if(lt(r,hi))hi=r;}
  }
  bool possible=same_count_ok && !lt(hi,lo);pass+=possible;
  out<<"INDEX "<<index<<" POSSIBLE "<<(possible?1:0)<<" K "<<kp<<" LO "<<fmt(lo)<<" HI "<<fmt(hi)<<" SAME_COUNT "<<(same_count_ok?1:0)<<"\n";
 }
 auto stop=chrono::steady_clock::now();
 out<<"TOTAL "<<paths.size()<<" POSSIBLE "<<pass<<" EXCLUDED "<<(paths.size()-pass)<<" DP_SECONDS "<<chrono::duration<double>(stop-start).count()<<"\n";
 cerr<<"TOTAL="<<paths.size()<<" POSSIBLE="<<pass<<" DP_SECONDS="<<chrono::duration<double>(stop-start).count()<<"\n";
}
