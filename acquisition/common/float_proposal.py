import time
import numpy as np
from safe_archive import fit

class FloatPELT:
    """Numerical proposal only; final certification uses directed evidence."""
    def __init__(self,X,y):self.X=X;self.y=y;self.cache={};self.calls=[]
    def cost(self,s,t):
        if (s,t) not in self.cache:
            theta,rec=fit(self.X,self.y,s,t,tol=1e-7)
            z=self.X[s:t]@theta
            value=float(np.logaddexp(0,z).sum()-self.y[s:t]@z+.05*(t-s)*(theta@theta))
            self.cache[s,t]=(value,theta)
        return self.cache[s,t][0]
    def solve(self,beta):
        start=time.process_time();n=len(self.y);h=8
        D=np.full(n+1,np.inf);D[0]=0.;pred=np.full(n+1,-1,np.int32)
        active=[0];expiry={0:n+1};evals=0
        for t in range(h,n+1):
            active=[s for s in active if expiry[s]>t]
            eligible=[s for s in active if s<=t-h]
            costs=[self.cost(s,t) for s in eligible]
            vals=[D[s]+c+beta for s,c in zip(eligible,costs)]
            j=int(np.argmin(vals));D[t]=vals[j];pred[t]=eligible[j];evals+=len(eligible)
            for s,c in zip(eligible,costs):
                if D[s]+c>D[t]+1e-9:expiry[s]=min(expiry[s],t+h)
            if t<n:active.append(t);expiry[t]=n+1
        P=[];t=n
        while t:
            s=int(pred[t]);P.append((s,t));t=s
        P=P[::-1]
        rec={'beta':beta,'path':P,'unpenalized':sum(self.cost(*e) for e in P),
             'ncp':len(P)-1,'edge_evals':evals,'cpu':time.process_time()-start}
        self.calls.append(rec);print('FLOAT_PELT',beta,len(P),len(self.cache),flush=True)
        return rec
