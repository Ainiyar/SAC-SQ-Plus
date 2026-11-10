import io,zipfile
from pathlib import Path
import numpy as np

def raw():
    import pandas as pd
    with zipfile.ZipFile('secom/uci_secom.zip') as z:
        dm=next(x for x in z.namelist() if Path(x).name=='secom.data')
        lm=next(x for x in z.namelist() if Path(x).name=='secom_labels.data')
        a=np.loadtxt(io.StringIO(z.read(dm).decode('utf8')))
        lines=[s.split(maxsplit=1) for s in z.read(lm).decode('utf8').splitlines() if s.strip()]
    labels=np.array([int(s[0]) for s in lines]);stamp=pd.to_datetime([s[1].strip().strip('"') for s in lines],format='%d/%m/%Y %H:%M:%S')
    order=np.argsort(stamp.asi8,kind='stable')
    missing=np.mean(~np.isfinite(a),axis=0);cols=[]
    for j in range(a.shape[1]):
        if missing[j]<=.05 and np.nanstd(a[:,j])>1e-12:cols.append(j)
        if len(cols)==8:break
    if len(cols)!=8:raise ValueError('fewer than eight eligible sensors')
    return a[order],(labels[order]==1).astype(np.int64),stamp[order],cols,order

def block(a,y,start,n,cols):
    b=a[start:start+n,cols].copy();med=np.nanmedian(b,axis=0)
    bad=~np.isfinite(b);b[bad]=np.broadcast_to(med,b.shape)[bad]
    sd=b.std(axis=0);sd[sd<1e-12]=1.
    X=np.column_stack([np.ones(len(b)),(b-b.mean(axis=0))/sd,np.linspace(-1,1,len(b))])
    return X,y[start:start+n]

def preprocess(start,n=128):
    a,y,stamp,cols,order=raw()
    return block(a,y,start,n,cols)
