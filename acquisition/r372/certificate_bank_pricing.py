"""R372 research pilot: exact *projected* archive-complement debt scoring.

All results are forecasts conditional on a menu's projected cost box;
actual Decimal70 evidence and real certification stay in the caller.
Preserves only current online archive, retirement and cover state.
"""
from __future__ import annotations
import time
import numpy as np
from safe_archive import SafeArchive

def shadow_current_archive(arch):
    out=SafeArchive(arch.n,arch.h,arch.b0,arch.b1,arch.mode,False)
    out.A=dict(arch.A)
    out.cover=dict(arch.cover)
    out.retired=dict(arch.retired)
    out.previous=set(arch.previous)
    # fresh node trie is required if no new path in projected hull
    return out

def priced_projected_certificate_debt(arch,alternatives):
    """Return full projected debt for every nested bank using one independent shadow per menu.

    Complexity per menu ≈ cost of one exact hull plus active archived-path
    complement/retirement calls. Worst case grows with archive size; no
    speed guarantee outside measured inputs.
    """
    t0=time.process_time(); out=[]
    for alt in alternatives:
        t1=time.process_time()
        s=shadow_current_archive(arch)
        C,new,_=s.refresh(alt['L'],alt['U'])
        st=s.check(alt['L'],alt['U'],C)
        debt=sum(max(0.,-float(v['margin'])) for v in st['failures'])
        out.append({'debt':debt,'failures':len(st['failures']),
                    'archived_paths':len(s.A),'new_paths':len(new),
                    'cpu':time.process_time()-t1})
    return out, time.process_time()-t0

def price_by_debt(arch,base_failures,alternatives,round_cpu=0.):
    base_debt=sum(max(0.,-float(v['margin'])) for v in base_failures)
    rows, cpu=priced_projected_certificate_debt(arch,alternatives)
    # Completion of current proof requirements is primary; debt-size marginal
    # progress per action cost is a meaningful tie-break for noncompletion.
    for i,(row,alt) in enumerate(zip(rows,alternatives)):
        row['gain']=max(0.,base_debt-row['debt'])
        row['score']=row['gain']/max(1e-10,float(alt['cost'])+round_cpu+cpu)
    choice=max(range(len(rows)),key=lambda i:(rows[i]['score'],-alternatives[i]['cost']))
    return choice, rows, {'policy':'certificate_debt','pricing_cpu':cpu,'base_debt':base_debt}
