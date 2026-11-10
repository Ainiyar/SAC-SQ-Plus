"""Exact possible-optimal penalty support relative to a finite path archive.

All pairwise inequalities use the SAME favorable scenario for R: L on R,
U outside R. Thus their intersection is exact, including singleton supports.
"""
from fractions import Fraction as F
from safe_archive import diff_lower

def support(R,archive,L,U,b0,b1):
    lo,hi=F(b0),F(b1);lower_witness=None;upper_witness=None
    for j,P in enumerate(archive):
        a,b=diff_lower(R,P,L,U)
        # In the R-favorable scenario, J_R-J_P = a+b beta.
        # R can win exactly when a+b beta <= 0 for every P.
        if b==0:
            if a>0:return {'empty':True,'reason':'same_count','witness':j,'gap':str(a)}
        elif b>0:
            x=F(-a,b)
            if x<hi:hi=x;upper_witness=j
        else:
            x=F(-a,b)
            if x>lo:lo=x;lower_witness=j
        if lo>hi:
            return {'empty':True,'reason':'crossed_support','lo':str(lo),'hi':str(hi),
                    'lower_witness':lower_witness,'upper_witness':upper_witness}
    return {'empty':False,'lo':str(lo),'hi':str(hi),'singleton':lo==hi,
            'lower_witness':lower_witness,'upper_witness':upper_witness}

def all_supports(archive,L,U,b0,b1):
    archive=list(archive)
    return [{'id':i,'segments':len(R),'support':support(R,archive,L,U,b0,b1)}
            for i,R in enumerate(archive)]
