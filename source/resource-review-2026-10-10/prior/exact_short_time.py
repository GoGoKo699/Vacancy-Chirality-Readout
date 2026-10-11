#!/usr/bin/env python3
"""Exact short-time calibration coefficients in Q(sqrt(2),sqrt(3),i)[u].

Integer powers of 2H avoid decimal trigonometry and square-root approximations.
The largest spatial/spin sector is the full 40-state N_up=2,N_down=1 block.
"""
from __future__ import annotations
import argparse,itertools,json,math
from fractions import Fraction
from pathlib import Path
import numpy as np

def product_index(i,j):
    common=i&j
    factor=(2 if common&1 else 1)*(3 if common&2 else 1)*(-1 if common&4 else 1)
    return i^j,factor

MUL=[[(product_index(i,j)) for j in range(8)] for i in range(8)]

def fieldmul(x,y):
    out=[0]*8
    for i,a in enumerate(x):
        if a:
            for j,b in enumerate(y):
                if b:
                    k,f=MUL[i][j];out[k]+=a*b*f
    return out

def fieldbasis_times(index,x):
    out=[0]*8
    for j,b in enumerate(x):
        k,f=MUL[index][j];out[k]+=b*f
    return out

def exact_coefficients(order=9):
    # Same occupation physics, built independently in the relevant fixed-spin block.
    occ=[]
    for ups in itertools.combinations_with_replacement(range(4),2):
        for down in range(4):
            n=[0]*8
            for j in ups:n[2*j]+=1
            n[2*down+1]=1;occ.append(tuple(n))
    lookup={n:i for i,n in enumerate(occ)}
    assert len(occ)==40
    starts=[]
    for ds in (3,2,1):
        n=[0]*8
        for j in (1,2,3):n[2*j+(1 if j==ds else 0)]=1
        starts.append(lookup[tuple(n)])
    edges=[];w=[]
    sqrts={1:(0,1),2:(1,1),3:(2,1),4:(0,2),6:(3,1)}
    for col,n in enumerate(occ):
        w.append(sum((n[2*j]+n[2*j+1])*(n[2*j]+n[2*j+1]-1)//2 for j in range(4)))
        for j in range(4):
            k=(j+1)%4
            for spin in range(2):
                # phi=-2pi/3: 2T_left=-sqrt(a)*sqrt3+i*sqrt(a).
                for src,dst,imag_sign in ((2*k+spin,2*j+spin,1),(2*j+spin,2*k+spin,-1)):
                    if not n[src]:continue
                    nn=list(n);nn[src]-=1;nn[dst]+=1
                    index,coef=sqrts[n[src]*(n[dst]+1)]
                    term=[0]*8
                    for basis,sgn in ((2,-1),(4,imag_sign)):
                        idx,fac=MUL[index][basis];term[idx]+=coef*sgn*fac
                    edges.append((lookup[tuple(nn)],col,[(i,x) for i,x in enumerate(term) if x]))
    # powers[n][degree][output,input,field coefficient] = (2H)^n S.
    init=np.zeros((1,40,3,8),dtype=object)
    for j,i in enumerate(starts):init[0,i,j,0]=1
    powers=[];poly=init
    for n in range(order+1):
        powers.append(poly[:,starts,:,:].copy())
        if n==order:break
        out=np.zeros((n+2,40,3,8),dtype=object)
        for row,col,terms in edges:
            for field,weight in terms:
                for j in range(8):
                    target,fac=MUL[field][j]
                    out[:n+1,row,:,target]+=weight*fac*poly[:,col,:,j]
        for state,weight in enumerate(w):
            if weight:out[1:,state,:,:]+=2*weight*poly[:,state,:,:]
        poly=out
    # Compute E_l*2^l*l! = sum binom(l,n) i^n (-i)^(l-n) G_n^dagger G_(l-n).
    def i_power(n):
        return ((0,1),(4,1),(0,-1),(4,-1))[n%4]
    records=[]
    for ell in range(order+1):
        coefs=np.zeros((ell+1,3,3,8),dtype=object)
        for n in range(ell+1):
            m=ell-n;phase_i,phase_s=i_power(n-m);factor=math.comb(ell,n)*phase_s
            for a,X in enumerate(powers[n]):
                for b,Y in enumerate(powers[m]):
                    for r in range(3):
                        for c in range(3):
                            value=[0]*8
                            for k in range(3):
                                xx=[(-v if j&4 else v) for j,v in enumerate(X[k,r])]
                                v=fieldmul(xx,list(Y[k,c]));value=[x+y for x,y in zip(value,v)]
                            v=fieldbasis_times(phase_i,value)
                            coefs[a+b,r,c,:]+=np.array(v,dtype=object)*factor
        # Reflection R swaps input/output indices 0 and 2.
        reverse=coefs[:,::-1,::-1,:]
        delta=coefs-reverse
        denom=2**ell*math.factorial(ell)
        hs=[];gs=[]
        for degree in range(ell+1):
            h=list(delta[degree,0,0])
            if any(h):hs.append({'u_power':degree,'field_coefficients':[str(Fraction(int(v),denom)) for v in h]})
            # chi_(0,1)=+i/4 in this basis (checked exactly by its 3-cycle).
            # Return both imaginary coefficients rather than rely on that sign convention.
            im=[delta[degree,0,1,j] for j in range(4,8)]
            if any(im):gs.append({'u_power':degree,'four_times_imag_coefficients':[str(Fraction(4*int(v),denom)) for v in im]})
        records.append({'tau_power':ell,'h_coefficients':hs,'g_coefficients':gs})
    return {'dimension':40,'max_order':order,'field_basis':['1','sqrt2','sqrt3','sqrt6','i','i*sqrt2','i*sqrt3','i*sqrt6'],'coefficients':records,'arithmetic':'Exact Python integers and Fraction; no floating-point amplitudes or exponentials.'}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.output.exists():p.error('Use a new output file.')
    record=exact_coefficients()
    # chi_(0,1) is +i/4 for the chosen descending down-spin positions.
    for item in record['coefficients']:
        if item['tau_power']<9:assert not item['h_coefficients'],item
    h=record['coefficients'][9]['h_coefficients']
    assert h==[{'u_power':1,'field_coefficients':['0','0','7/72','0','0','0','0','0']},
               {'u_power':3,'field_coefficients':['0','0','1/2520','0','0','0','0','0']}],h
    g4=record['coefficients'][4]['g_coefficients']
    assert g4==[{'u_power':0,'four_times_imag_coefficients':['0','0','4','0']}],g4
    for item in record['coefficients'][:4]:assert not item['g_coefficients'],item
    record.update(status='PASS',test_groups=1)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(record,f,indent=2,sort_keys=True);f.write('\n')
    print('Exact short-time coefficient identity: PASS')
if __name__=='__main__':main()
