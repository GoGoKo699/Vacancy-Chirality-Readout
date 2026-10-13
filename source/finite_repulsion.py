"""Finite-repulsion vacancy-return instrument, with no extra spin controls.

Units: J=1, tau=Jt, u=U/J. Numerical values are diagnostics,
not interval certificates.
"""
from __future__ import annotations
from functools import lru_cache
import itertools
import math
import numpy as np
from scipy.linalg import eigh

SPINS=tuple(itertools.product((0,1),repeat=3))
I8=np.eye(8,dtype=complex)
PHI=2*np.pi/3
TAU1=np.pi
TAU2=3*np.pi/4

def permutation(order):
    p=np.zeros((8,8),complex)
    for j,s in enumerate(SPINS):p[SPINS.index(tuple(s[i] for i in order)),j]=1
    return p

C=permutation((1,2,0))
CHI=(C-C.conj().T)/(4j)
R=permutation((2,1,0))
D=permutation((1,0,2))-permutation((0,2,1))
PPLUS=(I8+np.exp(-1j*PHI)*C+np.exp(-2j*PHI)*(C@C))/3
PMINUS=PPLUS.conj()
PSYM=I8-PPLUS-PMINUS
AA=(np.cos(np.sqrt(3)*np.pi)-1)/2
RR=float(AA**2)

def trace_real(a):return float(np.trace(a).real)

def ideal(tau=TAU1,phi=-PHI):
    def amp(f):return (np.cos(2*tau*np.cos(f/4))+np.cos(2*tau*np.sin(f/4)))/2
    return sum(amp(phi+t)*p for t,p in ((0.,PSYM),(PHI,PPLUS),(-PHI,PMINUS)))

@lru_cache(maxsize=32)
def model(phi=-PHI,bonds=(1.,1.,1.,1.),potential=(0.,0.,0.,0.)):
    """Full 120-state fixed-three-boson model. Return T, interaction, HC, input."""
    if len(bonds)!=4 or len(potential)!=4:raise ValueError('Four bonds and potentials required.')
    if not all(np.isfinite(x) for x in (*bonds,*potential,phi)):raise ValueError('Finite parameters required.')
    # Enumerate the three-boson occupation basis using multisets.
    modes=list(itertools.combinations_with_replacement(range(8),3))
    occ=[]
    for seq in modes:
        n=[0]*8
        for i in seq:n[i]+=1
        occ.append(tuple(n))
    lookup={n:i for i,n in enumerate(occ)}
    T=np.zeros((120,120),complex);interaction=np.zeros(120)
    for col,n in enumerate(occ):
        ns=[n[2*j]+n[2*j+1] for j in range(4)]
        interaction[col]=sum(v*(v-1)//2 for v in ns)
        T[col,col]=sum(potential[j]*ns[j] for j in range(4))
        for j in range(4):
            k=(j+1)%4
            for s in range(2):
                for src,dst,phase in ((2*k+s,2*j+s,np.exp(1j*phi/4)),
                                      (2*j+s,2*k+s,np.exp(-1j*phi/4))):
                    if n[src]:
                        new=list(n);new[src]-=1;new[dst]+=1
                        T[lookup[tuple(new)],col]+=-bonds[j]*phase*math.sqrt(n[src]*(n[dst]+1))
    hc=np.flatnonzero(interaction==0)
    ii=[]
    for s in SPINS:
        n=[0]*8
        for j,spin in enumerate(s,1):n[2*j+spin]=1
        ii.append(lookup[tuple(n)])
    return T,interaction,hc,np.array(ii),tuple(occ)

@lru_cache(maxsize=128)
def spectrum(u,phi=-PHI):
    if not np.isfinite(u) or u<0:raise ValueError('Use finite u>=0.')
    T,N,hc,ii,occ=model(phi)
    ev,V=eigh(T+np.diag(u*N))
    return ev,V

def return_operator(u,tau=TAU1,phi=-PHI):
    if not np.isfinite(tau) or tau<0:raise ValueError('Use finite tau>=0.')
    ev,V=spectrum(float(u),float(phi));ii=model(phi)[3]
    W=V[ii,:]
    return (W*np.exp(-1j*tau*ev))@W.conj().T

def contrast_coefficients(u,tau=TAU1,phi=PHI):
    kp=return_operator(u,tau,-phi);km=return_operator(u,tau,phi)
    delta=kp.conj().T@kp-km.conj().T@km
    g=float(np.vdot(CHI,delta).real/np.vdot(CHI,CHI).real)
    h=float(np.vdot(D,delta).real/np.vdot(D,D).real)
    residual=float(np.linalg.norm(delta-g*CHI-h*D,2))
    return g,h,residual,delta

def calibration(u,times=(TAU1,TAU2),singularity_guard=1e-10):
    vals=[contrast_coefficients(u,float(t)) for t in times]
    M=np.array([v[:2] for v in vals]);det=float(np.linalg.det(M))
    if abs(det)<singularity_guard:raise ArithmeticError('Numerically unresolved/degenerate two-time calibration.')
    weights=np.array([M[1,1],-M[0,1]])/det
    identity=weights[0]*vals[0][3]+weights[1]*vals[1][3]
    return {'u':float(u),'times':list(times),'coefficients':M.tolist(),'determinant':det,
            'weights':weights.tolist(),'operator_residual':float(np.linalg.norm(identity-CHI,2)),
            'variance_numerator':float(sum(abs(weights))**2),
            'scope':'Coefficients from complete finite matrices, not interval-certified values or device calibration.'}

def norm_bound(u,tau=TAU1):
    if not np.isfinite(u) or u<=0:raise ValueError('Use finite u>0.')
    if not np.isfinite(tau) or tau<0:raise ValueError('Use finite tau>=0.')
    return min(2.,20*np.sqrt(3)*tau/u)

def repeated(u,m,tau=TAU1):
    if not isinstance(m,int) or m<1:raise ValueError('Use integer m>=1.')
    k=return_operator(u,tau);km=np.linalg.matrix_power(k,m)
    out=km@km.conj().T/8;p=trace_real(out)
    if p<=np.finfo(float).tiny:raise ArithmeticError('Success probability underflow.')
    rho=out/p
    return {'u':float(u),'rounds':m,'total_success':p,'target_weight':trace_real(PPLUS@rho),
            'chirality':trace_real(CHI@rho),'return_operator_norm':float(np.linalg.norm(k,2)),
            'spectral_radius':float(max(abs(np.linalg.eigvals(k))))}


def first_correction(tau=TAU1,phi=-PHI):
    """K=K_HC+K1/u+O_tau(u^-2) via the single-doublon virtual paths.

    This finite sum computes the coefficient, not a certified remainder constant.
    """
    T,N,hc,ii,_=model(phi);q=np.flatnonzero(N>0)
    A=T[np.ix_(hc,hc)];B=T[np.ix_(q,hc)]
    S=B.conj().T@((1/N[q])[:,None]*B)
    e,V=eigh(A);s=V.conj().T@S@V
    de=e[:,None]-e[None,:]
    integ=tau*np.exp(1j*de*tau/2)*np.sinc(de*tau/(2*np.pi))
    k1=(V@(1j*np.exp(-1j*tau*e)[:,None]*s*integ))@V.conj().T
    loc=[list(hc).index(i) for i in ii]
    return k1[np.ix_(loc,loc)]
