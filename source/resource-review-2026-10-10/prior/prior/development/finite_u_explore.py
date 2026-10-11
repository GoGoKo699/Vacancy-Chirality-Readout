import itertools, numpy as np
from scipy.linalg import eigh
from math import pi
import explore

def compositions(n,m):
 if m==1:
  yield (n,)
 else:
  for k in range(n+1):
   for rest in compositions(n-k,m-1):yield (k,)+rest
states=list(compositions(3,8)); lookup={s:i for i,s in enumerate(states)}
inds=[]
for spins in explore.spins:
 occ=[0]*8
 for j,s in enumerate(spins,1):occ[2*j+s-1]=1
 inds.append(lookup[tuple(occ)])
base=np.zeros((120,120),complex);diag=np.zeros(120)
for i,s in enumerate(states):
 diag[i]=sum((s[2*j]+s[2*j+1])*(s[2*j]+s[2*j+1]-1)/2 for j in range(4))
 for j in range(4):
  k=(j+1)%4
  for sp in range(2):
   source=2*k+sp;dest=2*j+sp
   if s[source]:
    occ=list(s);amp=np.sqrt(occ[source]*(occ[dest]+1));occ[source]-=1;occ[dest]+=1
    row=lookup[tuple(occ)];val=-np.exp(1j*2*pi/3/4)*amp
    base[row,i]+=val;base[i,row]+=val.conjugate()
for U in (8,16,32,64,128,256):
 e,v=eigh(base+U*np.diag(diag));K=((v*np.exp(-1j*pi*e))@v.conj().T)[np.ix_(inds,inds)]
 out=K@K.conj().T/8;p=np.trace(out).real;rho=out/p
 print('finiteU',U,'p',p,'chi',-np.trace(explore.chi@rho).real,'instrument distance',np.linalg.norm(K-(explore.a*explore.I+(1-explore.a)*(explore.I+np.exp(2j*pi/3)*explore.C+np.exp(4j*pi/3)*explore.C@explore.C)/3),2))
