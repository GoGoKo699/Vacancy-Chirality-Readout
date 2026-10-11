import numpy as np,itertools
from scipy.linalg import expm
from math import pi,sqrt
states=[s for s in itertools.product((0,1,2),repeat=4) if s.count(0)==1]
index={s:i for i,s in enumerate(states)}
spins=list(itertools.product((1,2),repeat=3))
start=[index[(0,)+s] for s in spins]
C=np.zeros((8,8),complex)
for i,s in enumerate(spins):C[spins.index((s[1],s[2],s[0])),i]=1
chi=(C-C.conj().T)/(4j)
I=np.eye(8)
a=(np.cos(sqrt(3)*pi)-1)/2
r=a*a
print('a,r,p,chirality',a,r,(1+3*r)/4,sqrt(3)/4*(1-r)/(1+3*r))
def H(phi):
 h=np.zeros((32,32),complex)
 for i,s in enumerate(states):
  j=s.index(0);k=(j+1)%4
  t=list(s);t[j],t[k]=t[k],t[j]
  m=index[tuple(t)]
  h[m,i] += -np.exp(1j*phi/4)
  h[i,m] += -np.exp(-1j*phi/4)
 return h
for phi in (2*pi/3,-2*pi/3):
 u=expm(-1j*pi*H(phi));K=u[np.ix_(start,start)]
 theta=-phi
 Proj=(I+np.exp(-1j*theta)*C+np.exp(-2j*theta)*C@C)/3
 expected=a*I+(1-a)*Proj
 print(phi,'Kerror',np.linalg.norm(K-expected),'eig chi',[np.linalg.eigvalsh(chi)],'K',np.linalg.eigvalsh(K))
 print('spinpostchi',np.trace(chi@K@K.conj().T).real/np.trace(K@K.conj().T).real)
for n,m in ((1,2),(15,26),(209,362)):
 rN=(np.cos(sqrt(3)*n*pi)-1)**2/4
 print('Pell',n,m,rN,rN*n**4)
