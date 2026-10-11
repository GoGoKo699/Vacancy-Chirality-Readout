from pathlib import Path
import sys,itertools,numpy as np
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'prior'))
import finite_repulsion as f
P12=f.permutation((1,0,2));P23=f.permutation((0,2,1));P13=f.R
basis=[f.I8,P12,P23,P13,f.CHI]
B=np.column_stack([np.r_[b.real.ravel(),b.imag.ravel()] for b in basis])
for u in (0.,8.,16.,128.):
 T,N,hc,ii,occ=f.model();ev,W=f.spectrum(u,-f.PHI)
 for tau in (np.pi,3*np.pi/4):
  Ucol=(W*np.exp(-1j*tau*ev))@W[ii,:].conj().T
  groups={}
  for a,n in enumerate(occ):
   parity=tuple((n[2*j]+n[2*j+1])%2 for j in range(4));groups.setdefault(parity,[]).append(a)
  eff=[Ucol[inds,:].conj().T@Ucol[inds,:] for inds in groups.values()]
  E=np.column_stack([np.r_[e.real.ravel(),e.imag.ravel()] for e in eff])
  target=np.r_[f.CHI.real.ravel(),f.CHI.imag.ravel()]
  vals=np.linalg.svd(E,compute_uv=False)
  weights=np.linalg.lstsq(E,target,rcond=1e-10)[0]
  print('u,tau',u,tau,'rank',sum(vals>1e-9),'residual',np.linalg.norm(E@weights-target),'weights',list(zip(groups,weights)))
  probs=np.array([np.trace(e).real/8 for e in eff]);print('var maximixed',np.sum(probs*weights**2),'range',np.ptp(weights))
