#!/usr/bin/env python3
"""Reproducibility checks for a vacancy-mediated scalar-chirality filter.

All matrices are finite and built here; no earlier research checker is imported.
The 32-state hard-core Hamiltonian and 120-state finite-U Bose-Hubbard control
are independently constructed from occupation bases. Numeric tests supplement
rather than replace the derivation in PILOT.md.
"""
from __future__ import annotations
import argparse
from functools import lru_cache
import itertools
import json
import math
from pathlib import Path
import sys
import unittest
import numpy as np
from scipy.linalg import expm, eigh
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
REPORT: dict = {}
PI = np.pi
PHI = 2*PI/3
A = (np.cos(np.sqrt(3)*PI)-1)/2
R = float(A*A)
SPINS = list(itertools.product((1,2), repeat=3))
BASIS = [s for s in itertools.product((0,1,2), repeat=4) if s.count(0)==1]
INDEX = {s:i for i,s in enumerate(BASIS)}
START = [INDEX[(0,)+s] for s in SPINS]
I8 = np.eye(8,dtype=complex)
C = np.zeros((8,8),complex)
for j,s in enumerate(SPINS):
    C[SPINS.index((s[1],s[2],s[0])),j]=1
PAULI = [np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1.])]
def kron3(x,y,z): return np.kron(np.kron(x,y),z)
CHI = np.zeros((8,8),complex)
for perm in itertools.permutations(range(3)):
    parity=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
    CHI += parity*kron3(*(PAULI[j]/2 for j in perm))

def projector(theta: float) -> np.ndarray:
    return (I8+np.exp(-1j*theta)*C+np.exp(-2j*theta)*(C@C))/3
PROJECTORS = [projector(t) for t in (0.,PHI,-PHI)]

def hard_hamiltonian(phi: float, bonds=(1.,1.,1.,1.), potential=(0.,0.,0.,0.)) -> np.ndarray:
    """Positive phi is the phase accumulated by a clockwise vacancy."""
    if len(bonds)!=4 or len(potential)!=4: raise ValueError('Four bonds and potentials are required.')
    h=np.zeros((32,32),complex)
    for col,s in enumerate(BASIS):
        j=s.index(0);k=(j+1)%4
        h[col,col]=potential[j]
        target=list(s);target[j],target[k]=target[k],target[j]
        row=INDEX[tuple(target)]
        element=-bonds[j]*np.exp(1j*phi/4)
        h[row,col]+=element;h[col,row]+=element.conjugate()
    return h

def return_block(h: np.ndarray, t: float, indices=START) -> np.ndarray:
    return expm(-1j*t*h)[np.ix_(indices,indices)]

def charge_amplitude(flux: float,t: float,bonds=(1.,1.,1.,1.),potential=(0.,0.,0.,0.)) -> complex:
    h=np.diag(np.asarray(potential,complex))
    for j in range(4):
        k=(j+1)%4
        element=-bonds[j]*np.exp(1j*flux/4)
        h[k,j]+=element;h[j,k]+=element.conjugate()
    return complex(expm(-1j*t*h)[0,0])

def uniform_amplitude(flux: float,t: float) -> float:
    return float((np.cos(2*t*np.cos(flux/4))+np.cos(2*t*np.sin(flux/4)))/2)

def partial_trace(rho: np.ndarray,keep: tuple[int,...]) -> np.ndarray:
    a=rho.reshape((2,)*6);n=3
    for subsystem in reversed(range(3)):
        if subsystem not in keep:
            a=np.trace(a,axis1=subsystem,axis2=subsystem+n);n-=1
    return a.reshape((2**len(keep),)*2)

def partial_transpose_pair(rho: np.ndarray) -> np.ndarray:
    return rho.reshape(2,2,2,2).transpose(0,3,2,1).reshape(4,4)

def compositions(n: int,m: int):
    if m==1:
        yield (n,)
    else:
        for k in range(n+1):
            for tail in compositions(n-k,m-1):yield (k,)+tail

@lru_cache(None)
def finite_u_data(phi: float):
    states=list(compositions(3,8));lookup={s:i for i,s in enumerate(states)}
    hopping=np.zeros((120,120),complex);interaction=np.zeros(120)
    for col,s in enumerate(states):
        interaction[col]=sum((s[2*j]+s[2*j+1])*(s[2*j]+s[2*j+1]-1)/2 for j in range(4))
        for j in range(4):
            k=(j+1)%4
            for spin in range(2):
                src=2*k+spin;dst=2*j+spin
                if not s[src]:continue
                occupation=list(s)
                amp=math.sqrt(occupation[src]*(occupation[dst]+1))
                occupation[src]-=1;occupation[dst]+=1
                row=lookup[tuple(occupation)]
                element=-np.exp(1j*phi/4)*amp
                hopping[row,col]+=element;hopping[col,row]+=element.conjugate()
    hc_indices=[]
    for s in BASIS:
        n=[0]*8
        for j,spin in enumerate(s):
            if spin:n[2*j+spin-1]=1
        hc_indices.append(lookup[tuple(n)])
    return hopping,interaction,hc_indices,[hc_indices[i] for i in START]

class Checks(unittest.TestCase):
    def test_01_spin_algebra_and_witness(self):
        self.assertLess(np.linalg.norm(CHI-(C-C.conj().T)/(4j)),1e-14)
        self.assertLess(np.linalg.norm(C@C@C-I8),1e-14)
        self.assertLess(np.linalg.norm(sum(PROJECTORS)-I8),1e-14)
        for P,rank in zip(PROJECTORS,(4,2,2)):
            self.assertLess(np.linalg.norm(P@P-P),2e-14)
            self.assertAlmostEqual(np.trace(P).real,rank,places=13)
        cross=np.kron(PAULI[0]/2,PAULI[1]/2)-np.kron(PAULI[1]/2,PAULI[0]/2)
        self.assertAlmostEqual(np.linalg.norm(cross,2)/2,.25,places=14)
        v=np.kron(np.array([1,0]),np.array([0,1,1j,0])/np.sqrt(2))
        self.assertAlmostEqual(abs(np.vdot(v,CHI@v)),.25,places=14)
        rng=np.random.default_rng(161)
        peak=0.
        for _ in range(150):
            a=rng.normal(size=2)+1j*rng.normal(size=2);a/=np.linalg.norm(a)
            b=rng.normal(size=4)+1j*rng.normal(size=4);b/=np.linalg.norm(b)
            v=np.kron(a,b);value=abs(np.vdot(v,CHI@v))
            self.assertLessEqual(value,.25+1e-13);peak=max(peak,float(value))
        REPORT['spin_algebra']={'chirality_eigenvalues':np.linalg.eigvalsh(CHI).tolist(),
            'projector_ranks':[4,2,2],'biseparable_bound':.25,'random_product_bipartition_maximum':peak,
            'scope':'The bound follows analytically from rotational invariance and a two-spin norm; samples are only checks.'}

    def test_02_full_hamiltonian_instrument(self):
        errors=[]
        for phi in (0.,.43,PHI,-PHI,2.73):
            h=hard_hamiltonian(phi)
            self.assertLess(np.linalg.norm(h-h.conj().T),1e-14)
            for t in (.31,1.17,PI,2.4):
                actual=return_block(h,t)
                expected=sum(uniform_amplitude(phi+theta,t)*P for theta,P in zip((0.,PHI,-PHI),PROJECTORS))
                errors.append(float(np.linalg.norm(actual-expected,2)))
                self.assertLess(errors[-1],1e-12)
        kp=return_block(hard_hamiltonian(-PHI),PI)
        km=return_block(hard_hamiltonian(PHI),PI)
        Pplus=PROJECTORS[1]
        self.assertLess(np.linalg.norm(kp-(Pplus+A*(I8-Pplus))),1e-12)
        difference=kp.conj().T@kp-km.conj().T@km
        self.assertLess(np.linalg.norm(difference-4*(1-R)/np.sqrt(3)*CHI),1e-12)
        # Independent time integrator for all eight inputs, using sparse exp action.
        h=hard_hamiltonian(-PHI);inputs=np.eye(32)[:,START]
        direct=expm_multiply(csr_matrix(-1j*PI*h),inputs)
        error=float(np.linalg.norm(direct[START]-kp))
        self.assertLess(error,2e-12)
        REPORT['instrument']={'dimension':32,'tested_flux_time_pairs':len(errors),
            'maximum_operator_error':max(errors),'sparse_exponential_error':error,
            'return_amplitude_nonselected':float(A),'false_return_probability':R,
            'chirality_readout_coefficient':float(4*(1-R)/np.sqrt(3))}

    def test_03_heralding_repetition_and_local_marginals(self):
        k=return_block(hard_hamiltonian(-PHI),PI);rows=[]
        for m in (1,2,3,4):
            km=np.linalg.matrix_power(k,m);out=km@km.conj().T/8;p=float(np.trace(out).real);rho=out/p
            rm=R**m;expected_p=(1+3*rm)/4
            expected=(PROJECTORS[1]+rm*(I8-PROJECTORS[1]))/(2+6*rm)
            self.assertAlmostEqual(p,expected_p,places=13)
            self.assertLess(np.linalg.norm(rho-expected),2e-12)
            chirality=float(np.trace(CHI@rho).real)
            self.assertGreater(chirality,.25)
            self.assertAlmostEqual(chirality,np.sqrt(3)/4*(1-rm)/(1+3*rm),places=13)
            minima=[]
            for keep in ((0,1),(0,2),(1,2)):
                pair=partial_trace(rho,keep)
                mineig=float(np.linalg.eigvalsh(partial_transpose_pair(pair))[0]);minima.append(mineig)
                self.assertGreater(mineig,-2e-12)
            for j in range(3):self.assertLess(np.linalg.norm(partial_trace(rho,(j,))-np.eye(2)/2),1e-12)
            rows.append({'returns':m,'dimensionless_time':float(m*PI),'success_probability':p,
                'target_sector_weight':float(np.trace(PROJECTORS[1]@rho).real),'chirality':chirality,
                'pairwise_PT_min_eigenvalues':minima})
        # A polarized input has no target chiral weight: filtering cannot create it.
        polarized=np.zeros(8);polarized[0]=1
        out=k@polarized
        self.assertLess(np.linalg.norm(out-A*polarized),1e-12)
        REPORT['heralding']={'input':'I_8/8 on the three occupied spins; initial vacancy is known',
            'rows':rows,'postselection':'Every nonreturn terminates that trial; probabilities include failed trials.',
            'measurement_assumption':'Repeated filters require nondestructive spin-insensitive vacancy-location measurement.'}

    def test_04_scalar_disorder_and_open_chain_control(self):
        rng=np.random.default_rng(716);rows=[]
        for _ in range(12):
            bonds=tuple(rng.uniform(.5,1.5,4));onsite=tuple(rng.uniform(-.7,.7,4));t=float(rng.uniform(.2,4))
            phi=float(rng.uniform(.3,2.7))
            hplus=hard_hamiltonian(phi,bonds,onsite);hminus=hard_hamiltonian(-phi,bonds,onsite)
            kp=return_block(hminus,t);km=return_block(hplus,t)
            expected=sum(charge_amplitude(phi+theta,t,bonds,onsite)*P for theta,P in zip((0.,PHI,-PHI),PROJECTORS))
            reduction_error=float(np.linalg.norm(km-expected,2));self.assertLess(reduction_error,2e-12)
            d=abs(charge_amplitude(PHI-phi,t,bonds,onsite))**2-abs(charge_amplitude(PHI+phi,t,bonds,onsite))**2
            readout_error=float(np.linalg.norm(kp.conj().T@kp-km.conj().T@km-4*d/np.sqrt(3)*CHI,2))
            self.assertLess(readout_error,2e-12)
            rows.append({'time':t,'flux':phi,'bonds':bonds,'hole_potentials':onsite,
                         'reduction_error':reduction_error,'readout_operator_error':readout_error,
                         'calibration_contrast':float(d)})
        for phi in (0.,PHI,-PHI):
            k=return_block(hard_hamiltonian(phi,bonds=(1,1,1,0)),PI)
            self.assertLess(np.linalg.norm(k-np.trace(k)/8*I8),2e-12)
        # Static real reciprocal ring at zero flux cannot distinguish handedness.
        for t in (.23,1.51,PI):
            k=return_block(hard_hamiltonian(0),t);e=k.conj().T@k
            self.assertAlmostEqual(np.trace(PROJECTORS[1]@e).real,np.trace(PROJECTORS[2]@e).real,places=12)
        REPORT['scalar_disorder']={'tests':rows,'open_chain':'Return operator is spin scalar; no chiral filter',
            'scope':'Exact operator structure survives static spin-independent disorder. Contrast is not protected and may vanish.'}

    def test_05_finite_repulsion_control(self):
        hopping,interaction,hc,initial=finite_u_data(PHI)
        self.assertEqual(len(hc),32);self.assertEqual(len(initial),8)
        self.assertLess(np.linalg.norm(hopping-hopping.conj().T),1e-13)
        self.assertLess(np.linalg.norm(hopping[np.ix_(hc,hc)]-hard_hamiltonian(PHI)),1e-13)
        self.assertEqual(np.count_nonzero(interaction==0),32)
        rows=[]
        for u in (8.,16.,32.,64.,128.):
            h=hopping+np.diag(u*interaction);energies,vectors=eigh(h)
            k=((vectors*np.exp(-1j*PI*energies))@vectors.conj().T)[np.ix_(initial,initial)]
            out=k@k.conj().T/8;p=float(np.trace(out).real);rho=out/p
            chi=float(-np.trace(CHI@rho).real)
            # The opposite-flux return block is K^T. The POVM difference need not
            # be purely chirality once particles can pass through double occupancy.
            kr=k.T;diff=kr.conj().T@kr-k.conj().T@k
            coefficient=float(np.trace(CHI@diff).real/np.trace(CHI@CHI).real)
            remainder=float(np.linalg.norm(diff-coefficient*CHI,2))
            self.assertGreater(p,0);self.assertLess(p,1)
            self.assertGreater(np.linalg.eigvalsh(rho)[0],-1e-12)
            rows.append({'U_over_J':u,'success_probability':p,'selected_signed_chirality':chi,
                'GME_certified_by_ideal_spin_witness':bool(chi>.25),'chi_component_in_flux_contrast':coefficient,
                'nonchirality_operator_norm':remainder,
                'return_filter_C_commutator':float(np.linalg.norm(k@C-C@k,2))})
            if u==16.:
                initial_vectors=np.eye(120)[:,initial]
                check=expm_multiply(csr_matrix(-1j*PI*h),initial_vectors)[initial]
                independent_error=float(np.linalg.norm(check-k,2));self.assertLess(independent_error,5e-11)
        self.assertGreater(rows[1]['selected_signed_chirality'],.25)
        self.assertGreater(rows[1]['nonchirality_operator_norm'],.01)
        REPORT['finite_repulsion']={'dimension':120,'independent_exponential_error_U16':independent_error,'rows':rows,
           'scope':'Finite-system controls, not a uniform finite-U theorem or an apparatus model. Accepted occupancy is (0,1,1,1); a doublon is a failure. Hard-core calibration is not exact at finite U.'}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('Refusing to overwrite an existing report.')
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Checks))
    REPORT.update(status='PASS' if result.wasSuccessful() else 'FAIL',test_groups=result.testsRun,
                  random_seeds_fixed=True,prior_project_modules_imported=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:json.dump(REPORT,stream,indent=2,sort_keys=True,allow_nan=False);stream.write('\n')
    raise SystemExit(0 if result.wasSuccessful() else 1)
if __name__=='__main__':main()
