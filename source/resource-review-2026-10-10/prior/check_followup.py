#!/usr/bin/env python3
"""Five new finite-repulsion checks. Outputs are new files outside prior/.

No archived scientific checker is imported. Full first-quantized assembly,
Taylor time-evolution, spectator inputs, and symmetry controls supplement proofs.
"""
from __future__ import annotations
import argparse,itertools,json,math,sys
from pathlib import Path
import unittest
sys.dont_write_bytecode=True
import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply
import finite_repulsion as f
ROOT=Path(__file__).resolve().parent
REPORT={}

def taylor_unitary(H,tau):
    """Independent scaling/squaring Taylor polynomial and exact-arithmetic tail.

    The tail does not include floating-point roundoff, reported by comparison.
    """
    scaled=-1j*tau*H
    s=max(0,int(np.ceil(np.log2(max(1.,2*np.linalg.norm(scaled,2))))))
    X=scaled/(2**s);term=np.eye(len(H),dtype=complex);out=term.copy()
    for j in range(1,23):term=term@X/j;out+=term
    x=np.linalg.norm(X,2);tail=np.exp(x)*x**23/math.factorial(23)
    # Exact exponential unitary; polynomial powering error controlled by telescoping.
    error=(2**s)*tail*(1+tail)**(2**s-1)
    for _ in range(s):out=out@out
    return out,float(error)

class Checks(unittest.TestCase):
    def test_01_model_and_uniform_bound(self):
        T,N,hc,ii,occ=f.model();n=8
        one=np.zeros((8,8),complex)
        for j in range(4):
            k=(j+1)%4
            for spin in range(2):
                one[2*j+spin,2*k+spin]=-np.exp(-1j*f.PHI/4)
                one[2*k+spin,2*j+spin]=-np.exp(1j*f.PHI/4)
        eye=np.eye(8);dist=np.kron(np.kron(one,eye),eye)+np.kron(np.kron(eye,one),eye)+np.kron(np.kron(eye,eye),one)
        W=np.zeros((512,120),complex);interaction=np.zeros(512)
        for xyz in itertools.product(range(8),repeat=3):
            idx=np.ravel_multi_index(xyz,(8,8,8));sites=[v//2 for v in xyz]
            interaction[idx]=sum(sites[a]==sites[b] for a,b in ((0,1),(0,2),(1,2)))
        for col,counts in enumerate(occ):
            seq=tuple(i for i,count in enumerate(counts) for _ in range(count))
            perms=set(itertools.permutations(seq));a=1/np.sqrt(len(perms))
            for seq in perms:W[np.ravel_multi_index(seq,(8,8,8)),col]=a
        kinetic_error=float(np.linalg.norm(W.conj().T@dist@W-T,2))
        int_error=float(np.linalg.norm(W.conj().T@(interaction[:,None]*W)-np.diag(N),2))
        self.assertLess(kinetic_error,3e-14);self.assertLess(int_error,3e-14)
        # Exact integer Gram row count underlying ||QTP||^2<=12 J^2.
        T0,N0,p0,_,_=f.model(0.);q0=np.flatnonzero(N0>0);B=T0[np.ix_(q0,p0)]
        squares=np.rint(abs(B)**2).astype(int)
        self.assertLess(np.max(abs(abs(B)**2-squares)),1e-13)
        gram=np.zeros((32,32),dtype=int)
        for r in range(32):
            for c in range(32):
                for x,y in zip(squares[:,r],squares[:,c]):
                    prod=int(x*y);s=math.isqrt(prod);self.assertEqual(s*s,prod);gram[r,c]+=s
        self.assertTrue(np.all(gram.sum(axis=1)==12))
        rows=[]
        for u in (16.,64.,256.):
            H=T+np.diag(u*N)
            for tau in (.31,1.,np.pi,4.):
                U=expm(-1j*tau*H);q=np.flatnonzero(N>0)
                leak=float(np.linalg.norm(U[np.ix_(q,hc)],2));err=float(np.linalg.norm(U[np.ix_(ii,ii)]-f.ideal(tau),2))
                self.assertLessEqual(leak,min(1.,10/u)+1e-12)
                self.assertLessEqual(err,f.norm_bound(u,tau)+1e-12)
                rows.append({'u':u,'tau':tau,'leakage_norm':leak,'return_error':err,'bound':f.norm_bound(u,tau)})
        REPORT['uniform_bound']={'first_quantized_dimension':512,'bosonic_dimension':120,'kinetic_error':kinetic_error,'interaction_error':int_error,'exact_gram_row_sum':12,'cases':rows}

    def test_02_reference_input_and_repetitions(self):
        rng=np.random.default_rng(9217);rows=[]
        for u in (32.,128.,2048.):
            K=f.return_operator(u);K0=f.ideal();eta=float(np.linalg.norm(K-K0,2))
            for m in (1,2,4):
                X=np.linalg.matrix_power(K,m);Y=np.linalg.matrix_power(K0,m)
                eps=float(np.linalg.norm(X-Y,2));self.assertLessEqual(eps,m*eta+1e-12)
                for _ in range(4):
                    v=rng.normal(size=24)+1j*rng.normal(size=24);v/=np.linalg.norm(v)
                    x=np.kron(X,np.eye(3))@v;y=np.kron(Y,np.eye(3))@v
                    p=np.vdot(x,x).real;p0=np.vdot(y,y).real
                    D=np.sqrt(max(0.,1-abs(np.vdot(x,y))**2/(p*p0)))
                    self.assertLessEqual(abs(np.sqrt(p)-np.sqrt(p0)),eps+1e-12)
                    self.assertLessEqual(D,min(1.,eps/np.sqrt(p0))+1e-10)
                rows.append(f.repeated(u,m))
        longrows=[f.repeated(16.,m) for m in (1,4,16,64)]
        self.assertLess(longrows[-1]['total_success'],1e-6)
        self.assertGreater(longrows[1]['chirality'],.43)
        self.assertLess(longrows[1]['total_success'],.12)
        p0=(1+3*f.RR)/4;chi0=np.sqrt(3)/4*(1-f.RR)/(1+3*f.RR)
        eps=f.norm_bound(2048.);pmin=(np.sqrt(p0)-eps)**2;witness=chi0-np.sqrt(3)/2*eps/np.sqrt(p0)
        self.assertGreater(witness,.30);self.assertGreater(pmin,.218)
        REPORT['instrument_and_repetition']={'reference_dimension':3,'cases':rows,'fixed_u16_long_sequence':longrows,'analytical_sufficient_u2048':{'epsilon':eps,'success_lower_bound':pmin,'chirality_lower_bound':witness},'scope':'Reference inputs test the proven uniform statements; long sequences are finite-matrix diagnostics, not an all-U spectral-radius theorem.'}

    def test_03_symmetry_and_leading_correction(self):
        rows=[]
        for u in (0.,8.,16.,64.,128.):
            for tau in (.57,3*np.pi/4,np.pi,4.1):
                kp=f.return_operator(u,tau);km=f.return_operator(u,tau,f.PHI)
                self.assertLess(np.linalg.norm(km-f.R@kp@f.R,2),5e-12)
                self.assertLess(np.linalg.norm(km-kp.T,2),5e-12)
                g,h,res,_=f.contrast_coefficients(u,tau)
                self.assertLess(res,5e-12)
                rows.append({'u':u,'tau':tau,'g':g,'h':h,'residual':res})
        # Leading effective Hamiltonian correction, compared after cancelling K0.
        corrected=[]
        for tau in (3*np.pi/4,np.pi):
            K0=f.ideal(tau);K1=f.first_correction(tau)
            df=(K0.conj().T@K1+K1.conj().T@K0)
            df-=f.R@df@f.R
            g1=float(np.vdot(f.CHI,df).real/np.vdot(f.CHI,f.CHI).real)
            h1=float(np.vdot(f.D,df).real/np.vdot(f.D,f.D).real)
            self.assertLess(abs(g1),1e-10)
            vals=[]
            for u in (128.,256.,512.,1024.):
                K=f.return_operator(u,tau);scaled=np.linalg.norm(K-K0-K1/u,2)*u*u
                self.assertLess(scaled,180.)
                vals.append({'u':u,'u_squared_remainder':float(scaled)})
            corrected.append({'tau':tau,'g1':g1,'h1':h1,'cases':vals})
        # Breaking reflection while preserving spin independence defeats two-operator form.
        bonds=(1.,.7,1.2,.9)
        def dis(phi):
            T,N,_,ii,_=f.model(phi,bonds);U=expm(-1j*np.pi*(T+np.diag(16*N)));return U[np.ix_(ii,ii)]
        kp,km=dis(-f.PHI),dis(f.PHI);delta=kp.conj().T@kp-km.conj().T@km
        fit=sum(np.vdot(M,delta).real/np.vdot(M,M).real*M for M in (f.CHI,f.D))
        remainder=float(np.linalg.norm(delta-fit,2));self.assertGreater(remainder,.001)
        REPORT['finite_U_symmetry']={'cases':rows,'leading_virtual_correction':corrected,'reflection_broken_residual':remainder,'scope':'Two-observable finite-U calibration requires reflection plus spin-rotation symmetry; arbitrary hard-core disorder robustness does not transfer.'}

    def test_04_four_setting_readout(self):
        rows=[];rng=np.random.default_rng(9183)
        for u in (8.,16.,32.,64.,128.):
            cal=f.calibration(u);w=cal['weights'];self.assertLess(cal['operator_residual'],1e-10)
            errors=[]
            for _ in range(12):
                X=rng.normal(size=(8,8))+1j*rng.normal(size=(8,8));rho=X@X.conj().T;rho/=np.trace(rho)
                contrasts=[]
                for tau in (f.TAU1,f.TAU2):
                    kp=f.return_operator(u,tau);km=f.return_operator(u,tau,f.PHI)
                    contrasts.append(f.trace_real((kp.conj().T@kp-km.conj().T@km)@rho))
                errors.append(abs(np.dot(w,contrasts)-f.trace_real(f.CHI@rho)))
            self.assertLess(max(errors),5e-12);cal['random_state_error']=max(errors);rows.append(cal)
        c16=rows[1];self.assertGreater(abs(c16['determinant']),.04)
        self.assertLess(c16['variance_numerator'],1.04)
        # Zero flux cannot furnish either odd observable or an invertible calibration.
        g,h,res,_=f.contrast_coefficients(16.,np.pi,0.)
        self.assertLess(abs(g)+abs(h),1e-12)
        try:f.calibration(16.,times=(np.pi,np.pi))
        except ArithmeticError:pass
        else:self.fail('Identical times must reject singular calibration.')
        REPORT['calibration']={'rows':rows,'scope':'Exact identity conditional on nonzero determinant and calibrated coefficients. Numerical examples are not interval or device calibration certificates.'}

    def test_05_independent_evolution_and_noise_accounting(self):
        u=16.;T,N,_,ii,_=f.model();H=T+np.diag(u*N);reports=[]
        for tau in (f.TAU1,f.TAU2):
            V,tail=taylor_unitary(H,tau);actual=f.return_operator(u,tau)
            err=float(np.linalg.norm(V[np.ix_(ii,ii)]-actual,2));self.assertLess(err,2e-11)
            sparse=expm_multiply(csr_matrix(-1j*tau*H),np.eye(120)[:,ii])[ii]
            serr=float(np.linalg.norm(sparse-actual,2));self.assertLess(serr,2e-11)
            reports.append({'tau':tau,'Taylor_vs_spectral_norm':err,'sparse_vs_spectral_norm':serr,'exact_arithmetic_Taylor_tail':tail})
        c=f.calibration(16.);w=np.array(c['weights']);budgets=np.array([.002,.001]);bound=float(np.dot(abs(w),budgets))
        # Adversarial operator calibration errors obey the same triangle bound.
        X=.002*f.D/np.linalg.norm(f.D,2);Y=-.001*f.CHI/np.linalg.norm(f.CHI,2)
        self.assertLessEqual(np.linalg.norm(w[0]*X+w[1]*Y,2),bound+1e-14)
        coeff=np.repeat(abs(w),2);alloc=coeff/coeff.sum()
        var_bound=float(np.sum(coeff**2/(4*alloc)))
        self.assertAlmostEqual(var_bound,sum(abs(w))**2,places=12)
        REPORT['independent_checks']={'full_sector_dimension':120,'evolutions':reports,'calibration_error_bound_example':bound,'optimized_shot_variance_numerator':var_bound,'scope':'Taylor tail is an exact-arithmetic truncation bound; comparisons include floating-point roundoff and are not interval-certified.'}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    out=args.output.resolve()
    if out.exists() or (ROOT/'prior') in out.parents:p.error('Use a new output path outside preserved prior/.')
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Checks))
    REPORT.update(status='PASS' if result.wasSuccessful() else 'FAIL',test_groups=result.testsRun,archived_modules_imported=False)
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as s:json.dump(REPORT,s,indent=2,sort_keys=True,allow_nan=False);s.write('\n')
    raise SystemExit(0 if result.wasSuccessful() else 1)
if __name__=='__main__':main()
