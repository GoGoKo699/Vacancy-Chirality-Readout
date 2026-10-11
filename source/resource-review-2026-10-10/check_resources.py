#!/usr/bin/env python3
"""Five bounded resource/readout check groups; no apparatus or interval certification."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import unittest
import numpy as np
sys.dont_write_bytecode=True
from readout import (ROOT, fr, occupations, site_occupancy, parity_acceptance, effect,
                     phase_damping, damped_filter, damped_recurrence, sector_hamiltonian, loss_evolution)
REPORT={}


class Checks(unittest.TestCase):
    def test_01_exact_parity_projector(self):
        counts=[]
        for n in range(4):
            occ=occupations(n)
            parity=[parity_acceptance(v) for v in occ]
            exact=[float(site_occupancy(v)==(0,1,1,1)) for v in occ]
            self.assertEqual(parity,exact)
            counts.append({'particles':n,'dimension':len(occ),'accepted':sum(parity)})
        self.assertEqual(counts[3]['accepted'],8)
        # The number bound cannot be omitted: a five-atom configuration aliases.
        alias=(0,0,3,0,1,0,1,0)
        self.assertEqual(parity_acceptance(alias),1.)
        self.assertNotEqual(site_occupancy(alias),(0,1,1,1))
        residuals=[]
        for tau in (np.pi,3*np.pi/4):
            for phi in (-fr.PHI,fr.PHI):
                K=fr.return_operator(16.,tau,phi)
                err=np.linalg.norm(effect(16.,tau,phi)-K.conj().T@K,2)
                self.assertLess(err,2e-13)
                residuals.append(float(err))
        REPORT['parity']={'exact_counts':counts,'effect_residuals':residuals,
                          'five_particle_alias':list(alias),'scope':'Ideal total on-site parity and known particle number at most three; no spin-selective parity assumption.'}

    def test_02_loss_is_heralded(self):
        H,occ=sector_hamiltonian(3,16.,-fr.PHI)
        T,W,_,_,oldocc=fr.model(-fr.PHI)
        self.assertEqual(occ,oldocc)
        self.assertLess(np.linalg.norm(H.toarray()-T-np.diag(16*W)),1e-13)
        rng=np.random.default_rng(731)
        x=rng.normal(size=(8,8))+1j*rng.normal(size=(8,8));rho=x@x.conj().T;rho/=np.trace(rho)
        rows=[]
        for tau,rate in ((np.pi,.04),(3*np.pi/4,.1)):
            result=loss_evolution(16.,tau,rate,rho)
            K=fr.return_operator(16.,tau)
            expected=np.exp(-3*rate*tau)*K@rho@K.conj().T
            err=float(np.linalg.norm(result['accepted']-expected))
            self.assertLess(err,2e-12)
            self.assertAlmostEqual(result['total_trace'],1.,places=11)
            self.assertAlmostEqual(result['parity_probability'],np.trace(expected).real,places=11)
            self.assertGreater(sum(result['sector_traces'][:3]),.2)
            rows.append({k:v for k,v in result.items() if k!='accepted'}|{'tau':tau,'rate_over_J':rate,'accepted_error':err})
        REPORT['loss']={'checks':rows,'scope':'Independent complete number-block Lindblad diagnostic. Equal local spin-independent loss, no gains, ideal final parity. Not a measured atom lifetime.'}

    def test_03_instrument_dephasing_and_repetition(self):
        K=fr.ideal();chi0=fr.trace_real(fr.CHI@K@K.conj().T)/fr.trace_real(K@K.conj().T)
        threshold=float(np.sqrt(.25/chi0))
        self.assertAlmostEqual(threshold,.8022166239295353,places=12)
        rho=K@K.conj().T/8
        for eta in (0.,.8,.9,1.):
            out=phase_damping(rho,eta)
            self.assertAlmostEqual(np.trace(out).real,np.trace(rho).real,places=13)
            self.assertAlmostEqual(fr.trace_real(fr.CHI@out),eta*eta*fr.trace_real(fr.CHI@rho),places=13)
        complete=phase_damping(rho,0.)
        self.assertLess(np.linalg.norm(complete-np.diag(np.diag(complete))),1e-14)
        self.assertGreaterEqual(np.min(np.diag(complete).real),0.)
        # Diagonal in a local product basis is explicitly fully separable.
        rows=[]
        for eta in (0.,.8,.9,.99,1.):
            for m in (1,4,16):
                direct=damped_filter(eta,m);rec=damped_recurrence(eta,m)
                self.assertLess(abs(direct['probability']-rec['probability']),2e-13)
                self.assertLess(abs(direct['chirality']-rec['chirality']),2e-13)
                if eta<1:self.assertLess(rec['spectral_radius'],1.)
                rows.append(rec)
        # General input traces, including coherences, agree for all such channels.
        rng=np.random.default_rng(19)
        for _ in range(5):
            psi=rng.normal(size=8)+1j*rng.normal(size=8);psi/=np.linalg.norm(psi)
            accepted=K@np.outer(psi,psi.conj())@K.conj().T
            self.assertAlmostEqual(np.trace(phase_damping(accepted,.7)).real,np.trace(accepted).real,places=13)
        REPORT['dephasing']={'one_round_ideal_chirality':chi0,'witness_eta_threshold':threshold,'recurrence':rows,
                            'scope':'Explicit phase-damping diagnostic after each accepted outcome. Same outcome effects do not imply preserved output entanglement.'}

    def test_04_matched_sampling_costs(self):
        X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1.,-1.])
        kron=lambda a,b,c:np.kron(np.kron(a,b),c)
        expansion=(kron(X,Y,Z)+kron(Y,Z,X)+kron(Z,X,Y)-kron(X,Z,Y)-kron(Z,Y,X)-kron(Y,X,Z))/8
        self.assertLess(np.linalg.norm(expansion-fr.CHI),1e-14)
        self.assertLess(np.linalg.norm((fr.C-fr.C.conj().T)/(2j)-2*fr.CHI),1e-14)
        calibration=fr.calibration(16.)
        w=np.array(calibration['weights']);W=sum(abs(w))
        probs=[];variance=0.
        for tau,weight in zip(calibration['times'],w):
            p=[]
            for phi in (-fr.PHI,fr.PHI):
                E=effect(16.,tau,phi);p.append(float(np.trace(E).real/8))
                allocation=abs(weight)/(2*W)
                variance+=weight*weight*p[-1]*(1-p[-1])/allocation
            probs.append(p)
        self.assertAlmostEqual(variance,.8759834952660344,places=10)
        self.assertGreater(variance,9/16)
        hc=3/(16*(1-fr.RR)**2)
        q=(.5-fr.RR)/(1-fr.RR)
        self.assertTrue(0<=2*q<=1)
        rho=q*(fr.PPLUS+fr.PMINUS)/2+(1-2*q)*fr.PSYM/4
        Eplus=fr.RR*fr.I8+(1-fr.RR)*fr.PPLUS
        self.assertAlmostEqual(np.trace(rho).real,1.,places=13)
        self.assertAlmostEqual(fr.trace_real(Eplus@rho),.5,places=13)
        REPORT['sampling']={'finite_u_calibration':calibration,'mixed_input_return_probabilities':probs,
            'finite_u_mixed_input_n_variance':float(variance),'direct_pauli_worst_n_variance':9/16,
            'hadamard_worst_n_variance':1/4,'projective_worst_n_variance':3/16,
            'hard_core_return_worst_n_variance':float(hc),
            'scope':'Independent preparations with stated allocations; ideal detection and calibrated responses. No hardware-time comparison, all-protocol optimum, or sample advantage claim.'}

    def test_05_detector_misclassification_budget(self):
        u=16.;bit_error=.003;cal=fr.calibration(u);w=cal['weights']
        occ=fr.model(-fr.PHI)[4]
        ideal=np.array([parity_acceptance(v) for v in occ]);noisy=np.array([parity_acceptance(v,bit_error) for v in occ])
        b=float(max(abs(noisy-ideal)))
        contrasts=[];rows=[]
        for tau in cal['times']:
            pair=[]
            for phi in (-fr.PHI,fr.PHI):
                E=effect(u,tau,phi);En=effect(u,tau,phi,bit_error)
                err=float(np.linalg.norm(En-E,2));self.assertLessEqual(err,b+2e-13)
                pair.append(En);rows.append({'tau':tau,'phi':phi,'effect_error':err})
            contrasts.append(pair[0]-pair[1])
        reconstruction=sum(x*e for x,e in zip(w,contrasts))
        err=float(np.linalg.norm(reconstruction-fr.CHI,2));bound=2*sum(abs(x) for x in w)*b
        self.assertLessEqual(err,bound+2e-13)
        self.assertGreater(err,1e-4)
        REPORT['misclassification']={'bit_error_probability':bit_error,'parity_effect_error_bound':b,
            'individual_effects':rows,'chirality_operator_bias':err,'chirality_bias_bound':bound,
            'scope':'Independent classical parity-bit flips only. This diagnostic does not model coherent detection errors or certify an apparatus.'}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    out=args.output.expanduser().resolve()
    if out.exists() or out==ROOT or ROOT in out.parents:ap.error('Use a new output path outside this package.')
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Checks))
    REPORT.update(status='PASS' if result.wasSuccessful() else 'FAIL',test_groups=result.testsRun,
                  scope='Author-side analytical/numerical resource audit, not external proof review, apparatus validation, or interval certification.',
                  previous_checkers_imported=False,preserved_hamiltonian_helper_reused=True)
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as f:json.dump(REPORT,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n')
    raise SystemExit(0 if result.wasSuccessful() else 1)
if __name__=='__main__':main()
