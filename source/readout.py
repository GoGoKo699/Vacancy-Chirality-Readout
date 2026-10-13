"""Parity readout and retained-state resource diagnostics for the four-site ring.

Loss and dephasing maps are stated diagnostic models, not apparatus fits.
"""
from __future__ import annotations
import itertools
import math
from pathlib import Path
import sys
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import finite_repulsion as fr

TARGET = (0, 1, 1, 1)


def occupations(number: int) -> tuple[tuple[int, ...], ...]:
    if type(number) is not int or number < 0:
        raise ValueError('Require a nonnegative integer particle number.')
    out = []
    for labels in itertools.combinations_with_replacement(range(8), number):
        n = [0]*8
        for label in labels:
            n[label] += 1
        out.append(tuple(n))
    return tuple(out)


def site_occupancy(n: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(n[2*j]+n[2*j+1] for j in range(4))


def parity_acceptance(n: tuple[int, ...], bit_error: float = 0.) -> float:
    if not 0 <= bit_error <= 1:
        raise ValueError('A bit error probability must lie in [0,1].')
    parity = tuple(x % 2 for x in site_occupancy(n))
    return math.prod(1-bit_error if b == target else bit_error
                     for b, target in zip(parity, TARGET))


def evolved_inputs(u: float, tau: float, phi: float) -> np.ndarray:
    ev, vec = fr.spectrum(float(u), float(phi))
    ii = fr.model(phi)[3]
    return (vec*np.exp(-1j*tau*ev))@vec[ii, :].conj().T


def effect(u: float, tau: float, phi: float, bit_error: float = 0.) -> np.ndarray:
    out = evolved_inputs(u, tau, phi)
    weights = np.array([parity_acceptance(n, bit_error) for n in fr.model(phi)[4]])
    return out.conj().T@(weights[:, None]*out)


def phase_damping(rho: np.ndarray, eta: float) -> np.ndarray:
    if rho.shape != (8, 8) or not 0 <= eta <= 1:
        raise ValueError('Use an 8x8 matrix and 0<=eta<=1.')
    hamming = np.array([[sum(x != y for x, y in zip(a, b))
                         for b in fr.SPINS] for a in fr.SPINS])
    return rho*(eta**hamming)


def damped_filter(eta: float, rounds: int) -> dict:
    if not 0 <= eta <= 1 or type(rounds) is not int or rounds < 1:
        raise ValueError('Use 0<=eta<=1 and a positive integer round count.')
    rho = fr.I8/8
    K = fr.ideal()
    for _ in range(rounds):
        rho = phase_damping(K@rho@K.conj().T, eta)
    p = float(np.trace(rho).real)
    return {'eta':eta, 'rounds':rounds, 'probability':p,
            'chirality':float(np.trace(fr.CHI@rho).real/p), 'unnormalized':rho}


def damped_recurrence(eta: float, rounds: int) -> dict:
    if not 0 <= eta <= 1 or type(rounds) is not int or rounds < 1:
        raise ValueError('Use 0<=eta<=1 and a positive integer round count.')
    lam, r = eta*eta, fr.RR
    M = np.array([[(1+2*lam)/3, 2*r*(1-lam)/3],
                  [(1-lam)/3, r*(2+lam)/3]])
    x, y = np.linalg.matrix_power(M, rounds)@np.array([1/8, 1/8])
    p = 2*(x+2*y)+r**rounds/4
    return {'eta':eta, 'rounds':rounds, 'probability':float(p),
            'chirality':float(np.sqrt(3)/2*(x-y)/p),
            'spectral_radius':float(max(abs(np.linalg.eigvals(M))))}


def sector_hamiltonian(number: int, u: float, phi: float):
    """Independent occupation construction for the 0..3 loss diagnostic sectors."""
    occ = occupations(number)
    lookup = {n:i for i,n in enumerate(occ)}
    H = np.zeros((len(occ), len(occ)), complex)
    for col, n in enumerate(occ):
        sites = site_occupancy(n)
        H[col,col] = u*sum(x*(x-1)/2 for x in sites)
        for site in range(4):
            nxt = (site+1) % 4
            for spin in (0,1):
                for src,dst,p in ((2*nxt+spin,2*site+spin,np.exp(1j*phi/4)),
                                  (2*site+spin,2*nxt+spin,np.exp(-1j*phi/4))):
                    if n[src]:
                        new = list(n); new[src]-=1; new[dst]+=1
                        H[lookup[tuple(new)],col] -= p*np.sqrt(n[src]*(n[dst]+1))
    return sparse.csr_matrix(H), occ


def loss_evolution(u: float, tau: float, rate: float, rho_spin: np.ndarray):
    """Complete number-block Lindblad evolution through N=0,1,2,3.

    No-jump dynamics is not substituted in this check. Boson removal jumps
    populate all lower sectors. No gains or detector bit errors are added.
    """
    if rate < 0 or tau < 0 or rho_spin.shape != (8,8):
        raise ValueError('Invalid loss parameters or input matrix.')
    phi = -fr.PHI
    data = [sector_hamiltonian(n, u, phi) for n in range(4)]
    dims = [H.shape[0] for H,_ in data]
    sizes = [d*d for d in dims]
    blocks = [[sparse.csr_matrix((sizes[i],sizes[j]),dtype=complex) for j in range(4)] for i in range(4)]
    for n,(H,occ) in enumerate(data):
        I = sparse.eye(dims[n],format='csr')
        blocks[n][n] = -1j*(sparse.kron(I,H)-sparse.kron(H.T,I))-rate*n*sparse.eye(sizes[n])
        if n:
            lower = {v:i for i,v in enumerate(data[n-1][1])}
            feed = sparse.csr_matrix((sizes[n-1],sizes[n]),dtype=complex)
            for mode in range(8):
                rows,cols,values = [],[],[]
                for col,v in enumerate(occ):
                    if v[mode]:
                        new = list(v);new[mode]-=1
                        rows.append(lower[tuple(new)]);cols.append(col);values.append(np.sqrt(v[mode]))
                a = sparse.csr_matrix((values,(rows,cols)),shape=(dims[n-1],dims[n]))
                feed += rate*sparse.kron(a.conj(),a,format='csr')
            blocks[n-1][n] = feed
    L = sparse.bmat(blocks,format='csr')
    lookup = {v:i for i,v in enumerate(data[3][1])}
    ii=[]
    for s in fr.SPINS:
        n=[0]*8
        for j,spin in enumerate(s,1):n[2*j+spin]=1
        ii.append(lookup[tuple(n)])
    rho = np.zeros((dims[3],dims[3]),complex);rho[np.ix_(ii,ii)] = rho_spin
    offsets=np.cumsum([0]+sizes)
    vec=np.zeros(sum(sizes),complex);vec[offsets[3]:]=rho.reshape(-1,order='F')
    final=expm_multiply(tau*L,vec)
    states=[final[offsets[n]:offsets[n+1]].reshape((dims[n],dims[n]),order='F') for n in range(4)]
    accepted=states[3][np.ix_(ii,ii)]
    parity_p=sum(sum(parity_acceptance(v)*states[n][i,i].real for i,v in enumerate(data[n][1])) for n in range(4))
    return {'accepted':accepted,'total_trace':float(sum(np.trace(x).real for x in states)),
            'sector_traces':[float(np.trace(x).real) for x in states],
            'parity_probability':float(parity_p),'liouville_dimension':L.shape[0],
            'hilbert_sectors':dims}
