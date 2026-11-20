"""Gradient optimization of a fixed final-only POVM for kicked-top OE slopes.

The slope is an ordinary least-squares fit to S_obs(t), t=0,...,T; this
is not an instantaneous continuous-time derivative or a measured history.
Run: py -3.13 calculations/code/optimize_kicked_top_obsen_rate.py
"""
import json
from pathlib import Path

import numpy as np
import torch
from scipy.linalg import expm

J, M, SEED = 20, 4, 20261120
KAPPAS = (3.0, 6.0)
STEPS, RESTARTS, LR = 100, 2, 0.05
OUT = Path(__file__).resolve().parents[1] / 'results'


def trajectory(k, T):
    m = np.arange(-J, J+1, dtype=float)
    raising = np.diag(np.sqrt(J*(J+1)-m[:-1]*(m[:-1]+1)), -1)
    jy = (raising-raising.T)/(2j)
    rotation = expm(-1j*np.pi*jy/2)
    theta, phi = 1.1, 0.7
    north = np.zeros(len(m), dtype=complex)
    north[-1] = 1
    state = np.exp(-1j*phi*m) * (expm(-1j*theta*jy) @ north)
    U = np.exp(-1j*k*m**2/(2*J))[:, None] * rotation
    assert np.allclose(U.conj().T@U, np.eye(len(m)), atol=1e-11)
    states = [state]
    for _ in range(T):
        states.append(U @ states[-1])
    states = np.stack(states)
    assert np.allclose(np.linalg.norm(states, axis=1), 1, atol=1e-11)
    return states


def effects(A):
    B = A.conj().transpose(-2,-1) @ A
    matrix = B.sum(dim=0)
    eig, vectors = torch.linalg.eigh(matrix)
    inverse = (vectors * eig.clamp_min(1e-12).rsqrt().unsqueeze(0)) @ vectors.conj().T
    return inverse.unsqueeze(0) @ B @ inverse.unsqueeze(0)


def curve(A, states):
    E = effects(A)
    probs = torch.einsum('ta,iab,tb->ti', states.conj(), E, states).real
    volumes = E.diagonal(dim1=-2,dim2=-1).real.sum(dim=-1)
    assert torch.all(volumes > 0)
    return -(probs * (probs.clamp_min(1e-300).log()-volumes.log())).sum(dim=1)


def slope(entropies):
    t = torch.arange(len(entropies), dtype=entropies.dtype)
    weights = (t-t.mean()) / ((t-t.mean())**2).sum()
    return weights @ entropies


def optimize(k, T, rng):
    psi_np = trajectory(k, T)
    psi = torch.tensor(psi_np, dtype=torch.complex128)
    d = len(psi_np[0])
    candidates = []
    for restart in range(RESTARTS):
        noise = rng.normal(size=(M,d,d)) + 1j*rng.normal(size=(M,d,d))
        initial = np.eye(d)[None,:,:]/np.sqrt(M) + 0.14*noise/np.sqrt(d)
        A = torch.nn.Parameter(torch.tensor(initial, dtype=torch.complex128))
        opt = torch.optim.Adam([A], lr=LR)
        with torch.no_grad():
            baseline = float(slope(curve(A, psi)))
        best = (baseline, A.detach().clone())
        for _ in range(STEPS):
            opt.zero_grad()
            objective = -slope(curve(A, psi))
            objective.backward()
            opt.step()
            with torch.no_grad():
                value = float(slope(curve(A, psi)))
                if np.isfinite(value) and value > best[0]:
                    best = (value, A.detach().clone())
        candidates.append((best[0], best[1], baseline))
    winner = max(candidates, key=lambda v: v[0])
    A = winner[1]
    E = effects(A).detach().numpy()
    ent = curve(A,psi).detach().numpy()
    probs = np.einsum('ta,iab,tb->ti', psi_np.conj(), E, psi_np).real
    assert np.linalg.norm(E.sum(axis=0)-np.eye(d)) < 1e-9
    assert np.min(np.linalg.eigvalsh(E)) > -1e-9
    assert np.max(np.abs(probs.sum(axis=1)-1)) < 1e-9
    assert np.all(ent >= -1e-8) and np.all(ent <= np.log(d)+1e-8)
    assert np.isclose(np.polyfit(np.arange(T+1),ent,1)[0], winner[0], atol=1e-9)
    return {'k': k, 'j': J, 'dimension':d, 'outcomes':M, 'T':T,
            'tE_estimate':np.log(d)/np.log(k/2), 'slope':winner[0],
            'restart_slopes':[v[0] for v in candidates],
            'initial_slopes':[v[2] for v in candidates],
            'entropies':ent.tolist(), 'volumes':np.trace(E,axis1=1,axis2=2).real.tolist()}, E


def main():
    torch.set_num_threads(1)
    rng = np.random.default_rng(SEED)
    OUT.mkdir(parents=True, exist_ok=True)
    results = []
    for k in KAPPAS:
        T = int(np.floor(np.log(2*J+1)/np.log(k/2)))
        result, E = optimize(k,T,rng)
        filename = OUT / f'kicked_top_obsen_optimized_k{int(k)}_j{J}_m{M}.npz'
        np.savez_compressed(filename, effects=E, entropy=result['entropies'], k=k, j=J, seed=SEED)
        result['effects_file'] = str(filename.relative_to(OUT.parent.parent))
        print(json.dumps(result,indent=2), flush=True)
        results.append(result)
    summary = {'seed':SEED,'restarts':RESTARTS,'steps':STEPS,'optimizer':'Adam',
               'learning_rate':LR,'theta':1.1,'phi':0.7,'numpy':np.__version__,
               'torch':torch.__version__,'results':results}
    (OUT / f'kicked_top_obsen_optimized_j{J}_m{M}.json').write_text(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
