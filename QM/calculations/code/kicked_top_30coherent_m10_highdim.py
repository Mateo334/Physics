"""Thirty random spin-coherent states; one shared ten-effect POVM per K.

Dense projected-gradient reference implementation. Benchmark and choose j via
--spin; runtime-controlled and explicitly NOT a hardware-dimension maximum.
Use OPENBLAS_NUM_THREADS=1 to make runtime predictable.
"""
import argparse
import json
import time
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import kicked_top_instantaneous_povm_pga as core

N,SEED_PREP,SEED_OPT=30,20261123,20261124
OUT=Path(__file__).resolve().parents[1]/'results'


def states(j):
    core.J=j; core.D=2*j+1; core.M=10
    d=core.D
    m=np.arange(-j,j+1,dtype=float)
    raise_op=np.diag(np.sqrt(j*(j+1)-m[:-1]*(m[:-1]+1)),-1)
    Jy=(raise_op-raise_op.T)/(2j)
    eig,vec=np.linalg.eigh(Jy)
    north=vec.conj().T[:,-1]
    rng=np.random.default_rng(SEED_PREP)
    theta=np.arccos(rng.uniform(-1,1,N))
    phi=rng.uniform(0,2*np.pi,N)
    initial=vec@(np.exp(-1j*eig[:,None]*theta[None,:])*north[:,None])
    initial*=np.exp(-1j*m[:,None]*phi[None,:])
    assert np.allclose(np.sum(abs(initial)**2,axis=0),1,atol=1e-11)
    assert np.max(abs(np.sum(m[:,None]*abs(initial)**2,axis=0)-j*np.cos(theta)))<1e-8
    return m,Jy,theta,phi,initial


def evolved_states(k,m,Jy,initial):
    eig,V=np.linalg.eigh(Jy)
    phase=np.exp(-1j*np.pi*(.5+core.PHASE)*eig)
    rotated=V@(phase[:,None]*(V.conj().T@initial))
    # Phase .1 of second rotation acts after first twist, not before it.
    first=V@(np.exp(-1j*.5*np.pi*eig)[:,None]*(V.conj().T@initial))
    twist=np.exp(-1j*k*m*m/(2*core.J))[:,None]
    Rphase=V@(np.exp(-1j*np.pi*core.PHASE*eig)[:,None]
              *(V.conj().T@(twist*first)))
    assert np.allclose(np.sum(abs(Rphase)**2,axis=0),1,atol=1e-10)
    assert abs(np.linalg.norm(rotated)-np.linalg.norm(initial))<1e-9
    psi=Rphase.T.copy()
    rho=np.einsum('na,nb->nab',psi,psi.conj())
    H=np.pi*Jy
    drho=core.herm(-1j*(H@rho-rho@H))
    assert np.max(abs(np.trace(rho,axis1=1,axis2=2)-1))<1e-10
    return psi,rho,drho


def ensemble(E,rho,drho):
    p=np.einsum('iab,nba->ni',E,rho).real
    dp=np.einsum('iab,nba->ni',E,drho).real
    V=np.trace(E,axis1=1,axis2=2).real
    if np.min(p)<=1e-13 or np.min(V)<=1e-13:
        return None
    logratio=np.log(p)-np.log(V)[None,:]
    G=(-np.einsum('ni,nab->iab',logratio,drho)
       -np.einsum('ni,nab->iab',dp/p,rho)
       +(dp/V[None,:]).sum(axis=0)[:,None,None]*np.eye(core.D)[None,:,:])/N
    return float((-dp*logratio).sum(axis=1).mean()),float((-p*logratio).sum(axis=1).mean()),core.herm(G),p,dp


def optimize(rho,drho,rng,restarts,steps,budget):
    records=[]; best=None
    start=time.monotonic()
    for restart in range(restarts):
        if restart and time.monotonic()-start>budget:
            print('Time budget exhausted; no further restart launched',flush=True)
            break
        E=core.random_povm(rng)
        trace=[]; status='iteration limit'
        for n in range(steps):
            if n>0 and time.monotonic()-start>budget:
                status='time budget';break
            rate,_,G,*_=ensemble(E,rho,drho)
            trace.append(rate)
            eta=min(.4,2./max(np.linalg.norm(G),1e-12))
            for backtrack in range(18):
                trial,_=core.project(E+eta*G)
                value=ensemble(trial,rho,drho)
                if value is not None and value[0]>rate+1e-11:
                    break
                eta*=.5
            else:
                status='no improving step';break
            E=trial
            if (n+1)%20==0:
                print('restart',restart,'iteration',n+1,'rate',value[0],
                      'elapsed',round(time.monotonic()-start,1),flush=True)
        R,S,G,*_=ensemble(E,rho,drho)
        residual=core.residual(E,G)
        record={'restart':restart,'rate':R,'entropy':S,'projected_residual_alpha_0.05':residual,
                'iterations':len(trace),'status':status,'trace':trace+[R]}
        records.append(record)
        print('END restart',restart,'rate',R,'residual',residual,'status',status,flush=True)
        if best is None or R>best[0]['rate']:
            best=(record,E.copy())
    return records,best,time.monotonic()-start


def mean_entropy(E,psi,Jy,s):
    eig,V=np.linalg.eigh(Jy)
    vec=V@(np.exp(-1j*np.pi*s*eig)[:,None]*(V.conj().T@psi.T))
    p=np.einsum('an,iab,bn->ni',vec.conj(),E,vec).real
    vol=np.trace(E,axis1=1,axis2=2).real
    assert np.min(p)>0
    return float(np.mean(-np.sum(p*(np.log(p)-np.log(vol)[None,:]),axis=1)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spin',type=int,default=50)
    parser.add_argument('--restarts',type=int,default=3)
    parser.add_argument('--steps',type=int,default=120)
    parser.add_argument('--budget',type=float,default=180.,help='seconds per K; finish last started step')
    args=parser.parse_args()
    if args.spin<10 or args.restarts<1 or args.steps<1:
        parser.error('spin>=10, restarts>=1, steps>=1')
    prefix=f'kicked_top_30coherent_j{args.spin}_m10'
    OUT.mkdir(exist_ok=True)
    if any((OUT/f'{prefix}_{suffix}').exists() for suffix in ['summary.json','K3.npz','K6.npz']):
        parser.error('Output exists; preserve results and choose another run path')
    m,Jy,theta,phi,initial=states(args.spin)
    rng=np.random.default_rng(SEED_OPT)
    summary={'spin':args.spin,'d':core.D,'outcomes':10,'preparations':N,'time':core.T0,
             'rotation_phase':core.PHASE,'theta':theta.tolist(),'phi':phi.tolist(),
             'seed_preparations':SEED_PREP,'seed_optimizer':SEED_OPT,
             'requested_restarts':args.restarts,'max_steps':args.steps,'budget_seconds_per_K':args.budget,
             'solver':'CVXPY/SCS' if core.cp is not None else 'Dykstra Euclidean projection; CVXPY missing',
             'results':[]}
    fig,ax=plt.subplots(1,2,figsize=(11,4))
    for k in (3,6):
        psi,rho,drho=evolved_states(k,m,Jy,initial)
        check=core.gradient_check(rho,drho,rng,evaluator=ensemble)
        records,(winner,E),elapsed=optimize(rho,drho,rng,args.restarts,args.steps,args.budget)
        # The stopping tolerance for Dykstra permits tiny negative eigenvalues.
        # This strictly interior adjustment preserves completeness exactly.
        epsilon=1e-6
        E=(1-epsilon)*E+epsilon*np.eye(core.D)[None,:,:]/core.M
        R,S,G,p,dp=ensemble(E,rho,drho)
        mineig=np.linalg.eigvalsh(E)[:,0]
        normerr=np.linalg.norm(E.sum(axis=0)-np.eye(core.D))
        assert min(mineig)>0 and normerr<1e-8
        h=1e-6
        central=(mean_entropy(E,psi,Jy,h)-mean_entropy(E,psi,Jy,-h))/(2*h)
        assert abs(R-central)<2e-4,(R,central)
        resid=core.residual(E,G)
        grid=np.linspace(-.09,.09,61)
        curve=np.array([mean_entropy(E,psi,Jy,s) for s in grid])
        output=OUT/f'{prefix}_K{k}.npz'
        np.savez_compressed(output,effects=E,initial=initial,rho_t0=rho,
                            drho_t0=drho,theta=theta,phi=phi,individual_rates=(-dp*(np.log(p)-np.log(np.trace(E,axis1=1,axis2=2).real)[None,:])).sum(axis=1),
                            dt=grid,entropy=curve)
        item={'K':k,'rate':R,'entropy':S,'min_effect_eigenvalues':mineig.tolist(),
              'completeness_error':normerr,'gradient_directional_error':check,
              'entropy_central_difference':central,'projected_residual':resid,
              'best_restart':winner['restart'],'actual_attempted_restarts':len(records),
              'solver_time_seconds':elapsed,
              'all_restarts':records,'effects_file':str(output)}
        summary['results'].append(item)
        print('RESULT',json.dumps({key:value for key,value in item.items() if key!='all_restarts'},indent=2),flush=True)
        for r in records:
            ax[0].plot(r['trace'],color='C0' if k==3 else 'C1',alpha=.18)
        ax[0].plot(winner['trace'],color='C0' if k==3 else 'C1',lw=2,label=f'K={k}, best')
        ax[1].plot(core.T0+grid,curve,label=f'K={k}, slope={R:.3f}')
    ax[0].set(xlabel='Projected-gradient iteration',ylabel='30-state mean instantaneous rate (nats/time)')
    ax[1].set(xlabel='Time within rotation interval',ylabel='Mean of 30 individual OEs (nats)')
    for a in ax:
        a.grid(alpha=.3);a.legend()
    fig.tight_layout()
    fig.savefig(OUT/f'{prefix}_comparison.png',dpi=160)
    plt.close(fig)
    (OUT/f'{prefix}_summary.json').write_text(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
