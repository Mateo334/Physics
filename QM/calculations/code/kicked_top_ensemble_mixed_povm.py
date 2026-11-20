"""Twenty independent initial states (10 pure, 10 mixed): shared four-POVM.

Run: py -3.13 calculations/code/kicked_top_ensemble_mixed_povm.py
Requires NumPy/SciPy/Matplotlib and the adjacent kicked_top_instantaneous_povm_pga.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm

import kicked_top_instantaneous_povm_pga as core

COUNT=20
SEED_STATES=20261121
SEED_OPT=20261122
ITERATIONS=90
PREFIX='kicked_top_ensemble20_mixed_j10_m4'
OUT=Path(__file__).resolve().parents[1]/'results'


def preparations():
    rng=np.random.default_rng(SEED_STATES)
    z,jy,_=core.spin()
    north=np.zeros(core.D,complex)
    north[-1]=1
    theta=np.arccos(rng.uniform(-1,1,size=COUNT))
    phi=rng.uniform(0,2*np.pi,size=COUNT)
    eta=np.ones(COUNT)
    eta[10:]=rng.uniform(.25,.75,size=10)
    rho=[]
    for index in range(COUNT):
        psi=np.exp(-1j*phi[index]*z)*(expm(-1j*theta[index]*jy)@north)
        rho.append(eta[index]*np.outer(psi,psi.conj())
                   +(1-eta[index])*np.eye(core.D)/core.D)
    rho=np.stack(rho)
    assert np.allclose(np.trace(rho,axis1=1,axis2=2),1,atol=1e-12)
    assert np.min(np.linalg.eigvalsh(rho))>-1e-12
    assert np.max(abs(np.trace(rho[:10]@rho[:10],axis1=1,axis2=2)-1))<1e-11
    assert np.min(np.linalg.eigvalsh(rho[10:]))>0
    return z,jy,rho,theta,phi,eta


def evolve(rho,k,z,jy):
    rotation=expm(-1j*np.pi*jy/2)
    twist=np.exp(-1j*k*z*z/(2*core.J))
    U=twist[:,None]*rotation
    R=expm(-1j*np.pi*core.PHASE*jy)
    rho=core.herm(R@U@rho@U.conj().T@R.conj().T)
    H=np.pi*jy
    drho=core.herm(-1j*(H@rho-rho@H))
    assert np.max(abs(np.trace(rho,axis1=1,axis2=2)-1))<1e-11
    assert np.min(np.linalg.eigvalsh(rho))>-1e-10
    assert np.max(abs(np.trace(drho,axis1=1,axis2=2)))<1e-11
    return rho,drho


def individuals(E,rho,drho):
    p=np.einsum('iab,nba->ni',E,rho).real
    dp=np.einsum('iab,nba->ni',E,drho).real
    V=np.trace(E,axis1=1,axis2=2).real
    if np.min(p)<=1e-13 or np.min(V)<=1e-13:
        return None
    ratio=np.log(p)-np.log(V)[None,:]
    return -np.sum(dp*ratio,axis=1),-np.sum(p*ratio,axis=1),p,dp


def average_rate_gradient(E,rho,drho):
    per=individuals(E,rho,drho)
    if per is None:
        return None
    rates,entropy,p,dp=per
    V=np.trace(E,axis1=1,axis2=2).real
    logratio=np.log(p)-np.log(V)[None,:]
    G=(-np.einsum('ni,nab->iab',logratio,drho)
       -np.einsum('ni,nab->iab',dp/p,rho)
       +(dp/V[None,:]).sum(axis=0)[:,None,None]*np.eye(core.D)[None,:,:])/len(rho)
    return float(rates.mean()),float(entropy.mean()),core.herm(G),p.mean(axis=0),dp.mean(axis=0)


def entropy_trajectory(E,rho,jy,dt):
    H=np.pi*jy
    values=[]
    for delta in dt:
        rotation=expm(-1j*delta*H)
        moved=core.herm(rotation@rho@rotation.conj().T)
        deriv=core.herm(-1j*(H@moved-moved@H))
        values.append(individuals(E,moved,deriv)[1].mean())
    return np.array(values)


def main():
    assert core.M==4
    core.MAX_ITER=ITERATIONS
    OUT.mkdir(exist_ok=True)
    z,jy,rho0,theta,phi,eta=preparations()
    rng=np.random.default_rng(SEED_OPT)
    summary={'j':core.J,'dimension':core.D,'outcomes':core.M,
             'initial_state_count':COUNT,'pure_count':10,'mixed_count':10,
             'mixed_state_recipe':'eta |theta,phi><theta,phi| +(1-eta) I/21; eta uniform [0.25,0.75]',
             'sampling':'cos(theta) uniform [-1,1], phi uniform [0,2pi]',
             'same_states_for_both_K':True,'seed_states':SEED_STATES,'seed_optimizer':SEED_OPT,
             'theta':theta.tolist(),'phi':phi.tolist(),'eta':eta.tolist(),
             'objective':'mean of 20 individual instantaneous rates, not rate of average density',
             't0':core.T0,'phase':core.PHASE,'iterations':ITERATIONS,'restarts':core.RESTARTS,
             'projection':'CVXPY/SCS' if core.cp is not None else 'Dykstra (CVXPY unavailable)',
             'results':[]}
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for k in core.KAPPAS:
        rho,drho=evolve(rho0,k,z,jy)
        grad_error=core.gradient_check(rho,drho,rng,evaluator=average_rate_gradient)
        # Cross-check batch against the original single-state formula.
        trial=core.random_povm(rng)
        b=average_rate_gradient(trial,rho,drho)
        checks=[core.rate_entropy_gradient(trial,r,d) for r,d in zip(rho,drho)]
        assert np.allclose(b[0],np.mean([v[0] for v in checks]),atol=1e-11)
        assert np.linalg.norm(b[2]-np.mean([v[2] for v in checks],axis=0))<1e-10
        records,(winner,E)=core.optimize(rho,drho,rng,evaluator=average_rate_gradient)
        rate,S,G,*_=average_rate_gradient(E,rho,drho)
        per=individuals(E,rho,drho)
        eig=np.linalg.eigvalsh(E)[:,0]
        normerror=np.linalg.norm(E.sum(axis=0)-np.eye(core.D))
        assert min(eig)>-1e-8 and normerror<1e-8
        dt=np.linspace(-.09,.09,61)
        curve=entropy_trajectory(E,rho,jy,dt)
        h=1e-6
        fd=(entropy_trajectory(E,rho,jy,[h])[0]-entropy_trajectory(E,rho,jy,[-h])[0])/(2*h)
        assert abs(fd-rate)<1e-4,(fd,rate)
        path=OUT/f'{PREFIX}_K{int(k)}.npz'
        np.savez_compressed(path,effects=E,rho_initial=rho0,rho_t0=rho,drho_t0=drho,
                            theta=theta,phi=phi,eta=eta,individual_rates=per[0],
                            individual_entropy=per[1],dt=dt,mean_entropy_curve=curve)
        result={'K':k,'mean_rate':rate,'mean_entropy':S,
                'pure_mean_rate':float(per[0][:10].mean()),
                'mixed_mean_rate':float(per[0][10:].mean()),
                'individual_rate_min':float(per[0].min()),
                'individual_rate_max':float(per[0].max()),
                'mean_density_purity':float(np.trace(rho.mean(axis=0)@rho.mean(axis=0)).real),
                'gradient_directional_error':grad_error,'local_entropy_difference':float(fd),
                'projected_gradient_residual':winner['residual_alpha_0.05'],
                'min_eigenvalues':eig.tolist(),'normalization_error':float(normerror),
                'best_restart':winner['restart'],'restarts':records,'npz':str(path)}
        summary['results'].append(result)
        print('RESULT',json.dumps({a:b for a,b in result.items() if a!='restarts'},indent=2),flush=True)
        for rec in records:
            axes[0].plot(rec['trace'],color='C0' if k==3 else 'C1',alpha=.20)
        axes[0].plot(winner['trace'],color='C0' if k==3 else 'C1',lw=2,label=f'K={int(k)}, best')
        axes[1].plot(core.T0+dt,curve,label=f'K={int(k)}, slope={rate:.3f}')
    axes[0].set(xlabel='Projected-gradient iteration',ylabel='Mean instantaneous rate (nats/time)')
    axes[1].set(xlabel='Time in rotation interval',ylabel='Mean of individual entropies (nats)')
    for ax in axes:
        ax.grid(alpha=.3);ax.legend()
    fig.tight_layout()
    fig.savefig(OUT/f'{PREFIX}_comparison.png',dpi=170)
    plt.close(fig)
    (OUT/f'{PREFIX}_summary.json').write_text(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
