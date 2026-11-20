"""K=10: optimize ONE four-effect POVM, held fixed from t=0 through kick 3.

The objective is the least-squares slope of four post-kick stroboscopic OE
samples n=0,1,2,3. Every four-point set spans at most four dimensions, so
twirling on its orthogonal complement exactly reduces unrestricted d-dimensional
POVM optimization to four-by-four effects plus scalar complement weights.
"""
import json
import time
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import kicked_top_instantaneous_povm_pga as pg
from canonical_instantaneous_obsen_coefficient import simplex
from kicked_rotor_K10_postkick_instantaneous import HBAR,K,setup,free,kick

OUT=Path(__file__).resolve().parents[1]/'results'
PREFIX='kicked_rotor_K10_fixed_povm_kicks0to3_m4'
SEED=20261126
M=4


def make_states(d):
    q,p,H,psi0=setup(d)
    states=[psi0]
    for _ in range(3):
        states.append(kick(free(states[-1],H,1),q))
    states=np.stack(states)
    assert np.max(abs(np.linalg.norm(states,axis=1)-1))<1e-11
    F,_=np.linalg.qr(states.T)
    r=F.shape[1]
    assert np.linalg.norm(F.conj().T@F-np.eye(r))<1e-10
    z=states@F.conj()
    assert np.linalg.norm(z@F.T-states)<1e-9
    rho=np.einsum('ta,tb->tab',z,z.conj())
    return q,p,H,states,F,rho


def objective(A,alpha,d,rho):
    r=A.shape[-1];n=d-r
    p=np.einsum('iab,tba->ti',A,rho).real
    V=np.trace(A,axis1=1,axis2=2).real+n*alpha
    if np.min(p)<=1e-14 or np.min(V)<=1e-14:
        return None
    logs=np.log(p)-np.log(V[None,:])
    ent=-np.sum(p*logs,axis=1)
    weights=np.arange(len(ent),dtype=float)
    weights=(weights-weights.mean())/sum((weights-weights.mean())**2)
    score=float(weights@ent)
    grad_A=(-np.einsum('t,ti,tab->iab',weights,logs+1,rho)
            +np.sum(weights[:,None]*p/V[None,:],axis=0)[:,None,None]*np.eye(r)[None,:,:])
    grad_alpha=np.sum(weights[:,None]*p/V[None,:],axis=0)
    return score,ent,pg.herm(grad_A),grad_alpha


def projection(A,alpha,eta,G,g):
    return pg.dykstra_projection(A+eta*G,tol=1e-9,maxsteps=4000)[0],simplex(alpha+eta*g)


def residual(A,alpha,G,g,d,step=.03):
    E,a=projection(A,alpha,step,G,g)
    return float(np.sqrt(np.linalg.norm((E-A)/step)**2
                         +(d-len(A[0]))*np.linalg.norm((a-alpha)/step)**2))


def optimize(d,rho,restarts,iterations,seed):
    pg.M=M;pg.D=rho.shape[-1]
    rng=np.random.default_rng(seed)
    best=None;records=[]
    for restart in range(restarts):
        A=pg.random_povm(rng)
        alpha=rng.dirichlet(np.ones(M))
        trace=[];status='iteration limit'
        for step in range(iterations):
            score,ent,G,g=objective(A,alpha,d,rho)
            trace.append(score)
            norm=np.sqrt(np.linalg.norm(G)**2+(d-pg.D)*np.linalg.norm(g)**2)
            eta=min(.6,3/max(norm,1e-12))
            for _ in range(25):
                candidate,at=projection(A,alpha,eta,G,g)
                obj=objective(candidate,at,d,rho)
                if obj is not None and obj[0]>score+1e-11:
                    break
                eta*=.5
            else:
                status='no improving step';break
            A,alpha=candidate,at
            if (step+1)%100==0:
                print('d',d,'restart',restart,'step',step+1,'slope',obj[0],flush=True)
        final=objective(A,alpha,d,rho)
        score,ent,G,g=final
        res=residual(A,alpha,G,g,d)
        record={'restart':restart,'slope':score,'entropy':ent.tolist(),
                'projected_residual':res,'iterations':len(trace),
                'status':status,'trace':trace+[score]}
        records.append(record)
        print('d',d,'restart',restart,'slope',score,'residual',res,status,flush=True)
        if best is None or score>best[0]['slope']:
            best=record,A.copy(),alpha.copy()
    return records,best


def entropy(Eblock,alpha,F,psi,d):
    z=F.conj().T@psi
    weight=np.vdot(z,z).real
    p=np.einsum('a,iab,b->i',z.conj(),Eblock,z).real+alpha*(1-weight)
    v=np.trace(Eblock,axis1=1,axis2=2).real+(d-F.shape[1])*alpha
    assert min(p)>0 and abs(sum(p)-1)<1e-9
    return float(-p@(np.log(p)-np.log(v)))


def main():
    OUT.mkdir(exist_ok=True)
    target=OUT/f'{PREFIX}.json';figure=OUT/f'{PREFIX}.png';archive=OUT/f'{PREFIX}.npz'
    if any(x.exists() for x in (target,figure,archive)):
        raise FileExistsError('Refusing to overwrite 0-to-3 fixed-POVM optimization')
    config={'K':K,'hbar_eff':HBAR,'outcomes':M,'seed':SEED,
            'objective':'OLS slope of one fixed POVM on S(n+) at n=0,1,2,3',
            'time_convention':'prepare t=0; free to t=1; delta kicks at t=1,2,3',
            'solver':'Dykstra 4x4 POVM + complement simplex; CVXPY unavailable',
            'results':[]}
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    output={}
    for ax,d in zip(axes,(64,768)):
        q,p,H,states,F,rho=make_states(d)
        start=time.monotonic()
        records,(winner,A,alpha)=optimize(d,rho,12,1100,SEED+d)
        score,ent,G,g=objective(A,alpha,d,rho)
        # Ensure strict PSD despite Dykstra's finite tolerance.
        eps=1e-7
        A=(1-eps)*A+eps*np.eye(F.shape[1])[None,:,:]/M
        alpha=(1-eps)*alpha+eps*np.ones(M)/M
        score,ent,G,g=objective(A,alpha,d,rho)
        min_eig=np.linalg.eigvalsh(A)[:,0]
        completeness=np.linalg.norm(A.sum(axis=0)-np.eye(F.shape[1]))
        assert min(min_eig)>0 and completeness<1e-9 and abs(sum(alpha)-1)<1e-9
        assert np.allclose(ent,[entropy(A,alpha,F,s,d) for s in states],atol=1e-10)
        slopesec=(ent[-1]-ent[0])/3
        assert abs(score-np.polyfit(np.arange(4),ent,1)[0])<1e-10
        ts=[];curve=[]
        for interval in range(3):
            for delta in np.linspace(0,1,151,endpoint=False):
                ts.append(interval+delta)
                curve.append(entropy(A,alpha,F,free(states[interval],H,delta),d))
        ts=np.array(ts);curve=np.array(curve)
        ax.plot(ts,curve,label='One fixed POVM (free intervals)')
        ax.plot(np.arange(4),ent,'ko',ms=4,label='States after kicks 0, 1, 2, 3')
        ax.plot(np.arange(4),np.polyval(np.polyfit(np.arange(4),ent,1),np.arange(4)),
                color='C3',ls='--',label=f'Optimized four-point slope {score:.3f}')
        for ki in (1,2,3):ax.axvline(ki,color='0.4',ls=':',lw=1)
        ax.set(xlabel='Time (delta kicks at 1, 2, 3)',ylabel='Observational entropy (nats)',title=f'K=10, d={d}')
        ax.grid(alpha=.3);ax.legend(fontsize=7)
        edge=[]
        for s in states:
            prob=abs(s)**2
            edge.append(float(prob[np.abs(p)>.8*max(abs(p))].sum()))
        item={'d':d,'best_OLS_slope':score,'endpoint_secant':slopesec,
              'entropies_after_kicks_0_1_2_3':ent.tolist(),
              'min_eigenvalues_reduced_effects':min_eig.tolist(),
              'reduced_completeness_error':completeness,
              'best_projected_residual':residual(A,alpha,G,g,d),
              'edge_mass_after_kicks_0_1_2_3':edge,'best_restart':winner['restart'],
              'seconds':time.monotonic()-start,'restarts':records}
        config['results'].append(item)
        print('RESULT',json.dumps({k:v for k,v in item.items() if k!='restarts'},indent=2),flush=True)
        for name,value in [('states',states),('frame',F),('effects_4x4',A),('complement_weights',alpha),
                           ('momentum',p),('time',ts),('entropy_curve',curve),('entropy_at_kicks',ent)]:
            output[f'd{d}_{name}']=value
    fig.tight_layout();fig.savefig(figure,dpi=170);plt.close(fig)
    target.write_text(json.dumps(config,indent=2))
    np.savez_compressed(archive,**output)


if __name__=='__main__':
    main()
