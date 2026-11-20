"""Fixed four-effect POVM: maximize secant minus derivative-variance penalty.

For K=10,d=64, prepare at t=0, free to first kick at t=1; sample 50
equally spaced free-time points from 1+ through 2- and choose one POVM.
The free-energy eigenspace orbit gives an exact lower-dimensional representation
of all full-d POVMs after twirling the orthogonal complement.
"""
import argparse
import json
import time
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import kicked_top_instantaneous_povm_pga as pg
from canonical_instantaneous_obsen_coefficient import simplex
from kicked_rotor_K10_postkick_instantaneous import K,HBAR,setup,free,kick

OUT=Path(__file__).resolve().parents[1]/'results'
D=64
M=4
SEED=20261128
PREFIX='kicked_rotor_K10_d64_fixed_povm_linearity_mu1_50'


def orbit_data(num_samples):
    q,p,H,initial=setup(D)
    at_start=kick(free(initial,H,1),q)
    labels=np.round(H/HBAR*2).astype(int) # n², exact integer energies on FFT momentum grid.
    vectors=[]
    for label in np.unique(labels):
        z=np.where(labels==label,at_start,0)
        norm=np.linalg.norm(z)
        if norm>1e-14:
            vectors.append(z/norm)
    F=np.column_stack(vectors)
    r=F.shape[1]
    assert np.linalg.norm(F.conj().T@F-np.eye(r))<1e-11
    assert np.linalg.norm(F@(F.conj().T@at_start)-at_start)<1e-11
    t=np.linspace(0,1,num_samples)
    states=np.exp(-1j*t[:,None]*H[None,:])*at_start[None,:]
    z=states@F.conj()
    assert np.max(np.linalg.norm(z@F.T-states,axis=1))<1e-10
    rho=np.einsum('ta,tb->tab',z,z.conj())
    h=F.conj().T@(H[:,None]*F)
    assert np.linalg.norm(h-h.conj().T)<1e-11
    drho=pg.herm(-1j*(h@rho-rho@h))
    weights=np.ones(num_samples)/(num_samples-1)
    weights[0]/=2;weights[-1]/=2
    assert np.isclose(weights.sum(),1)
    return t,weights,rho,drho,F,at_start,H,initial


def statistics(A,alpha,d,rho,drho):
    p=np.einsum('iab,tba->ti',A,rho).real
    dp=np.einsum('iab,tba->ti',A,drho).real
    V=np.trace(A,axis1=1,axis2=2).real+(d-A.shape[-1])*alpha
    if min(p.min(),V.min())<=1e-14:
        return None
    logs=np.log(p)-np.log(V)[None,:]
    entropy=-np.sum(p*logs,axis=1)
    ds=-np.sum(dp*logs,axis=1)
    return entropy,ds,p,dp,V,logs


def objective(A,alpha,mu,weights,rho,drho,derivative=True):
    stats=statistics(A,alpha,D,rho,drho)
    if stats is None:
        return None
    entropy,ds,p,dp,V,logs=stats
    secant=float(entropy[-1]-entropy[0])
    diff=ds-secant
    penalty=float(weights@(diff*diff))
    score=secant-mu*penalty
    if not derivative:
        return score,secant,penalty,entropy,ds
    r=rho.shape[-1];n=D-r
    eye=np.eye(r)[None,:,:]
    # Matrix derivatives of S and dot S, in the Hilbert--Schmidt pairing.
    # Use endpoint differences, rather than sums, for the secant gradient.
    gR=(-(logs[-1,:,None,None]+1)*rho[-1,None,:,:]
         +(p[-1]/V)[:,None,None]*eye
         +(logs[0,:,None,None]+1)*rho[0,None,:,:]
         -(p[0]/V)[:,None,None]*eye)
    gR_alpha=(p[-1]-p[0])/V
    gdot=(-logs[:,:,None,None]*drho[:,None,:,:]
          -(dp/p)[:,:,None,None]*rho[:,None,:,:]
          +(dp/V[None,:])[:,:,None,None]*eye[None,:,:,:])
    gdot_alpha=dp/V[None,:]
    coeff=2*mu*weights*diff
    gA=pg.herm((1+coeff.sum())*gR-np.einsum('t,tiab->iab',coeff,gdot))
    ga=(1+coeff.sum())*gR_alpha-np.einsum('t,ti->i',coeff,gdot_alpha)
    assert gA.shape==(M,r,r) and ga.shape==(M,)
    return score,secant,penalty,entropy,ds,gA,ga


def project(A,alpha,eta,G,g):
    nextA,_=pg.dykstra_projection(A+eta*G,tol=1e-8,maxsteps=3500)
    return nextA,simplex(alpha+eta*g)


def residual(A,alpha,G,g,step=.03):
    E,b=project(A,alpha,step,G,g)
    return float(np.sqrt(np.linalg.norm((E-A)/step)**2
                         +(D-A.shape[-1])*np.linalg.norm((b-alpha)/step)**2))


def check_gradient(rho,drho,weights,mu,rng):
    pg.M=M;pg.D=rho.shape[-1]
    A=pg.random_povm(rng);alpha=rng.dirichlet(np.ones(M))
    DA=pg.herm(rng.normal(size=A.shape)+1j*rng.normal(size=A.shape))
    DA/=np.linalg.norm(DA)
    da=rng.normal(size=M);da/=np.linalg.norm(da)
    current=objective(A,alpha,mu,weights,rho,drho)
    G,g=current[-2:]
    directional=np.einsum('iab,iba->',G,DA).real+(D-pg.D)*np.dot(g,da)
    h=1e-6
    numerical=(objective(A+h*DA,alpha+h*da,mu,weights,rho,drho,False)[0]
               -objective(A-h*DA,alpha-h*da,mu,weights,rho,drho,False)[0])/(2*h)
    error=abs(numerical-directional)/max(abs(numerical),1e-12)
    assert error<5e-6,(numerical,directional,error)
    return float(error)


def optimize(rho,drho,weights,mu,restarts,steps,rng):
    pg.M=M;pg.D=rho.shape[-1]
    results=[];best=None
    for restart in range(restarts):
        A=pg.random_povm(rng)
        alpha=rng.dirichlet(np.ones(M))
        history=[];status='iteration limit'
        for iteration in range(steps):
            score,rate,penalty,_,_,G,g=objective(A,alpha,mu,weights,rho,drho)
            history.append(score)
            norm=np.sqrt(np.linalg.norm(G)**2+(D-pg.D)*np.linalg.norm(g)**2)
            eta=min(.4,1.5/max(norm,1e-12))
            for _ in range(20):
                X,b=project(A,alpha,eta,G,g)
                val=objective(X,b,mu,weights,rho,drho,False)
                if val is not None and val[0]>score+1e-11:
                    break
                eta*=.5
            else:
                status='no improving step';break
            A,alpha=X,b
            if (iteration+1)%100==0:
                print('restart',restart,'step',iteration+1,'objective',val[0],flush=True)
        final=objective(A,alpha,mu,weights,rho,drho)
        score,rate,penalty,ent,ds,G,g=final
        res=residual(A,alpha,G,g)
        item={'restart':restart,'objective':score,'secant_rate':rate,'penalty_integral':penalty,
              'rms_derivative_deviation':np.sqrt(penalty),
              'projected_residual':res,'steps':len(history),'status':status,
              'history':history+[score]}
        results.append(item)
        print('END restart',restart,'score',score,'R',rate,'rms',np.sqrt(penalty),
              'residual',res,status,flush=True)
        if best is None or score>best[0]['objective']:
            best=(item,A.copy(),alpha.copy())
    return results,best


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mu',type=float,default=1.)
    parser.add_argument('--restarts',type=int,default=10)
    parser.add_argument('--steps',type=int,default=450)
    args=parser.parse_args()
    if args.mu!=1.:
        prefix=f'kicked_rotor_K10_d64_fixed_povm_linearity_mu{args.mu:g}_50'
    else:
        prefix=PREFIX
    if args.mu<0 or args.restarts<2 or args.steps<1:
        parser.error('mu>=0, restarts>=2, steps>=1 required')
    OUT.mkdir(exist_ok=True)
    file=OUT/f'{prefix}.npz';summary=OUT/f'{prefix}.json';plot=OUT/f'{prefix}.png'
    if any(f.exists() for f in (file,summary,plot)):
        parser.error('Output exists; refusing to overwrite')
    t,weights,rho,drho,F,state,H,initial=orbit_data(50)
    rng=np.random.default_rng(SEED)
    grad_error=check_gradient(rho,drho,weights,args.mu,rng)
    start=time.monotonic()
    records,(winner,A,alpha)=optimize(rho,drho,weights,args.mu,args.restarts,args.steps,rng)
    # Strict positive effects after finite-tolerance projection.
    eps=1e-7
    A=(1-eps)*A+eps*np.eye(A.shape[-1])[None,:,:]/M
    alpha=(1-eps)*alpha+eps*np.ones(M)/M
    score,R,P,S,ds,G,g=objective(A,alpha,args.mu,weights,rho,drho)
    residual_final=residual(A,alpha,G,g)
    assert np.linalg.norm(A.sum(axis=0)-np.eye(A.shape[-1]))<1e-9
    assert np.min(np.linalg.eigvalsh(A))>0
    assert np.all(alpha>0) and abs(alpha.sum()-1)<1e-10
    quadrature_error=float(np.trapezoid(ds,t)-R) # Sampling diagnostic, not an identity at finite quadrature.
    t_fine,w_fine,rho_fine,drho_fine,*_=orbit_data(501)
    _,Rf,Pfine,Sfine,dsfine=objective(A,alpha,args.mu,w_fine,rho_fine,drho_fine,False)
    assert np.isclose(R,Rf,atol=1e-9)
    fig,ax=plt.subplots(1,2,figsize=(11,4.3))
    ax[0].plot(1+t_fine,Sfine,label='One fixed POVM: entropy')
    ax[0].plot(1+t_fine,Sfine[0]+R*t_fine,ls='--',label=f'Endpoint line, R={R:.4f}')
    ax[0].scatter(1+t,S,color='C3',s=10,label='50 objective samples')
    ax[1].plot(1+t_fine,dsfine,label='Fixed-POVM derivative')
    ax[1].axhline(R,color='C3',ls='--',label=f'Interval mean R={R:.4f}')
    ax[0].set(xlabel='Time, first kick at 1 and second at 2',ylabel='Observational entropy (nats)')
    ax[1].set(xlabel='Time within kick-free interval',ylabel='Derivative (nats/free time)')
    for a in ax:
        a.grid(alpha=.3);a.legend(fontsize=8)
    fig.suptitle(f'K=10, d=64, four fixed outcomes, mu={args.mu:g}')
    fig.tight_layout();fig.savefig(plot,dpi=170);plt.close(fig)
    data={'K':K,'d':D,'outcomes':M,'mu':args.mu,'sample_count':len(t),'quadrature':'trapezoid',
          'kick_interval':'1+ through 2-','free_hbar_eff':HBAR,
          'subspace_dimension':F.shape[-1],'gradient_directional_relative_error':grad_error,
          'best_score':score,'best_R_secant':R,'penalty_integral_50':P,
          'derivative_quadrature_minus_exact_secant_50':quadrature_error,
          'rms_derivative_deviation_50':float(np.sqrt(P)),
          'penalty_integral_501_check':Pfine,'rms_derivative_deviation_501_check':float(np.sqrt(Pfine)),
          'entropy_start':S[0],'entropy_end':S[-1],
          'max_derivative_deviation_501':float(np.max(abs(dsfine-R))),
          'projected_gradient_residual':residual_final,
          'effect_minimum_eigenvalues_4x4':np.linalg.eigvalsh(A)[:,0].tolist(),
          'actual_independent_starts':len(records),'iterations_cap':args.steps,
          'optimizer_seconds':time.monotonic()-start,
          'lyapunov_reference_quoted':1.619,
          'edge_probability_after_first_kick':float(np.sum(abs(state[np.abs(np.fft.fftfreq(D,d=1/D)*HBAR)>
                                                       .8*np.max(abs(np.fft.fftfreq(D,d=1/D)*HBAR))])**2)),
          'objective_not_globally_certified':True,'restarts':records}
    np.savez_compressed(file,frame=F,two_by_two_effects=A,complement_weights=alpha,
                        time_samples=1+t,entropy_samples=S,derivative_samples=ds,
                        time_fine=1+t_fine,entropy_fine=Sfine,derivative_fine=dsfine,
                        psi_after_first=state,initial_state=initial)
    summary.write_text(json.dumps(data,indent=2))
    print('RESULT',json.dumps({k:v for k,v in data.items() if k!='restarts'},indent=2),flush=True)


if __name__=='__main__':
    main()
