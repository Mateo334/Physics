"""K=10 and K=0: one fixed four-POVM over three complete kicks, d=21 torus.

J=OLS_slope-mu*sum_n(S_n-S_0-n*OLS_slope)^2-gamma*S_0^2
  -beta*max(0,S_T-(log(d)-margin))^2 (optional saturation safeguard).
No measurement backaction. Output files never overwritten. CLI options allow
additional starts and longer runs with separate prefixes.
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import kicked_top_instantaneous_povm_pga as pg

OUT=Path(__file__).resolve().parents[1]/'results'
D=21;M=4;HBAR=2*np.pi/D;Q0=np.pi/2;P0=0.;SIGMA=np.sqrt(HBAR/2)
SEED=20261129


def prep(d):
    q=2*np.pi*np.arange(d)/d
    mode=np.fft.fftfreq(d,d=1/d)
    p=HBAR*mode
    wave=sum(np.exp(-(q-Q0+2*np.pi*j)**2/(4*SIGMA**2)+1j*P0*(q-Q0+2*np.pi*j)/HBAR)
             for j in range(-3,4))
    wave=np.asarray(wave,complex)/np.linalg.norm(wave)
    psi=np.fft.fft(wave,norm='ortho')
    return q,mode,p,psi


def evolve(d,K,T):
    q,mode,p,psi=prep(d)
    kinetic=np.exp(-1j*p*p/(2*HBAR))
    kick=np.exp(-1j*K*np.cos(q)/HBAR)
    states=[psi]
    for _ in range(T):
        psi=np.fft.fft(kick*np.fft.ifft(kinetic*psi,norm='ortho'),norm='ortho')
        states.append(psi)
    states=np.array(states)
    assert np.max(np.abs(np.linalg.norm(states,axis=1)-1))<1e-12
    rho=np.einsum('ta,tb->tab',states,states.conj())
    return q,mode,p,states,rho


def evaluate(E,rho,mu,gamma,beta=0.,margin=.35,gradient=True):
    d=rho.shape[-1]
    p=np.einsum('iab,tba->ti',E,rho).real
    v=np.trace(E,axis1=1,axis2=2).real
    if min(p.min(),v.min())<=1e-13:
        return None
    L=np.log(p)-np.log(v)[None,:]
    S=-np.sum(p*L,axis=1)
    t=np.arange(len(S),dtype=float)
    c=(t-t.mean())/np.sum((t-t.mean())**2)
    R=float(c@S)
    e=S-S[0]-t*R
    excess=max(0.,float(S[-1]-(np.log(d)-margin)))
    J=float(R-mu*e@e-gamma*S[0]**2-beta*excess**2)
    if not gradient:
        return J,R,S,e,p,v
    # dS_n/dE_i = -(log(p_ni/v_i)+1)rho_n+(p_ni/v_i)I.
    dS=pg.herm(-(L[:,:,None,None]+1)*rho[:,None,:,:]
                +(p/v[None,:])[:,:,None,None]*np.eye(d)[None,None,:,:])
    dR=np.einsum('t,tiab->iab',c,dS)
    de=dS-dS[0,None,:,:,:]-t[:,None,None,None]*dR[None,:,:,:]
    G=pg.herm(dR-2*mu*np.einsum('t,tiab->iab',e,de)-2*gamma*S[0]*dS[0]
              -2*beta*excess*dS[-1])
    return J,R,S,e,p,v,G


def aligned_start(psi,rng):
    # Exact zero-S0 seed: rank-one projector on psi, random PVM on complement.
    d=len(psi)
    Q,_=np.linalg.qr(np.column_stack([psi,rng.normal(size=(d,d-1))+
                                     1j*rng.normal(size=(d,d-1))]))
    # QR's first column equals psi up to a phase.
    E=np.zeros((M,d,d),complex)
    E[0]=np.outer(Q[:,0],Q[:,0].conj())
    perm=rng.permutation(np.arange(1,d))
    for i,indices in enumerate(np.array_split(perm,3),start=1):
        E[i]=Q[:,indices]@Q[:,indices].conj().T
    return (1-1e-5)*E+1e-5*np.eye(d)[None,:,:]/M


def optimize(rho,mu,gamma,beta,margin,restarts,steps,seed,warm=None):
    pg.D=D;pg.M=M
    rng=np.random.default_rng(seed)
    records=[];best=None
    for start in range(restarts):
        if warm is not None:
            E=warm.copy() if start==0 else (1-.03*start)*warm+.03*start*pg.random_povm(rng)
        else:
            E=aligned_start(np.linalg.eigh(rho[0])[1][:,-1],rng) if start%2==0 else pg.random_povm(rng)
        trace=[];reason='iteration limit'
        for k in range(steps):
            J,R,S,e,_,_,G=evaluate(E,rho,mu,gamma,beta,margin)
            trace.append(J)
            step=min(.5,1.2/max(np.linalg.norm(G),1e-12))
            for _ in range(22):
                trial,_=pg.dykstra_projection(E+step*G,tol=3e-9,maxsteps=2000)
                trial_val=evaluate(trial,rho,mu,gamma,beta,margin,False)
                if trial_val is not None and trial_val[0]>J+1e-12:
                    break
                step*=.5
            else:
                reason='line search stalled';break
            E=trial
            if (k+1)%200==0:
                print('start',start,'iteration',k+1,'J',trial_val[0],'S0',trial_val[2][0],flush=True)
        # Strict positivity for well-defined gradients and reproducible feasibility.
        E=(1-1e-8)*E+1e-8*np.eye(D)[None,:,:]/M
        J,R,S,e,p,v,G=evaluate(E,rho,mu,gamma,beta,margin)
        proj,_=pg.dykstra_projection(E+.01*G,tol=3e-9,maxsteps=2000)
        residual=float(np.linalg.norm(proj-E)/.01)
        vals={'start':start,'objective':J,'rate':R,'entropy0':float(S[0]),
              'fit_rms':float(np.sqrt(np.mean(e**2))),'fit_max':float(np.max(abs(e))),
              'final_entropy':float(S[-1]),'entropies':S.tolist(),
              'projected_residual_step_0.01':residual,
              'iterations':len(trace),'status':reason,'trace_tail':trace[-10:]}
        records.append(vals)
        print('FINISH',start,'J',J,'R',R,'S0',S[0],'rms',vals['fit_rms'],
              'residual',residual,reason,flush=True)
        if best is None or J>best[0]['objective']:
            best=(vals,E.copy())
    return records,best


def test_gradient(rho,mu,gamma,beta,margin,rng):
    pg.D=D;pg.M=M
    E=pg.random_povm(rng)
    direction=pg.herm(rng.normal(size=E.shape)+1j*rng.normal(size=E.shape))
    direction/=np.linalg.norm(direction)
    G=evaluate(E,rho,mu,gamma,beta,margin)[-1]
    h=1e-6
    fd=(evaluate(E+h*direction,rho,mu,gamma,beta,margin,False)[0]-
        evaluate(E-h*direction,rho,mu,gamma,beta,margin,False)[0])/(2*h)
    analytic=np.einsum('iab,iba->',G,direction).real
    rel=abs(fd-analytic)/max(1e-10,abs(fd))
    assert rel<1e-5,(fd,analytic,rel)
    return float(rel)


def cutoff_diagnostics():
    report={}
    for K in (10.,0.):
        cases={}
        for d in (21,257,513):
            q,mode,p,states,_=evolve(d,K,3)
            weights=abs(states)**2
            edge=np.abs(mode)>.8*np.max(np.abs(mode))
            # Physical momentum probabilities that a 21-mode truncation would discard.
            beyond=np.abs(mode)>10
            cases[str(d)]={'edge_mass_each_kick':weights[:,edge].sum(axis=1).tolist(),
                           'mass_outside_21_momentum_modes':weights[:,beyond].sum(axis=1).tolist(),
                           'p2_each_kick':(weights@(p*p)).tolist()}
        report[str(K)]=cases
    return report


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--T',type=int,default=3)
    parser.add_argument('--mu',type=float,default=1.)
    parser.add_argument('--gamma',type=float,default=10.)
    parser.add_argument('--beta',type=float,default=20.)
    parser.add_argument('--margin',type=float,default=.35)
    parser.add_argument('--restarts',type=int,default=8)
    parser.add_argument('--steps',type=int,default=1100)
    parser.add_argument('--tag',default='pilot')
    parser.add_argument('--warm-start',type=str,default=None,
                        help='Path to previous NPZ with K10_effects and K0_effects')
    args=parser.parse_args()
    if args.T<3 or args.mu<0 or args.gamma<0 or args.beta<0 or args.margin<0 or args.restarts<2:
        parser.error('Require T>=3, nonnegative penalties/margin, at least two restarts')
    OUT.mkdir(exist_ok=True)
    stem=f'kicked_rotor_d21_K10_K0_lowS0_T{args.T}_{args.tag}'
    paths=[OUT/f'{stem}.{ext}' for ext in ('json','png','npz')]
    if any(p.exists() for p in paths):
        parser.error('Refusing to overwrite prior data')
    result={'d':D,'hbar_eff':HBAR,'boundary':'periodic position/momentum torus',
            'momentum_modes':[-10,10],'K': [10,0],
            'preparation':{'q0':Q0,'p0':P0,'sigma_q':SIGMA},
            'T':args.T,'mu':args.mu,'gamma':args.gamma,'beta':args.beta,
            'saturation_margin':args.margin,'seed':SEED,
            'objective':'OLS slope R - mu*sum_n(Sn-S0-nR)^2 - gamma*S0^2 - beta*max(0,S_T-(ln(d)-margin))^2',
            'classical_reference_given_by_user':1.619,
            'cutoff_check':cutoff_diagnostics(),'runs':[]}
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),sharey=True)
    stored={}
    for ax,K in zip(axes,(10.,0.)):
        q,mode,p,states,rho=evolve(D,K,args.T)
        grad_error=test_gradient(rho,args.mu,args.gamma,args.beta,args.margin,
                                 np.random.default_rng(SEED+int(K)))
        warm=None
        if args.warm_start:
            with np.load(args.warm_start) as old:
                warm=old[f'K{int(K)}_effects'].copy()
        records,(winner,E)=optimize(rho,args.mu,args.gamma,args.beta,args.margin,args.restarts,args.steps,
                                    SEED+100+int(K),warm)
        J,R,S,e,prob,v,G=evaluate(E,rho,args.mu,args.gamma,args.beta,args.margin)
        assert abs(prob.sum(axis=1)-1).max()<1e-9
        assert np.linalg.norm(E.sum(axis=0)-np.eye(D))<1e-7
        assert min(np.linalg.eigvalsh(E)[:,0])>0
        t=np.arange(args.T+1)
        assert np.isclose(R,np.polyfit(t,S,1)[0],atol=1e-10)
        ax.plot(t,S,'o-',label='Fixed-POVM S_obs')
        ax.plot(t,S[0]+R*t,'--',label=f'S0+Rn, R={R:.3f}')
        ax.plot(t,S[0]+1.619*t,':',color='C3',label='S0+1.619n (reference)')
        ax.axhline(np.log(D),color='0.5',ls='-.',label='ln 21 saturation')
        ax.set(xlabel='Completed kicks n',title=f'K={K:g}, d=21')
        ax.set_xticks(t);ax.grid(alpha=.3);ax.legend(fontsize=7)
        ax.set_ylim(0,max(np.log(D)+.3,float(S[0]+1.619*args.T)+.2))
        record={'K':K,'best':winner,'gradient_check_relative_error':grad_error,
                'min_effect_eigenvalue':float(np.linalg.eigvalsh(E)[:,0].min()),
                'completeness_error':float(np.linalg.norm(E.sum(axis=0)-np.eye(D))),
                'outcome_probabilities':prob.tolist(),'effect_volumes':v.tolist(),
                'all_starts':records,'global_optimum_certified':False}
        result['runs'].append(record)
        stored[f'K{int(K)}_effects']=E
        stored[f'K{int(K)}_states']=states
        stored[f'K{int(K)}_probabilities']=prob
        stored[f'K{int(K)}_entropies']=S
    axes[0].set_ylabel('Observational entropy (nats)')
    fig.suptitle('One fixed four-outcome POVM per K (not changed between kicks)')
    fig.tight_layout();fig.savefig(paths[1],dpi=170);plt.close(fig)
    paths[0].write_text(json.dumps(result,indent=2))
    np.savez_compressed(paths[2],**stored)
    print('SAVED',*[str(p) for p in paths],sep='\n',flush=True)


if __name__=='__main__':
    main()
