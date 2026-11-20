"""Fixed-phase kicked top: projected-gradient ascent of instantaneous OE rate.

CVXPY/SCS used for projection if installed; otherwise Dykstra alternating
projections computes the SAME Frobenius-norm Euclidean projection. No packages
are installed by this program. Run from workspace root with Python 3.13.
"""
import json
import time
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm

try:
    import cvxpy as cp
except ImportError:
    cp = None

J, D, M = 10, 21, 4
THETA, PHI = 1.1, 0.7
KAPPAS = (3., 6.)
T0, PHASE = 1.1, .1   # After first full kick; 0.1 into second rotation.
SEED, RESTARTS, MAX_ITER = 20261120, 10, 55
OUT = Path(__file__).resolve().parents[1] / 'results'
PREFIX = 'kicked_top_instantaneous_povm_j10'


def herm(a):
    return (a+a.conj().swapaxes(-1,-2))*.5


def spin():
    q = np.arange(-J, J+1, dtype=float)
    raise_op = np.diag(np.sqrt(J*(J+1)-q[:-1]*(q[:-1]+1)), -1)
    jy = (raise_op-raise_op.T)/(2j)
    north = np.zeros(D,complex)
    north[-1] = 1.
    psi0 = np.exp(-1j*PHI*q)*(expm(-1j*THETA*jy)@north)
    assert abs(np.linalg.norm(psi0)-1) < 1e-12
    return q,jy,psi0


def state_at(k, phase, q, jy, psi0):
    # Cycle [n,n+1/2): H_R=pi Jy; [n+1/2,n+1): H_T=k Jz²/j.
    rotation = expm(-1j*np.pi*jy/2)
    twist = np.exp(-1j*k*q*q/(2*J))
    psi = twist*(rotation@psi0)    # one complete cycle, at t=1
    psi = expm(-1j*np.pi*phase*jy)@psi
    return psi


def classical(k, steps=20000):
    x,y,z = np.sin(THETA)*np.cos(PHI),np.sin(THETA)*np.sin(PHI),np.cos(THETA)
    v = np.array([np.cos(THETA)*np.cos(PHI),np.cos(THETA)*np.sin(PHI),-np.sin(THETA)])
    values = []
    total = 0.
    for n in range(steps):
        c,s = np.cos(k*x),np.sin(k*x)
        jac = np.array([[-k*z*s+k*y*c,s,c],[-k*z*c-k*y*s,c,-s],[-1.,0,0]])
        nxt = np.array([z*c+y*s,-z*s+y*c,-x])
        v = jac@v
        v -= np.dot(v,nxt)*nxt
        scale = np.linalg.norm(v)
        assert scale>0
        total += np.log(scale)
        v /= scale
        x,y,z = nxt
        if n+1 in (100,1000,10000,steps):
            values.append((n+1,total/(n+1)))
    return values


def rate_entropy_gradient(effects,rho,drho):
    p = np.einsum('iab,ba->i',effects,rho).real
    dp = np.einsum('iab,ba->i',effects,drho).real
    volumes = np.trace(effects,axis1=1,axis2=2).real
    # No clipping: if probability is too small to evaluate the gradient, a
    # trial step is rejected, rather than silently changing the objective.
    if np.min(p)<=1e-13 or np.min(volumes)<=1e-13:
        return None
    logratio = np.log(p)-np.log(volumes)
    G = (-logratio[:,None,None]*drho[None,:,:]
         -(dp/p)[:,None,None]*rho[None,:,:]
         +(dp/volumes)[:,None,None]*np.eye(D)[None,:,:])
    return -dp@logratio,-p@logratio,herm(G),p,dp


def psd(a):
    eig,v = np.linalg.eigh(herm(a))
    return (v*np.maximum(eig,0)[None,:])@v.conj().T


def affine(a):
    return a + (np.eye(D)-a.sum(axis=0))[None,:,:]/M


def dykstra_projection(target,tol=1e-7,maxsteps=1600):
    """Euclidean projection onto intersection of PSD cone and sum_i E_i=I."""
    x=herm(target.copy()); p=np.zeros_like(x);q=np.zeros_like(x)
    for step in range(maxsteps):
        y=np.stack([psd(z) for z in x+p])
        p=x+p-y
        next_x=affine(y+q)
        q=y+q-next_x
        x=herm(next_x)
        if step%10==9 or step==maxsteps-1:
            if np.min(np.linalg.eigvalsh(x)) >= -tol and np.linalg.norm(x-y)<tol:
                return x,step+1
    raise RuntimeError(f'Dykstra failed at {maxsteps} sweeps: min eig {np.min(np.linalg.eigvalsh(x))}')


def cvxpy_projection(target):
    if cp is None:
        raise RuntimeError('CVXPY not installed')
    E=[cp.Variable((D,D),hermitian=True) for _ in range(M)]
    problem=cp.Problem(cp.Minimize(sum(cp.sum_squares(cp.abs(E[i]-target[i]))
                                        for i in range(M))),
                       [e >> 0 for e in E]+[sum(E)==np.eye(D)])
    problem.solve(solver=cp.SCS,eps=1e-7,max_iters=30000,verbose=False)
    if problem.status not in ('optimal','optimal_inaccurate'):
        raise RuntimeError(f'CVXPY: {problem.status}')
    result=herm(np.stack([e.value for e in E]))
    if (np.min(np.linalg.eigvalsh(result)) < -1e-6
            or np.linalg.norm(result.sum(axis=0)-np.eye(D))>1e-6):
        raise RuntimeError('CVXPY inaccurate projection')
    return result,problem.solver_stats.num_iters


def project(target):
    return cvxpy_projection(target) if cp is not None else dykstra_projection(target)


def random_povm(rng):
    A=rng.standard_normal((M,D,D))+1j*rng.standard_normal((M,D,D))
    B=A.conj().transpose(0,2,1)@A
    eig,V=np.linalg.eigh(B.sum(axis=0))
    inv=(V*(eig**(-.5))[None,:])@V.conj().T
    E=herm(inv[None]@B@inv[None])
    # Strictly interior start: numerical derivatives remain regular.
    return .97*E+.03*np.eye(D)[None,:,:]/M


def lift_four_to_ten(effects):
    """Proportionally split one outcome; retain exact entropy and rate."""
    assert effects.shape==(4,D,D)
    return np.concatenate([np.repeat(effects[:1]/7,7,axis=0),effects[1:]])


def residual(E,G,alpha=0.05):
    projected,_=project(E+alpha*G)
    return np.linalg.norm(projected-E)/alpha


def gradient_check(rho,drho,rng,evaluator=None):
    if evaluator is None:
        evaluator=rate_entropy_gradient
    E=random_povm(rng)
    DIRECTION=herm(rng.standard_normal((M,D,D))+1j*rng.standard_normal((M,D,D)))
    DIRECTION/=np.linalg.norm(DIRECTION)
    R,S,G,*_=evaluator(E,rho,drho)
    h=1e-5
    finite=(evaluator(E+h*DIRECTION,rho,drho)[0]
           -evaluator(E-h*DIRECTION,rho,drho)[0])/(2*h)
    analytic=np.einsum('iab,iba->',G,DIRECTION).real
    error=abs(finite-analytic)/max(abs(finite),1e-12)
    assert error<1e-6,(finite,analytic,error)
    return error


def optimize(rho,drho,rng,evaluator=None):
    if evaluator is None:
        evaluator=rate_entropy_gradient
    records=[]; best=None
    for restart in range(RESTARTS):
        E=random_povm(rng)
        history=[]; evals=0; status='iteration limit'
        for iteration in range(MAX_ITER):
            R,S,G,p,dp=evaluator(E,rho,drho)
            history.append(float(R))
            # Bound the raw trial displacement: rare outcomes can produce
            # large true gradients and stall Euclidean cone projection.
            eta=min(0.4,1./max(np.linalg.norm(G),1e-12))
            for backtrack in range(18):
                trial,steps=project(E+eta*G)
                evals+=1
                result=evaluator(trial,rho,drho)
                if result is not None and result[0] > R+1e-11:
                    break
                eta*=.5
            else:
                status='no improving step';break
            improvement=result[0]-R
            E=trial
            if improvement < 2e-7:
                res=residual(E,evaluator(E,rho,drho)[2])
                if res<3e-3:
                    status='small improvement and residual';break
        final=evaluator(E,rho,drho)
        res=residual(E,final[2])
        # All starts retained, including imperfectly converged ones.
        record={'restart':restart,'rate':float(final[0]),'entropy':float(final[1]),
                'residual_alpha_0.05':float(res),'iterations':len(history),
                'projections':evals,'status':status,'trace':history+[float(final[0])]}
        records.append(record)
        print('restart',restart,'rate',f'{final[0]:.9g}','residual',f'{res:.3g}',status,flush=True)
        if best is None or record['rate'] > best[0]['rate']:
            best=(record,E.copy())
    return records,best


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(SEED)
    q,jy,psi0=spin()
    information={'seed':SEED,'j':J,'d':D,'m':M,'t0':T0,'rotation_phase':PHASE,
                 'rotation_hamiltonian':'pi Jy', 'twist_hamiltonian':'K Jz^2/j',
                 'projection':'CVXPY/SCS' if cp is not None else 'Dykstra (CVXPY unavailable)',
                 'numpy':np.__version__,'restarts':RESTARTS,'max_iterations':MAX_ITER,
                 'results':[]}
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    start=time.monotonic()
    for k in KAPPAS:
        classical_values=classical(k)
        exponent=classical_values[-1][1]
        tE=np.log(J)/exponent
        psi=state_at(k,PHASE,q,jy,psi0)
        rho=np.outer(psi,psi.conj())
        H=np.pi*jy
        drho=herm(-1j*(H@rho-rho@H))
        assert abs(np.trace(rho)-1)<1e-11
        assert np.min(np.linalg.eigvalsh(rho))>-1e-10
        assert np.linalg.norm(drho)>1e-8
        gc=gradient_check(rho,drho,rng)
        records,(winner,E)=optimize(rho,drho,rng)
        R,S,_,p,dp=rate_entropy_gradient(E,rho,drho)
        min_eig=np.linalg.eigvalsh(E)[:,0]
        normerr=np.linalg.norm(E.sum(axis=0)-np.eye(D))
        assert min(min_eig)>-5e-8 and normerr<5e-8
        # Keep the optimized measurement fixed along the rotation interval.
        dt=np.linspace(-.09,.09,61)
        entropy=[]
        for s in dt:
            state=expm(-1j*np.pi*s*jy)@psi
            matrix=np.outer(state,state.conj())
            entropy.append(rate_entropy_gradient(E,matrix,herm(-1j*(H@matrix-matrix@H)))[1])
        entropy=np.array(entropy)
        h=1e-6
        def local_entropy(s):
            state=expm(-1j*np.pi*s*jy)@psi
            matrix=np.outer(state,state.conj())
            return rate_entropy_gradient(E,matrix,herm(-1j*(H@matrix-matrix@H)))[1]
        central=(local_entropy(h)-local_entropy(-h))/(2*h)
        assert abs(central-R)<1e-3,(central,R)
        filename=OUT/f'{PREFIX}_K{int(k)}.npz'
        np.savez_compressed(filename,effects=E,rho=rho,drho=drho,probabilities=p,
                            probability_derivatives=dp,dt=dt,entropy=entropy,
                            rates=np.array([r['rate'] for r in records]))
        summary={'k':k,'classical_exponent':exponent,'exponent_convergence':classical_values,
                 'ehrenfest_estimate':tE,'gradient_test_relative_error':gc,
                 'best_rate':float(R),'best_entropy':float(S),'min_eigenvalues':min_eig.tolist(),
                 'normalization_error':float(normerr),'central_difference_rate':float(central),
                 'best_restart':winner['restart'],'restarts':records,'npz':str(filename)}
        information['results'].append(summary)
        print('RESULT',json.dumps({key:value for key,value in summary.items() if key!='restarts'},indent=2),flush=True)
        for rec in records:
            axes[0].plot(rec['trace'],color=('C0' if k==3 else 'C1'),alpha=.23)
        axes[0].plot(winner['trace'],color=('C0' if k==3 else 'C1'),lw=2,label=f'K={int(k)}, best')
        axes[1].plot(T0+dt,entropy,label=f'K={int(k)}, slope={R:.3f}')
        axes[1].plot([T0],[S],marker='o',color=('C0' if k==3 else 'C1'))
    axes[0].set(xlabel='Projected-gradient iteration',ylabel='Instantaneous rate (nats/time)')
    axes[1].set(xlabel='Time (rotation interval)',ylabel='Observational entropy (nats)')
    for ax in axes:
        ax.grid(alpha=.3);ax.legend()
    fig.tight_layout()
    fig.savefig(OUT/f'{PREFIX}_comparison.png',dpi=170)
    plt.close(fig)
    information['elapsed_seconds']=time.monotonic()-start
    (OUT/f'{PREFIX}_summary.json').write_text(json.dumps(information,indent=2))
    print('seconds',information['elapsed_seconds'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outcomes',type=int,default=4)
    args=parser.parse_args()
    if args.outcomes<2:
        parser.error('--outcomes must be at least two')
    M=args.outcomes
    if M!=4:
        PREFIX=f'kicked_top_instantaneous_povm_j10_m{M}'
    main()
