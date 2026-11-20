"""Full-POVM pure-state instantaneous OE coefficient via exact 2D reduction.

For d-dimensional pure rho and norm-one tangent, twirl every effect under all
unitaries on the orthogonal complement of span{|psi>,|v>}. The objective and
feasibility stay unchanged: E_i=A_i(2x2) direct-sum alpha_i I_(d-2).
Run: py -3.13 calculations/code/canonical_instantaneous_obsen_coefficient.py
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

OUT=Path(__file__).resolve().parents[1]/'results'


def herm(x):
    return (x+x.conj().swapaxes(-1,-2))*.5


def simplex(v):
    x=np.sort(v)[::-1]
    rho=np.flatnonzero(x*np.arange(1,len(x)+1)>np.cumsum(x)-1)
    threshold=(np.sum(x[:rho[-1]+1])-1)/(rho[-1]+1)
    return np.maximum(v-threshold,0)


def psd(a):
    e,V=np.linalg.eigh(herm(a))
    return (V*np.maximum(e,0)[None,:])@V.conj().T


def project_A(y,tol=1e-10,max_sweeps=5000):
    """Hilbert--Schmidt projection onto m two-by-two PSD effects summing to I2."""
    m=len(y);x=herm(y.copy());p=np.zeros_like(x);q=np.zeros_like(x)
    for n in range(max_sweeps):
        z=np.stack([psd(a) for a in x+p]);p=x+p-z
        next_x=z+q+(np.eye(2)-(z+q).sum(axis=0))[None,:,:]/m
        q=z+q-next_x;x=herm(next_x)
        if n%10==9 and np.min(np.linalg.eigvalsh(x))>=-tol and np.linalg.norm(x-z)<tol:
            return x
    raise RuntimeError('2x2 POVM projection failed')


def objective(A,alpha,d):
    m=len(A);p=A[:,0,0].real; dp=-2*A[:,0,1].real
    V=np.trace(A,axis1=1,axis2=2).real+(d-2)*alpha
    if min(p.min(),V.min())<=1e-14:
        return None
    logs=np.log(p)-np.log(V)
    rate=-dp@logs
    rho=np.diag([1.,0.]);drho=np.array([[0.,-1.],[-1.,0.]])
    G=-logs[:,None,None]*drho[None,:,:]-(dp/p)[:,None,None]*rho[None,:,:]
    G+=(dp/V)[:,None,None]*np.eye(2)[None,:,:]
    galpha=dp/V  # Full HS metric: ||delta(alpha I_n)||²=n delta(alpha)².
    return float(rate),herm(G),galpha,p,dp,V


def random_povm(m,rng):
    Z=rng.normal(size=(m,2,2))+1j*rng.normal(size=(m,2,2))
    B=Z.conj().transpose(0,2,1)@Z
    e,U=np.linalg.eigh(B.sum(axis=0))
    inv=(U*(e**(-.5))[None,:])@U.conj().T
    A=herm(inv[None]@B@inv[None])
    return .98*A+.02*np.eye(2)[None,:,:]/m, rng.dirichlet(np.ones(m))


def binary_seed(d,m):
    fun=lambda p:2*np.sqrt(p*(1-p))*np.log((d-1)*p/(1-p))
    result=minimize_scalar(lambda p:-fun(p),bounds=(1/d,1-1e-12),
                           method='bounded',options={'xatol':1e-14})
    p=result.x
    u=np.array([np.sqrt(p),np.sqrt(1-p)])
    P=np.outer(u,u)
    A=np.concatenate([P[None,:,:],np.repeat(((np.eye(2)-P)/(m-1))[None,:,:],m-1,axis=0)])
    alpha=np.array([0]+[1/(m-1)]*(m-1))
    return A.astype(complex),alpha,float(-result.fun),float(p)


def residual(A,alpha,G,g_alpha,d,h=.03):
    X=project_A(A+h*G)
    a=simplex(alpha+h*g_alpha)
    return float(np.sqrt(np.linalg.norm((X-A)/h)**2
                         +(d-2)*np.linalg.norm((a-alpha)/h)**2))


def run(d,m,restarts,iterations,seed):
    rng=np.random.default_rng(seed)
    baselineA,baselineAlpha,baselineRate,pstar=binary_seed(d,m)
    results=[];best=None
    for start in range(restarts):
        A,alpha=(baselineA.copy(),baselineAlpha.copy()) if start==0 else random_povm(m,rng)
        trace=[];status='iteration limit'
        for step in range(iterations):
            obj=objective(A,alpha,d)
            if obj is None:
                status='invalid probabilities';break
            R,G,g,*_=obj;trace.append(R)
            # Full-space gradient: norm in metric ||dA||²+n||dalpha||².
            norm=np.sqrt(np.linalg.norm(G)**2+(d-2)*np.linalg.norm(g)**2)
            eta=min(.4,2/max(norm,1e-12))
            for _ in range(25):
                X=project_A(A+eta*G)
                a=simplex(alpha+eta*g)
                trial=objective(X,a,d)
                if trial is not None and trial[0]>R+1e-12:
                    break
                eta*=.5
            else:
                status='no improving step';break
            A,alpha=X,a
            if step%200==199:
                print('d',d,'restart',start,'step',step+1,'rate',trial[0],flush=True)
        final=objective(A,alpha,d)
        if final is None:
            raise RuntimeError('Invalid final POVM')
        R,G,g,p,dp,V=final
        assert abs(np.sum(alpha)-1)<1e-10
        assert np.linalg.norm(A.sum(axis=0)-np.eye(2))<1e-9
        assert np.min(np.linalg.eigvalsh(A))>-1e-8
        res=residual(A,alpha,G,g,d)
        item={'restart':start,'rate':R,'residual_h0.03':res,'iterations':len(trace),
              'status':status,'trace':trace+[R]}
        results.append(item)
        print('d',d,'restart',start,'rate',R,'residual',res,status,flush=True)
        if best is None or R>best[0]['rate']:
            best=(item,A.copy(),alpha.copy())
    row,A,alpha=best
    return {'d':d,'m':m,'best_rate':row['rate'],'best_residual':row['residual_h0.03'],
            'binary_lower_bound':baselineRate,'binary_optimal_p':pstar,
            'restarts':results,'seed':seed},A,alpha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dimension',type=int,default=1024)
    parser.add_argument('--outcomes',type=int,default=10)
    parser.add_argument('--restarts',type=int,default=16)
    parser.add_argument('--iterations',type=int,default=1500)
    parser.add_argument('--seed',type=int,default=20261125)
    args=parser.parse_args()
    assert args.dimension>=3 and args.outcomes>=2
    OUT.mkdir(exist_ok=True)
    prefix=f'canonical_obsen_coefficient_d{args.dimension}_m{args.outcomes}'
    if any((OUT/f'{prefix}{suffix}').exists() for suffix in ['.json','.npz']):
        parser.error('Results already exist; refusing to overwrite')
    summary,A,alpha=run(args.dimension,args.outcomes,args.restarts,args.iterations,args.seed)
    # Directional gradient tests for two distinct interior configurations.
    rng=np.random.default_rng(args.seed+1)
    tests=[]
    for _ in range(2):
        X,a=random_povm(args.outcomes,rng)
        DX=herm(rng.normal(size=X.shape)+1j*rng.normal(size=X.shape))
        da=rng.normal(size=a.shape)
        DX/=np.linalg.norm(DX);da/=np.linalg.norm(da)
        f,G,g,*_=objective(X,a,args.dimension)
        analytic=np.einsum('iab,iba->',G,DX).real+(args.dimension-2)*np.dot(g,da)
        h=1e-6
        numerical=(objective(X+h*DX,a+h*da,args.dimension)[0]
                   -objective(X-h*DX,a-h*da,args.dimension)[0])/(2*h)
        err=abs(analytic-numerical)/max(abs(numerical),1e-12)
        assert err<1e-6
        tests.append(err)
    summary['gradient_checks_relative_error']=tests
    summary['finite_HS_projection']='A 2x2 POVM Euclidean projection plus scalar simplex'
    summary['rho']='|0><0|; dotrho=-(|0><1|+|1><0|); DeltaH=1'
    (OUT/f'{prefix}.json').write_text(json.dumps(summary,indent=2))
    np.savez_compressed(OUT/f'{prefix}.npz',two_by_two_effects=A,complement_weights=alpha,
                        complement_dimension=args.dimension-2)
    print('RESULT',json.dumps({k:v for k,v in summary.items() if k!='restarts'},indent=2),flush=True)


if __name__=='__main__':
    main()
