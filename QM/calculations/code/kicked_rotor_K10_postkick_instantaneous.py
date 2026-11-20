"""K=10 rotor post-kick instantaneous OE: d=64 torus and d=512 cutoff.

The packet is prepared at t=0, free evolution precedes first kick at t=1,
then delta kicks at t=1,2,3. Evaluate at t=1+ and t=2+, where H=p²/2.
Each instant gets its own 4-outcome optimized POVM; distinguish from one fixed
POVM optimized for growth across the entire two-kick window.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

OUT=Path(__file__).resolve().parents[1]/'results'
PREFIX='kicked_rotor_K10_postkick_d64_d512_m4'
K=10.
HBAR=2*np.pi/64
Q0=np.pi/2
SIGMA=np.sqrt(HBAR/2)


def setup(d):
    q=2*np.pi*np.arange(d)/d
    p=np.fft.fftfreq(d,d=1/d)*HBAR
    H=p*p/(2*HBAR)
    wave=sum(np.exp(-(q-Q0+2*np.pi*n)**2/(4*SIGMA**2)) for n in range(-2,3))
    wave=np.asarray(wave,dtype=complex)/np.linalg.norm(wave)
    psi=np.fft.fft(wave,norm='ortho')
    assert abs(np.linalg.norm(psi)-1)<1e-12
    return q,p,H,psi


def free(psi,H,t):
    return np.exp(-1j*H*t)*psi


def kick(psi,q):
    return np.fft.fft(np.exp(-1j*K*np.cos(q)/HBAR)*np.fft.ifft(psi,norm='ortho'),norm='ortho')


def frame_and_uncertainty(psi,H,p):
    P=abs(psi)**2
    p1=float(np.dot(P,p));p2=float(np.dot(P,p*p));p4=float(np.dot(P,p**4))
    mu=float(np.dot(P,H))
    dh=np.linalg.norm((H-mu)*psi)
    assert abs(dh**2-(p4-p2**2)/(4*HBAR**2))<1e-9
    v=1j*(H-mu)*psi/dh # Canonical dotrho=-dh(|1><0|+|0><1|).
    frame=np.stack((psi,v),axis=1)
    assert np.linalg.norm(frame.conj().T@frame-np.eye(2))<1e-10
    near_edge=float(P[np.abs(p)>.8*np.max(abs(p))].sum())
    return frame,{'momentum_mean':p1,'momentum_variance':p2-p1*p1,
                  'physical_energy_uncertainty':dh*HBAR,
                  'generator_uncertainty':dh,'edge_probability_outer_20_percent':near_edge}


def entropy(frame,A,alpha,psi,d):
    z=frame.conj().T@psi
    zweight=float(np.vdot(z,z).real)
    probabilities=np.einsum('a,iab,b->i',z.conj(),A,z).real+alpha*(1-zweight)
    volumes=np.trace(A,axis1=1,axis2=2).real+(d-2)*alpha
    if np.min(probabilities)<=0:
        raise ValueError('Nonpositive probability along free trajectory')
    assert abs(probabilities.sum()-1)<1e-10
    return float(-np.dot(probabilities,np.log(probabilities)-np.log(volumes)))


def main():
    OUT.mkdir(exist_ok=True)
    files=[OUT/f'{PREFIX}.png',OUT/f'{PREFIX}.json',OUT/f'{PREFIX}.npz']
    if any(file.exists() for file in files):
        raise FileExistsError('Refusing to overwrite post-kick results')
    summary={'kick_strength':K,'hbar_eff_fixed':HBAR,'q0':Q0,'p0':0,
             'position_gaussian_sigma':SIGMA,'time_schedule':'prepare t=0; free [0,1); delta kicks t=1,2,3; free H=p²/2 between kicks',
             'outcomes':4,'results':[]}
    all_data={}
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    for ax,d in zip(axes,(64,512)):
        coeff=json.loads((OUT/f'canonical_obsen_coefficient_d{d}_m4.json').read_text())
        compact=np.load(OUT/f'canonical_obsen_coefficient_d{d}_m4.npz')
        A=compact['two_by_two_effects'];alpha=compact['complement_weights']
        q,p,H,initial=setup(d)
        first=kick(free(initial,H,1),q)
        second=kick(free(first,H,1),q)
        after=[first,second]
        frames=[];entries=[]
        for index,state in enumerate(after,1):
            assert abs(np.linalg.norm(state)-1)<1e-11
            frame,moments=frame_and_uncertainty(state,H,p)
            frames.append(frame)
            coefficient=coeff['best_rate']
            tangent=-1j*(H-np.vdot(state,H*state).real)*state
            dotp=-2*moments['generator_uncertainty']*A[:,0,1].real
            vols=np.trace(A,axis1=1,axis2=2).real+(d-2)*alpha
            probabilities=A[:,0,0].real
            direct=-np.dot(dotp,np.log(probabilities)-np.log(vols))
            rate=coefficient*moments['generator_uncertainty']
            assert abs(direct-rate)<1e-9
            h=1e-7
            fd=(entropy(frame,A,alpha,free(state,H,h),d)
                -entropy(frame,A,alpha,free(state,H,-h),d))/(2*h)
            assert abs(fd-rate)<2e-4,(d,index,fd,rate)
            row={'after_kick':index,**moments,'candidate_max_instant_rate':rate,
                 'local_time_difference':fd,'S_at_left_endpoint':entropy(frame,A,alpha,state,d)}
            entries.append(row)
            print('d',d,'kick',index,json.dumps(row),flush=True)
        first_time=np.linspace(1,2,151,endpoint=False)
        second_time=np.linspace(2,3,151,endpoint=False)
        first_curve=np.array([entropy(frames[0],A,alpha,free(first,H,t-1),d) for t in first_time])
        second_curve=np.array([entropy(frames[1],A,alpha,free(second,H,t-2),d) for t in second_time])
        fixed_first_on_second=np.array([entropy(frames[0],A,alpha,free(second,H,t-2),d)
                                         for t in second_time])
        ax.plot(first_time,first_curve,color='C0',label='POVM optimized at 1+: held fixed (1,2)')
        ax.plot(second_time,second_curve,color='C1',label='POVM optimized at 2+: held fixed (2,3)')
        ax.plot(second_time,fixed_first_on_second,color='C0',ls='--',alpha=.7,
                label='1+ POVM retained after second kick')
        ax.axvline(2,color='gray',ls=':',label='Second delta kick')
        ax.axvline(3,color='gray',ls=':')
        ax.set(xlabel='Time (first, second and third kicks at 1, 2, 3)',
               ylabel='Observational entropy (nats)',title=f'd={d}, K=10')
        ax.grid(alpha=.3);ax.legend(fontsize=7)
        summary['results'].append({'d':d,'coefficient_attained':coefficient,
                                   'coefficient_residual':coeff['best_residual'],
                                   'after_kicks':entries})
        for name,val in [('initial_state',initial),('first_state',first),('second_state',second),
                         ('momentum',p),('first_frame',frames[0]),('second_frame',frames[1]),
                         ('first_time',first_time),('second_time',second_time),
                         ('first_entropy',first_curve),('second_entropy',second_curve),
                         ('fixed_first_on_second',fixed_first_on_second),('two_by_two_effects',A),
                         ('complement_weights',alpha)]:
            all_data[f'd{d}_{name}']=val
    fig.tight_layout()
    fig.savefig(files[0],dpi=175)
    plt.close(fig)
    files[1].write_text(json.dumps(summary,indent=2))
    np.savez_compressed(files[2],**all_data)


if __name__=='__main__':
    main()
