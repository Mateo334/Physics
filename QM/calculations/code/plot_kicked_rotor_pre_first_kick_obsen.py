"""Plot fixed-POVM d=64 rotor OE during the free interval before first kick."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

OUT=Path(__file__).resolve().parents[1]/'results'
PREFIX='kicked_rotor_gaussian_free_d64_m4'


def main():
    filename=OUT/f'{PREFIX}_pre_first_kick.png'
    table=OUT/f'{PREFIX}_pre_first_kick.npz'
    if filename.exists() or table.exists():
        raise FileExistsError('Refusing to overwrite existing pre-first-kick plot/data')
    saved=np.load(OUT/f'{PREFIX}.npz')
    metadata=json.loads((OUT/f'{PREFIX}.json').read_text())
    E=saved['effects'];psi0=saved['psi_momentum'];p=saved['momentum']
    hbar=metadata['hbar_eff_assumed']
    H=p*p/(2*hbar)
    # Clock convention: prepare at t=0, first delta kick at t=1.
    t=np.linspace(0,1,501,endpoint=False)
    states=np.exp(-1j*t[:,None]*H[None,:])*psi0[None,:]
    prob=np.einsum('ta,iab,tb->ti',states.conj(),E,states).real
    V=np.trace(E,axis1=1,axis2=2).real
    assert np.min(prob)>0
    assert np.max(abs(prob.sum(axis=1)-1))<1e-11
    S=-np.sum(prob*(np.log(prob)-np.log(V[None,:])),axis=1)
    slope=metadata['attained_instantaneous_rate_nats_per_free_time']
    tangent=S[0]+slope*t
    assert abs(S[0]-(-np.dot(saved['outcome_probabilities'],np.log(saved['outcome_probabilities'])-np.log(V))))<1e-11
    assert np.min(S)>-1e-11 and np.max(S)<np.log(len(psi0))+1e-11
    fig,ax=plt.subplots(1,2,figsize=(11,4))
    ax[0].plot(t,S,label='Fixed optimized POVM, free evolution')
    ax[0].axvline(1,color='0.4',ls=':',label='First kick at t=1 (not included)')
    ax[0].set(xlabel='Free-evolution time',ylabel='Observational entropy (nats)',title='Before the first kick')
    mask=t<=.06
    ax[1].plot(t[mask],S[mask],label='Actual entropy')
    ax[1].plot(t[mask],tangent[mask],ls='--',label=f'Tangent at 0: {slope:.4f} nats/time')
    ax[1].set(xlabel='Free-evolution time',ylabel='Observational entropy (nats)',title='Initial derivative (zoom)')
    for a in ax:
        a.grid(alpha=.3);a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(filename,dpi=180)
    plt.close(fig)
    np.savez_compressed(table,t=t,entropy=S,tangent=tangent,probabilities=prob)
    print('S(0)=',S[0],'S(last pre-kick point)=',S[-1],
          'S(.05) interpolation=',np.interp(.05,t,S),'slope at 0=',slope)
    print(filename)


if __name__=='__main__':
    main()
