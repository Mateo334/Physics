"""Plot t=0..3 OE using the saved best mu=1 interval-objective POVM, unchanged."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from kicked_rotor_K10_postkick_instantaneous import setup,free,kick
from kicked_rotor_fixed_povm_0_to_3_optimize import entropy

RESULTS=Path(__file__).resolve().parents[1]/'results'
SOURCE=RESULTS/'kicked_rotor_K10_d64_fixed_povm_linearity_mu1_50.npz'
TARGET=RESULTS/'kicked_rotor_K10_d64_fixed_povm_linearity_mu1_full_evolution.png'

def main():
    if TARGET.exists():
        raise FileExistsError(f'Refusing to overwrite {TARGET}')
    with np.load(SOURCE) as data:
        F=data['frame'];A=data['two_by_two_effects'];alpha=data['complement_weights']
        saved_S=data['entropy_fine'];saved_t=data['time_fine']
        initial=data['initial_state'];saved_first=data['psi_after_first']
    d=64;q,p,H,psi=setup(d)
    assert np.linalg.norm(initial-psi)<1e-12
    times=[];entropies=[];at_kicks=[]
    n_per_interval=301
    for k in range(3):
        for delta in np.linspace(0,1,n_per_interval):
            times.append(k+delta)
            entropies.append(entropy(A,alpha,F,free(psi,H,delta),d))
        before=entropies[-1]
        psi=kick(free(psi,H,1),q)
        after=entropy(A,alpha,F,psi,d)
        at_kicks.append((before,after))
        if k==0:
            assert np.linalg.norm(psi-saved_first)<1e-11
            assert abs(after-saved_S[0])<1e-10
        if k==1:
            assert abs(before-saved_S[-1])<1e-10
        times.extend([k+1,k+1]);entropies.extend([before,after])
    # The very last point is immediately after the third kick.
    fig,ax=plt.subplots(figsize=(10,4.9))
    ax.plot(times,entropies,color='C0',lw=1.7,label='S_obs with one unchanged POVM')
    ax.axvspan(1,2,color='C1',alpha=.12,label='Optimized interval 1+ to 2-')
    ax.scatter([1,2],[at_kicks[0][1],at_kicks[1][0]],color='C3',zorder=5,
               label='Endpoints used for R')
    for k in (1,2,3):
        ax.axvline(k,color='0.5',ls=':',lw=1)
    ax.set(xlabel='Time (delta kicks at 1, 2, 3)',ylabel='Observational entropy (nats)',
           title='K=10, d=64: one fixed four-outcome POVM, optimized for interval 1+ to 2-')
    ax.grid(alpha=.25);ax.legend(loc='best',fontsize=8)
    fig.tight_layout();fig.savefig(TARGET,dpi=170);plt.close(fig)
    print('Kick before/after OE:',at_kicks)
    print('S(0):',entropies[0], 'R:',at_kicks[1][0]-at_kicks[0][1], 'S(3+):',entropies[-1])
    print('Plot:',TARGET)

if __name__=='__main__':
    main()
