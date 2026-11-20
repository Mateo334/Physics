"""Pure-state time-translation QFI of the kicked top for the fixed generator Jy.

All eight K trajectories use one spin-coherent state. QFI is 4 Var_psi(Jy)
after each COMPLETE Floquet kick, not QFI for changes in kicking strength.
"""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm

from kicked_top_instantaneous_povm_pga import spin,J,THETA,PHI

OUT=Path(__file__).resolve().parents[1]/'results'
STEM='kicked_top_j10_Jy_time_translation_qfi_50kicks'
KS=(.5,1.,2.,4.,6.,8.,10.,12.)
N=50


def main():
    files=[OUT/f'{STEM}.{ext}' for ext in ('png','csv','npz','json')]
    if any(f.exists() for f in files):
        raise FileExistsError('Refusing to overwrite prior QFI output')
    OUT.mkdir(exist_ok=True)
    m,Jy,psi0=spin()
    d=2*J+1
    assert len(psi0)==d and np.linalg.norm(Jy-Jy.conj().T)<1e-13
    rotation=expm(-1j*np.pi*Jy/2)
    assert np.linalg.norm(rotation.conj().T@rotation-np.eye(d))<1e-12
    Jy2=Jy@Jy
    n0_expected=2*J*(1-(np.sin(THETA)*np.sin(PHI))**2)
    curves=np.empty((len(KS),N+1))
    diagnostics=[]
    for index,K in enumerate(KS):
        twist=np.exp(-1j*K*m*m/(2*J))
        state=psi0.copy()
        norm_error=0.
        imag_error=0.
        for n in range(N+1):
            norm_error=max(norm_error,abs(np.vdot(state,state).real-1))
            y=np.vdot(state,Jy@state)
            y2=np.vdot(state,Jy2@state)
            imag_error=max(imag_error,abs(y.imag),abs(y2.imag))
            fq=4*(y2.real-y.real**2)
            assert -1e-10<=fq<=4*J*J+1e-10
            curves[index,n]=fq
            if n<N:
                state=twist*(rotation@state) # U=e^{-i K Jz²/(2j)} e^{-i pi Jy/2}
        assert abs(curves[index,0]-n0_expected)<1e-10
        assert norm_error<1e-11
        diagnostics.append({'K':K,'norm_error_max':norm_error,'moment_imag_error_max':imag_error,
                            'qfi_min':float(curves[index].min()),
                            'qfi_max':float(curves[index].max())})
    fig,ax=plt.subplots(figsize=(11,6))
    kicks=np.arange(N+1)
    for K,curve in zip(KS,curves):
        ax.plot(kicks,curve,lw=1.35,label=f'K={K:g}')
    ax.set(xlabel='Complete kicks n',ylabel=r'$F_Q^{(t)}=4\mathrm{Var}_{\psi_n}(J_y)$',
           title=r'Kicked top: fixed-$J_y$ time-translation QFI ($j=10$, $\hbar=1$)',
           xlim=(0,N))
    ax.grid(alpha=.3);ax.legend(ncol=2,fontsize=9)
    fig.tight_layout();fig.savefig(files[0],dpi=180);plt.close(fig)
    with files[1].open('w',newline='') as out:
        writer=csv.writer(out)
        writer.writerow(['n']+[f'K_{K:g}' for K in KS])
        for n in kicks:
            writer.writerow([n]+[f'{value:.12g}' for value in curves[:,n]])
    np.savez_compressed(files[2],kicks=kicks,Ks=np.array(KS),qfi=curves)
    report={'j':J,'dimension':d,'theta':THETA,'phi':PHI,'hbar':1,
            'generator':'fixed Jy (auxiliary time-translation family, not Floquet Hamiltonian)',
            'floquet_operator':'exp(-i K Jz²/(2j)) exp(-i pi Jy/2)',
            'kicks_inclusive':[0,N],'initial_coherent_expected_qfi':float(n0_expected),
            'qfi_at_selected_kicks':{str(n):{f'{K:g}':float(curves[i,n]) for i,K in enumerate(KS)}
                                      for n in (0,1,5,10,20,30,40,50)},
            'diagnostics':diagnostics}
    files[3].write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    print('Saved:',*[str(p) for p in files],sep='\n')


if __name__=='__main__':
    main()
