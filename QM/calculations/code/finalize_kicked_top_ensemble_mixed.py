"""Make refined effects strictly positive and plot fixed-POVM ensemble curves."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import kicked_top_instantaneous_povm_pga as core
from kicked_top_ensemble_mixed_povm import average_rate_gradient,individuals,entropy_trajectory,PREFIX

OUT=Path(__file__).resolve().parents[1]/'results'
EPS=1e-6


def main():
    result={'strict_feasibility_mixture':EPS,'results':[]}
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    _,jy,_=core.spin()
    for k in (3,6):
        saved=np.load(OUT/f'{PREFIX}_refined_K{k}.npz')
        old=np.load(OUT/f'{PREFIX}_K{k}.npz')
        E=(1-EPS)*saved['effects']+EPS*np.eye(core.D)[None,:,:]/core.M
        rho,drho=saved['rho_t0'],saved['drho_t0']
        rate,S,G,*_=average_rate_gradient(E,rho,drho)
        per=individuals(E,rho,drho)
        eigen=np.linalg.eigvalsh(E)[:,0]
        completeness=np.linalg.norm(E.sum(axis=0)-np.eye(core.D))
        residual=core.residual(E,G)
        assert np.min(eigen)>0 and completeness<1e-10
        assert np.allclose(rate,per[0].mean(),atol=1e-11)
        dt=np.linspace(-.09,.09,61)
        curve=entropy_trajectory(E,rho,jy,dt)
        h=1e-6
        fd=(entropy_trajectory(E,rho,jy,[h])[0]-entropy_trajectory(E,rho,jy,[-h])[0])/(2*h)
        assert abs(fd-rate)<1e-4
        file=OUT/f'{PREFIX}_final_K{k}.npz'
        assert not file.exists()
        np.savez_compressed(file,effects=E,rho_initial=saved['rho_initial'],rho_t0=rho,
                            drho_t0=drho,eta=saved['eta'],theta=old['theta'],phi=old['phi'],
                            individual_rates=per[0],individual_entropies=per[1],
                            dt=dt,mean_entropy_curve=curve,refinement_trace=saved['trace'])
        record={'K':k,'rate':rate,'mean_entropy':S,'pure_mean_rate':float(per[0][:10].mean()),
                'mixed_mean_rate':float(per[0][10:].mean()),'min_individual_rate':float(per[0].min()),
                'max_individual_rate':float(per[0].max()),'minimum_eigenvalues':eigen.tolist(),
                'normalization_error':completeness,'projected_gradient_residual':residual,
                'central_difference_rate':fd,'original_90step_rate':float(old['individual_rates'].mean()),
                'refined_unmixed_rate':float(saved['individual_rates'].mean()),
                'file':str(file)}
        result['results'].append(record)
        print(json.dumps(record,indent=2),flush=True)
        color='C0' if k==3 else 'C1'
        initial=json.loads((OUT/f'{PREFIX}_summary.json').read_text())['results'][0 if k==3 else 1]
        axtr=np.array(initial['restarts'][initial['best_restart']]['trace'])
        traces=np.concatenate([axtr[:-1],saved['trace']])
        axes[0].plot(traces,color=color,label=f'K={k}')
        axes[1].plot(core.T0+dt,curve,color=color,label=f'K={k}, slope={rate:.3f}')
    axes[0].set(xlabel='Initial + refinement iteration',ylabel='Mean instantaneous rate (nats/time)')
    axes[1].set(xlabel='Time in rotation interval',ylabel='Mean individual OE (nats)')
    for ax in axes:
        ax.grid(alpha=.3)
        ax.legend()
    fig.tight_layout()
    fig.savefig(OUT/f'{PREFIX}_final_comparison.png',dpi=170)
    plt.close(fig)
    file=OUT/f'{PREFIX}_final_summary.json'
    assert not file.exists()
    file.write_text(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
