"""Continue the saved 20-preparation common-POVM gradient runs without overwrite."""
import json
from pathlib import Path

import numpy as np

import kicked_top_instantaneous_povm_pga as core
from kicked_top_ensemble_mixed_povm import average_rate_gradient, individuals, entropy_trajectory, PREFIX

OUT=Path(__file__).resolve().parents[1]/'results'


def refine(E,rho,drho,max_iterations=550):
    trace=[]
    for n in range(max_iterations):
        R,S,G,*_=average_rate_gradient(E,rho,drho)
        trace.append(R)
        eta=min(.4,3./max(np.linalg.norm(G),1e-12))
        for _ in range(22):
            trial,_=core.project(E+eta*G)
            eval_trial=average_rate_gradient(trial,rho,drho)
            if eval_trial is not None and eval_trial[0]>R+1e-11:
                break
            eta*=.5
        else:
            return E,trace,'no improving step'
        E=trial
        if (n+1)%50==0:
            resid=core.residual(E,average_rate_gradient(E,rho,drho)[2])
            print(n+1,'rate',round(eval_trial[0],8),'residual',round(resid,7),flush=True)
        if (n+1)%10==0 and eval_trial[0]-R<1e-7:
            resid=core.residual(E,average_rate_gradient(E,rho,drho)[2])
            if resid<.015:
                return E,trace,'small improvement and residual'
    return E,trace,'iteration limit'


def main():
    summary={}
    for k in (3,6):
        if (OUT/f'{PREFIX}_refined_K{k}.npz').exists():
            print('Existing refined result preserved, skipping K',k,flush=True)
            continue
        saved=np.load(OUT/f'{PREFIX}_K{k}.npz')
        E=saved['effects'];rho=saved['rho_t0'];drho=saved['drho_t0']
        before=average_rate_gradient(E,rho,drho)[0]
        E,trace,status=refine(E,rho,drho)
        rate,S,G,*_=average_rate_gradient(E,rho,drho)
        err=core.residual(E,G)
        eig=np.linalg.eigvalsh(E)[:,0]
        completeness=np.linalg.norm(E.sum(axis=0)-np.eye(core.D))
        print('feasibility before saving',k,'before',before,'after',rate,
              'min eigenvalue',min(eig),'completeness',completeness,flush=True)
        assert rate>=before and min(eig)>-1e-7 and completeness<1e-8
        per=individuals(E,rho,drho)
        h=1e-6
        z,jy,_=core.spin()
        fd=(entropy_trajectory(E,rho,jy,[h])[0]-entropy_trajectory(E,rho,jy,[-h])[0])/(2*h)
        assert abs(fd-rate)<1e-4
        file=OUT/f'{PREFIX}_refined_K{k}.npz'
        assert not file.exists()
        np.savez_compressed(file,effects=E,rho_t0=rho,drho_t0=drho,
                            rho_initial=saved['rho_initial'],eta=saved['eta'],
                            individual_rates=per[0],individual_entropy=per[1],
                            trace=np.array(trace+[rate]))
        summary[str(k)]={'before':before,'rate':rate,'mean_entropy':S,
                         'pure_mean_rate':float(per[0][:10].mean()),
                         'mixed_mean_rate':float(per[0][10:].mean()),
                         'residual_alpha_0.05':err,'iterations':len(trace),
                         'status':status,'minimum_eigenvalues':eig.tolist(),
                         'normalization_error':completeness,'finite_difference':fd,
                         'effects_file':str(file)}
        print(k,json.dumps(summary[str(k)],indent=2),flush=True)
    path=OUT/f'{PREFIX}_refined_summary.json'
    assert not path.exists()
    path.write_text(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
