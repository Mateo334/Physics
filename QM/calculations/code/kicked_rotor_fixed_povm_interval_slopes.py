"""Secant and least-squares OE slopes for the same POVM across K=10 kicks."""
import json
from pathlib import Path

import numpy as np

from kicked_rotor_K10_postkick_instantaneous import OUT,PREFIX,HBAR,entropy,free


def main():
    target=OUT/f'{PREFIX}_fixed_interval_slopes.json'
    if target.exists():
        raise FileExistsError('Refusing to overwrite interval-slope results')
    archive=np.load(OUT/f'{PREFIX}.npz')
    answer={'slope_definition':'endpoint entropy difference per unit free interval; also sampled OLS fit',
            'same_povm_retained_across_kick':True,
            'instant_optimized_povms_not_interval_slope_optimized':True,'results':[]}
    for d in (64,512):
        state1=archive[f'd{d}_first_state'];state2=archive[f'd{d}_second_state']
        p=archive[f'd{d}_momentum'];H=p*p/(2*HBAR)
        A=archive[f'd{d}_two_by_two_effects'];alpha=archive[f'd{d}_complement_weights']
        table={'d':d,'measurements':[]}
        for chosen in (1,2):
            frame=archive[f'd{d}_{"first" if chosen==1 else "second"}_frame']
            S1=entropy(frame,A,alpha,state1,d)
            S2minus=entropy(frame,A,alpha,free(state1,H,1),d)
            S2plus=entropy(frame,A,alpha,state2,d)
            S3minus=entropy(frame,A,alpha,free(state2,H,1),d)
            t=np.linspace(0,1,401)
            curve1=np.array([entropy(frame,A,alpha,free(state1,H,float(s)),d) for s in t])
            curve2=np.array([entropy(frame,A,alpha,free(state2,H,float(s)),d) for s in t])
            assert np.allclose([curve1[0],curve1[-1],curve2[0],curve2[-1]],
                               [S1,S2minus,S2plus,S3minus],atol=1e-10)
            entry={'povm_optimized_instantaneously_at':f'{chosen}+',
                   'S_1_plus':S1,'S_2_minus':S2minus,'S_2_plus':S2plus,'S_3_minus':S3minus,
                   'first_free_interval_secant':S2minus-S1,
                   'second_free_interval_secant':S3minus-S2plus,
                   'first_free_interval_OLS':np.polyfit(t,curve1,1)[0],
                   'second_free_interval_OLS':np.polyfit(t,curve2,1)[0],
                   'second_kick_entropy_jump':S2plus-S2minus,
                   'two_kick_window_secant_including_jump':(S3minus-S1)/2}
            table['measurements'].append(entry)
        answer['results'].append(table)
        print(json.dumps(table,indent=2))
    target.write_text(json.dumps(answer,indent=2))


if __name__=='__main__':
    main()
