"""d=64 torus rotor: Gaussian packet, momentum/energy variance and POVM rate.

Assumed torus convention hbar_eff=2pi/d; free Hamiltonian H_phys=p²/2;
time generator H_phys/hbar_eff. Sample immediately in a kick-free interval.
Run after canonical_instantaneous_obsen_coefficient.py --dimension 64 --outcomes 4.
"""
import json
from pathlib import Path

import numpy as np

OUT=Path(__file__).resolve().parents[1]/'results'
D=64
Q0=np.pi/2
P0=0.
PREFIX='kicked_rotor_gaussian_free_d64_m4'


def entropy(E,psi):
    p=np.einsum('a,iab,b->i',psi.conj(),E,psi).real
    V=np.trace(E,axis1=1,axis2=2).real
    assert np.min(p)>1e-12
    return float(-np.dot(p,np.log(p)-np.log(V)))


def main():
    OUT.mkdir(exist_ok=True)
    data_file=OUT/f'{PREFIX}.npz'
    summary_file=OUT/f'{PREFIX}.json'
    if data_file.exists() or summary_file.exists():
        raise FileExistsError('Refusing to overwrite rotor results')
    coeff_file=OUT/'canonical_obsen_coefficient_d64_m4.npz'
    coeff_summary=json.loads((OUT/'canonical_obsen_coefficient_d64_m4.json').read_text())
    coeff=np.load(coeff_file)
    A=coeff['two_by_two_effects'];alpha=coeff['complement_weights']
    hbar=2*np.pi/D
    sigma_q=np.sqrt(hbar/2)
    q=2*np.pi*np.arange(D)/D
    displacement=(q-Q0+np.pi)%(2*np.pi)-np.pi
    # Periodized continuous Gaussian sampled on an exactly aligned torus grid.
    psi_q=sum(np.exp(-(q-Q0+2*np.pi*shift)**2/(4*sigma_q**2))
              for shift in range(-2,3)).astype(complex)
    psi_q/=np.linalg.norm(psi_q)
    psi=np.fft.fft(psi_q,norm='ortho')
    n=np.fft.fftfreq(D,d=1/D)
    p=hbar*n
    prob=abs(psi)**2
    assert abs(np.sum(prob)-1)<1e-12
    p_mean=float(np.dot(prob,p))
    p2=float(np.dot(prob,p**2))
    p4=float(np.dot(prob,p**4))
    var_p=p2-p_mean**2
    width_p=np.sqrt(var_p)
    q_var=float(np.dot(abs(psi_q)**2,displacement**2))
    width_q=np.sqrt(q_var)
    hbar_inferred=2*width_q*width_p
    energy_unc=np.sqrt(max(p4-p2**2,0))/2
    gaussian_energy_unc=var_p/np.sqrt(2) # p0=0, assuming Gaussian momentum.
    ideal_energy_unc=hbar/(2*np.sqrt(2))
    energy=p**2/(2*hbar) # Generator H_phys/hbar, not the physical kinetic energy.
    mu=float(np.dot(prob,energy))
    tangent=-1j*(energy-mu)*psi
    delta_gen=np.linalg.norm(tangent)
    assert np.isclose(delta_gen,energy_unc/hbar,rtol=1e-11)
    u=psi
    v=-tangent/delta_gen
    frame=np.column_stack([u,v])
    assert np.linalg.norm(frame.conj().T@frame-np.eye(2))<1e-11
    proj=frame@frame.conj().T
    E=np.stack([frame@a@frame.conj().T+weight*(np.eye(D)-proj)
                for a,weight in zip(A,alpha)])
    E=(E+E.conj().transpose(0,2,1))/2
    min_eigs=np.linalg.eigvalsh(E)[:,0]
    completeness=np.linalg.norm(E.sum(axis=0)-np.eye(D))
    assert min(min_eigs)>-1e-10 and completeness<1e-10
    p_i=np.einsum('a,iab,b->i',psi.conj(),E,psi).real
    dotpsi=-1j*energy*psi
    dp_i=2*np.einsum('a,iab,b->i',dotpsi.conj(),E,psi).real
    V=np.trace(E,axis1=1,axis2=2).real
    rate=-np.dot(dp_i,np.log(p_i)-np.log(V))
    coeff_rate=coeff_summary['best_rate']*delta_gen
    assert np.isclose(rate,coeff_rate,rtol=1e-11)
    step=1e-5
    forward=np.exp(-1j*energy*step)*psi
    backward=np.exp(1j*energy*step)*psi
    finite=(entropy(E,forward)-entropy(E,backward))/(2*step)
    assert abs(finite-rate)<2e-6,(finite,rate)
    summary={'d':D,'outcomes':4,'q0':Q0,'p0':P0,'hbar_eff_assumed':hbar,
             'initial_position_width_parameter':sigma_q,'position_variance_measured':q_var,
             'momentum_mean_measured':p_mean,'momentum_variance_measured':var_p,
             'momentum_width_measured':width_p,'hbar_from_measured_widths':hbar_inferred,
             'energy_uncertainty_exact_discrete_phys':energy_unc,
             'energy_uncertainty_gaussian_from_varp':gaussian_energy_unc,
             'energy_uncertainty_ideal_continuous':ideal_energy_unc,
             'generator_uncertainty':delta_gen,'four_outcome_coefficient_attained':coeff_summary['best_rate'],
             'attained_instantaneous_rate_nats_per_free_time':rate,
             'time_central_difference':finite,
             'minimum_effect_eigenvalues':min_eigs.tolist(),
             'completeness_hs_error':completeness,'povm_file':str(data_file),
             'conventions':'torus momentum p_n=hbar_eff*n; norm=ortho FFT; kinetic free interval H_phys=p²/2; generator H_phys/hbar_eff; before the first delta kick',
             'coefficient_is_not_globally_certified':True}
    np.savez_compressed(data_file,effects=E,psi_position=psi_q,psi_momentum=psi,
                        momentum=p,momentum_probabilities=prob,q=q,
                        position_probabilities=abs(psi_q)**2,
                        outcome_probabilities=p_i,outcome_derivatives=dp_i)
    summary_file.write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
