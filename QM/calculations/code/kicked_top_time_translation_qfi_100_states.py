"""Average individual fixed-Jy QFIs over sphere-uniform coherent preparations."""
import csv
import json
import platform

import numpy as np
import scipy
from kicked_top_time_translation_qfi_j10 import OUT, KS, N, J, spin, expm, plt, matplotlib

COUNT = 100
SEED = 20261122
STEM = 'kicked_top_j10_Jy_time_translation_qfi_50kicks_mean100'


def main():
    files = [OUT / f'{STEM}.{ext}' for ext in ('png', 'csv', 'npz', 'json')]
    if any(p.exists() for p in files):
        raise FileExistsError('Refusing to overwrite prior ensemble QFI output')
    OUT.mkdir(exist_ok=True)
    m, jy, _ = spin()
    rng = np.random.default_rng(SEED)
    theta = np.arccos(rng.uniform(-1, 1, COUNT))
    phi = rng.uniform(0, 2*np.pi, COUNT)
    north = np.zeros(2*J+1, dtype=complex)
    north[-1] = 1
    initial = np.column_stack([
        np.exp(-1j*p*m)*(expm(-1j*t*jy) @ north)
        for t, p in zip(theta, phi)
    ])
    assert np.max(np.abs(np.sum(abs(initial)**2, axis=0)-1)) < 1e-12
    assert np.max(np.abs(np.sum(m[:, None]*abs(initial)**2, axis=0)
                         - J*np.cos(theta))) < 1e-11
    rotation = expm(-1j*np.pi*jy/2)
    assert np.linalg.norm(rotation.conj().T @ rotation-np.eye(2*J+1)) < 1e-12
    expected = 2*J*(1-(np.sin(theta)*np.sin(phi))**2)
    curves = np.empty((len(KS), N+1, COUNT))
    norm_error = 0.
    for index, k in enumerate(KS):
        state = initial.copy()
        twist = np.exp(-1j*k*m*m/(2*J))
        for n in range(N+1):
            norm_error = max(norm_error, float(np.max(
                np.abs(np.sum(abs(state)**2, axis=0)-1))))
            ystate = jy @ state
            y = np.sum(state.conj()*ystate, axis=0)
            assert np.max(abs(y.imag)) < 1e-11
            curves[index, n] = 4*(np.sum(abs(ystate)**2, axis=0)-y.real**2)
            if n < N:
                state = twist[:, None]*(rotation @ state)
        assert np.max(abs(curves[index, 0]-expected)) < 1e-10
    assert norm_error < 1e-11
    assert np.isfinite(curves).all()
    assert curves.min() >= -1e-10 and curves.max() <= 4*J*J+1e-10
    # Zero twist is an exactly conserved-Jy control, independent of the plotted K values.
    control = np.linalg.matrix_power(rotation, N) @ initial
    control_y = jy @ control
    control_qfi = 4*(np.sum(abs(control_y)**2, axis=0)
                     - np.sum(control.conj()*control_y, axis=0).real**2)
    assert np.max(abs(control_qfi-expected)) < 1e-9
    mean = curves.mean(axis=2)
    sem = curves.std(axis=2, ddof=1)/np.sqrt(COUNT)
    kicks = np.arange(N+1)
    fig, ax = plt.subplots(figsize=(11, 6))
    for k, curve in zip(KS, mean):
        ax.plot(kicks, curve, lw=1.35, label=f'K={k:g}')
    ax.set(xlabel='Complete kicks n',
           ylabel=r'$\overline{F_Q}=100^{-1}\sum_s4\mathrm{Var}_{\psi_{n,s}}(J_y)$',
           title='Kicked top: fixed-Jy QFI, mean over 100 coherent states (j=10)',
           xlim=(0, N))
    ax.grid(alpha=.3)
    ax.legend(ncol=2, fontsize=9)
    fig.tight_layout()
    fig.savefig(files[0], dpi=180)
    plt.close(fig)
    with files[1].open('w', newline='') as out:
        writer = csv.writer(out)
        writer.writerow(['n']+[f'K_{k:g}_mean' for k in KS]+[f'K_{k:g}_sem' for k in KS])
        for n in kicks:
            writer.writerow([n]+mean[:, n].tolist()+sem[:, n].tolist())
    np.savez_compressed(files[2], kicks=kicks, Ks=np.array(KS), qfi=curves,
                        mean=mean, sem=sem, theta=theta, phi=phi)
    report = {
        'j': J, 'dimension': 2*J+1, 'hbar': 1, 'count': COUNT, 'seed': SEED,
        'sampling': 'uniform cos(theta) in [-1,1], uniform phi in [0,2pi)',
        'theta': theta.tolist(), 'phi': phi.tolist(), 'Ks': KS,
        'kicks_inclusive': [0, N], 'generator': 'fixed Jy encoded AFTER n complete kicks',
        'statistic': 'mean of individual pure-state QFIs, NOT QFI of mixture',
        'norm_error_max': norm_error, 'initial_mean_qfi': float(mean[0, 0]),
        'zero_twist_control_error': float(np.max(abs(control_qfi-expected))),
        'precision': 'complex128', 'assertion_tolerances': '1e-12 to 1e-9',
        'versions': {'python': platform.python_version(), 'numpy': np.__version__,
                     'scipy': scipy.__version__, 'matplotlib': matplotlib.__version__}
    }
    files[3].write_text(json.dumps(report, indent=2))
    print(f'Checks passed; maximum norm error: {norm_error:.3g}')
    print('Saved:', *files, sep='\n')


if __name__ == '__main__':
    main()
