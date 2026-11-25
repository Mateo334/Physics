"""Physical-time QFI for the workspace's two-half-period kicked-top drive."""
import json
import csv

import numpy as np
from kicked_top_time_translation_qfi_j10 import OUT, KS, N, J, spin, expm, plt
from kicked_top_time_translation_qfi_100_states import STEM as SOURCE

STEM = 'kicked_top_j10_physical_time_qfi_50periods_mean100'


def qfi(state, hstate):
    expectation = np.sum(state.conj()*hstate, axis=0)
    assert np.max(abs(expectation.imag)) < 1e-9
    return 4*(np.sum(abs(hstate)**2, axis=0)-expectation.real**2)


def main():
    files = [OUT/f'{STEM}.{ext}' for ext in ('png', 'csv', 'npz', 'json')]
    if any(p.exists() for p in files):
        raise FileExistsError('Refusing to overwrite physical-time QFI outputs')
    with np.load(OUT/f'{SOURCE}.npz') as saved:
        theta, phi = saved['theta'], saved['phi']
        angular_qfi = saved['qfi']
    m, jy, _ = spin()
    north = np.zeros(2*J+1, complex)
    north[-1] = 1
    initial = np.column_stack([np.exp(-1j*p*m)*(expm(-1j*t*jy)@north)
                               for t, p in zip(theta, phi)])
    rotation = expm(-1j*np.pi*jy/2)
    values = np.empty((len(KS), 2*N, len(theta)))
    norm_error = 0.
    for i, k in enumerate(KS):
        state = initial.copy()
        ht = k*m*m/J
        twist = np.exp(-.5j*ht)
        for n in range(N):
            values[i, 2*n] = qfi(state, np.pi*(jy@state))
            state = rotation@state
            values[i, 2*n+1] = qfi(state, ht[:, None]*state)
            state = twist[:, None]*state
            norm_error = max(norm_error, float(np.max(abs(np.sum(abs(state)**2, axis=0)-1))))
    assert norm_error < 1e-11
    assert np.isfinite(values).all() and values.min() >= -1e-8
    assert np.max(abs(values[:, ::2]-np.pi**2*angular_qfi[:, :N])) < 1e-8
    # Energy variance is conserved within each constant-H segment.
    probe = expm(-.17j*np.pi*jy)@initial
    assert np.max(abs(qfi(probe, np.pi*(jy@probe))
                      - qfi(initial, np.pi*(jy@initial)))) < 1e-9
    edges = np.arange(2*N+1)/2
    mean = values.mean(axis=2)
    sem = values.std(axis=2, ddof=1)/np.sqrt(len(theta))
    fig, ax = plt.subplots(figsize=(12, 6))
    for k, curve in zip(KS, mean):
        ax.stairs(curve, edges, baseline=None, lw=1, label=f'K={k:g}')
    ax.set(xlabel='Physical time t (drive period = 1)',
           ylabel=r'$\overline{F_Q^{(t)}}=\overline{4\mathrm{Var}(H(t))}$',
           title='Time QFI: 100 coherent states, j=10\n'
                 r'$H=\pi J_y$ first half-period; $H=KJ_z^2/j$ second half-period',
           xlim=(0, N))
    ax.grid(alpha=.3)
    ax.legend(ncol=2, fontsize=9)
    fig.tight_layout()
    fig.savefig(files[0], dpi=180)
    plt.close(fig)
    with files[1].open('w', newline='') as out:
        writer = csv.writer(out)
        writer.writerow(['t_start', 't_end']+[f'K_{k:g}_mean' for k in KS]
                        +[f'K_{k:g}_sem' for k in KS])
        for a in range(2*N):
            writer.writerow([edges[a], edges[a+1]]+mean[:, a].tolist()+sem[:, a].tolist())
    np.savez_compressed(files[2], edges=edges, Ks=KS, qfi=values,
                        mean=mean, sem=sem, theta=theta, phi=phi)
    report = json.loads((OUT/f'{SOURCE}.json').read_text())
    report.update(generator='physical H(t): pi Jy first half, K Jz²/j second half',
                  period=1, time_interval=[0, N], norm_error_max=norm_error,
                  statistic='mean of individual physical-time QFIs; not mixture QFI',
                  switches='one-sided limits; derivative generally undefined at switches',
                  source_ensemble=f'{SOURCE}.npz')
    report.pop('zero_twist_control_error', None)
    report.pop('initial_mean_qfi', None)
    files[3].write_text(json.dumps(report, indent=2))
    print(f'Checks passed, maximum norm error {norm_error:.3g}')
    print(*files, sep='\n')


if __name__ == '__main__':
    main()
