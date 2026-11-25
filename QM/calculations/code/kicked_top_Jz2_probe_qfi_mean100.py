"""Hypothetical exp(-i theta Jz²) probe after complete Floquet kicks."""
import json
import numpy as np
from kicked_top_time_translation_qfi_j10 import OUT, KS, J, spin, expm, plt
from kicked_top_time_translation_qfi_100_states import STEM as SOURCE, SEED


def main():
    stem = 'kicked_top_j10_Jz2_probe_qfi_log_rate_first5_mean100'
    files = [OUT/f'{stem}.{ext}' for ext in ('png', 'csv', 'npz', 'json')]
    if any(p.exists() for p in files):
        raise FileExistsError('Refusing to overwrite Jz² probe results')
    with np.load(OUT/f'{SOURCE}.npz') as saved:
        theta, phi = saved['theta'], saved['phi']
    m, jy, _ = spin()
    north = np.zeros(2*J+1, complex)
    north[-1] = 1
    initial = np.column_stack([np.exp(-1j*p*m)*(expm(-1j*t*jy)@north)
                               for t, p in zip(theta, phi)])
    rotation = expm(-1j*np.pi*jy/2)
    g = m*m
    values = np.empty((len(KS), 6, len(theta)))
    norm_error = 0.
    control_error = 0.
    for a, k in enumerate(KS):
        state = initial.copy()
        twist = np.exp(-1j*k*g/(2*J))
        for n in range(6):
            norm_error = max(norm_error, float(np.max(abs(np.sum(abs(state)**2, axis=0)-1))))
            mean_g = np.sum(g[:, None]*abs(state)**2, axis=0)
            # Center the generator to avoid subtracting two large fourth moments.
            tangent = (g[:, None]-mean_g)*state
            fq = 4*np.sum(abs(tangent)**2, axis=0)
            direct = 4*(np.sum(g[:, None]**2*abs(state)**2, axis=0)-mean_g**2)
            assert np.max(abs(fq-direct)) < 1e-8
            values[a, n] = fq
            if n < 5:
                rotated = rotation@state
                state = twist[:, None]*rotated
                control_error = max(control_error, float(np.max(
                    abs(abs(state)**2-abs(rotated)**2))))
    assert norm_error < 1e-11 and control_error < 1e-12
    assert np.isfinite(values).all() and np.all(values > 0)
    assert np.max(abs(values[:, 0]-values[0, 0])) < 1e-10
    times = np.arange(1, 6)
    rates = np.log(values[:, 1:]/values[:, :1])/times[None, :, None]
    assert np.isfinite(rates).all()
    mean, sem = rates.mean(axis=2), rates.std(axis=2, ddof=1)/np.sqrt(len(theta))
    # First-step Jz² moments cannot depend on K: the final twist commutes with Jz².
    assert np.max(abs(values[:, 1]-values[0, 1])) < 1e-8
    fig, ax = plt.subplots(figsize=(11, 6))
    for k, curve in zip(KS, mean):
        ax.plot(times, curve, marker='o', markersize=4, lw=1.3, label=f'K={k:g}')
    ax.set(xlabel='Complete kicks n',
           ylabel=r'$\frac{1}{100}\sum_s\ln[F_{Q,s}(n)/F_{Q,s}(0)]/n$ (per kick)',
           title='Full kicked top; hypothetical Jz² probe (j=10)\n'
                 'Mean over the same 100 sphere-uniform coherent states',
           xlim=(1, 5), xticks=times)
    ax.grid(alpha=.3)
    ax.legend(ncol=2)
    fig.tight_layout()
    fig.savefig(files[0], dpi=180)
    plt.close(fig)
    header = 'n,'+','.join([f'K_{k:g}_mean' for k in KS]+[f'K_{k:g}_sem' for k in KS])
    np.savetxt(files[1], np.column_stack((times, mean.T, sem.T)),
               delimiter=',', header=header, comments='')
    np.savez_compressed(files[2], qfi=values, rates=rates, mean=mean, sem=sem,
                        theta=theta, phi=phi, Ks=KS, kicks=np.arange(6))
    report = {'generator': 'Jz², hypothetical probe AFTER complete kick; not applied',
              'family': 'exp(-i parameter Jz²) U_K^n |coherent>',
              'j': J, 'count': len(theta), 'seed': SEED,
              'statistic': 'mean_s log(FQ_s(n)/FQ_s(0))/n, n=1..5',
              'norm_error_max': norm_error, 'twist_population_error': control_error,
              'precision': 'complex128', 'source_angles': f'{SOURCE}.npz'}
    files[3].write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    print(files[0])


if __name__ == '__main__':
    main()
