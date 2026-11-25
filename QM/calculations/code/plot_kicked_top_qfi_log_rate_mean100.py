"""Average individual log QFI ratios per elapsed kick, not log of mean QFI."""
import numpy as np
from kicked_top_time_translation_qfi_j10 import OUT, plt
from kicked_top_time_translation_qfi_100_states import STEM


def main():
    stem = 'kicked_top_j10_Jy_probe_qfi_log_rate_first5_mean100'
    png, csv = OUT/f'{stem}.png', OUT/f'{stem}.csv'
    if png.exists() or csv.exists():
        raise FileExistsError('Refusing to overwrite ensemble log-rate results')
    with np.load(OUT/f'{STEM}.npz') as data:
        kicks, ks, fq = data['kicks'], data['Ks'], data['qfi']
    assert fq.shape == (len(ks), len(kicks), 100)
    mask = (kicks >= 1) & (kicks <= 5)
    times = kicks[mask]
    assert np.all(fq[:, :1] > 0) and np.all(fq[:, mask] > 0)
    rates = np.log(fq[:, mask]/fq[:, :1])/times[None, :, None]
    assert np.isfinite(rates).all()
    mean = rates.mean(axis=2)
    sem = rates.std(axis=2, ddof=1)/10
    assert np.allclose(mean, np.sum(rates, axis=2)/100, rtol=0, atol=1e-14)
    fig, ax = plt.subplots(figsize=(11, 6))
    for k, curve in zip(ks, mean):
        ax.plot(times, curve, marker='o', markersize=4, lw=1.3, label=f'K={k:g}')
    ax.set(xlabel='Complete kicks n',
           ylabel=r'$\frac{1}{100}\sum_s\ln[F_{Q,s}(n)/F_{Q,s}(0)]/n$ (per kick)',
           title='Full kicked top; hypothetical Jy probe (j=10)\n'
                 'Mean over 100 sphere-uniform coherent initial states',
           xlim=(1, 5), xticks=times)
    ax.grid(alpha=.3)
    ax.legend(ncol=2)
    fig.tight_layout()
    fig.savefig(png, dpi=180)
    plt.close(fig)
    header = 'n,'+','.join([f'K_{k:g}_mean' for k in ks]+[f'K_{k:g}_sem' for k in ks])
    np.savetxt(csv, np.column_stack((times, mean.T, sem.T)),
               delimiter=',', header=header, comments='')
    print('Checks passed; reused seed-20261122 ensemble, no dynamics recomputation.')
    print(png)
    print(csv)


if __name__ == '__main__':
    main()
