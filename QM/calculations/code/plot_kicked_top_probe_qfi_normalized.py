"""Normalize saved single-state QFI for a hypothetical final Jy probe."""
import numpy as np
from kicked_top_time_translation_qfi_j10 import OUT, STEM, plt


def main(log_first_five=False, per_kick=False):
    log_first_five = log_first_five or per_kick
    suffix = '_log_first5_per_kick' if per_kick else ('_log_first5' if log_first_five else '')
    target = OUT/f'kicked_top_j10_Jy_probe_qfi_normalized{suffix}.png'
    if target.exists():
        raise FileExistsError('Refusing to overwrite normalized probe-QFI plot')
    with np.load(OUT/f'{STEM}.npz') as data:
        kicks, ks, fq = data['kicks'], data['Ks'], data['qfi']
    assert fq.shape == (len(ks), len(kicks))
    assert np.isfinite(fq).all() and np.all(fq[:, 0] > 0)
    normalized = fq/fq[:, :1]
    assert np.all(normalized[:, 0] == 1)
    if log_first_five:
        mask = kicks <= 5
        kicks, normalized = kicks[mask], normalized[:, mask]
        assert np.all(normalized > 0)
        normalized = np.log(normalized)
        assert np.isfinite(normalized).all() and np.all(normalized[:, 0] == 0)
    if per_kick:
        mask = kicks > 0
        kicks, normalized = kicks[mask], normalized[:, mask]/kicks[mask]
        assert np.isfinite(normalized).all()
    fig, ax = plt.subplots(figsize=(11, 6))
    for k, curve in zip(ks, normalized):
        ax.plot(kicks, curve, lw=1.3, label=f'K={k:g}')
    ax.set(xlabel='Complete kicks n',
           ylabel=(r'$\ln[F_Q(n)/F_Q(0)]/n$ (per kick)' if per_kick else
                   (r'$\ln[F_Q(n)/F_Q(0)]$' if log_first_five else r'$F_Q(n)/F_Q(0)$')),
           title='Full kicked-top evolution; hypothetical Jy probe (j=10)\n'
                 'One coherent initial state: polar angle 1.1, azimuth 0.7',
           xlim=(kicks[0], kicks[-1]))
    ax.grid(alpha=.3)
    ax.legend(ncol=2)
    fig.tight_layout()
    fig.savefig(target, dpi=180)
    plt.close(fig)
    print('Initial QFI:', fq[:, 0])
    print(target)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--log-first-five', action='store_true')
    parser.add_argument('--per-kick', action='store_true')
    args = parser.parse_args()
    main(args.log_first_five, args.per_kick)
