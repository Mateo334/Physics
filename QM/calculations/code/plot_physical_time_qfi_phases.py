"""Readable phase-separated view of the saved physical-time QFI calculation."""
import numpy as np
from kicked_top_time_translation_qfi_j10 import OUT, plt
from kicked_top_physical_time_qfi_100_states import STEM


def main():
    target = OUT/f'{STEM}_phases.png'
    if target.exists():
        raise FileExistsError('Refusing to overwrite phase-separated figure')
    with np.load(OUT/f'{STEM}.npz') as data:
        mean, ks, edges = data['mean'], data['Ks'], data['edges']
    assert mean.shape == (len(ks), len(edges)-1)
    assert mean.shape[1] % 2 == 0 and np.isfinite(mean).all()
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    for phase, (ax, title) in enumerate(zip(axes, (
        r'Rotation phase: $H=\pi J_y$',
        r'Twist phase: $H=KJ_z^2/j$',
    ))):
        times = edges[:-1][phase::2]+.25
        for k, curve in zip(ks, mean):
            ax.plot(times, curve[phase::2], marker='.', markersize=4,
                    linewidth=1.2, label=f'K={k:g}')
        ax.set(title=title, ylabel=r'Mean time QFI $4\mathrm{Var}(H)$')
        ax.grid(alpha=.3)
        ax.set_ylim(bottom=0)
    axes[0].legend(ncol=4, fontsize=9)
    axes[1].set(xlabel='Physical time t (drive period = 1)', xlim=(0, edges[-1]))
    fig.suptitle('Physical-time QFI: same 100 coherent states, j=10')
    fig.text(.5, .012, 'Points sample phase midpoints; connecting lines guide the eye, not the within-period evolution.',
             ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .035, 1, .96))
    fig.savefig(target, dpi=180)
    plt.close(fig)
    print(target)


if __name__ == '__main__':
    main()
