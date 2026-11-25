"""Propagate a preparation family and its parameter derivative under common U_K."""
import json
import numpy as np
from kicked_top_time_translation_qfi_j10 import OUT, KS, N, J, spin, expm, plt

STEM = 'kicked_top_preparation_Jy_qfi_pair_normalized'
PARAMETERS = (0., .07)


def qfi(psi, tangent):
    return 4*(np.vdot(tangent, tangent).real-abs(np.vdot(psi, tangent))**2)


def main():
    files = [OUT/f'{STEM}.{ext}' for ext in ('png', 'npz', 'json')]
    if any(p.exists() for p in files):
        raise FileExistsError('Refusing to overwrite preparation-QFI outputs')
    m, jy, coherent = spin()
    rotation = expm(-1j*np.pi*jy/2)
    initial = [expm(-1j*theta*jy)@coherent for theta in PARAMETERS]
    initial_qfi = np.array([qfi(p, -1j*jy@p) for p in initial])
    assert np.all(initial_qfi > 0)
    values = np.empty((2, len(KS), N+1))
    norm_error = 0.
    derivative_error = 0.
    h = 1e-5
    for a, theta in enumerate(PARAMETERS):
        for b, k in enumerate(KS):
            twist = np.exp(-1j*k*m*m/(2*J))
            psi = initial[a].copy()
            tangent = -1j*jy@psi
            plus = expm(-1j*(theta+h)*jy)@coherent
            minus = expm(-1j*(theta-h)*jy)@coherent
            for n in range(N+1):
                values[a, b, n] = qfi(psi, tangent)
                norm_error = max(norm_error, abs(np.vdot(psi, psi).real-1))
                derivative_error = max(derivative_error,
                                       np.linalg.norm((plus-minus)/(2*h)-tangent))
                if n < N:
                    psi = twist*(rotation@psi)
                    tangent = twist*(rotation@tangent)
                    plus = twist*(rotation@plus)
                    minus = twist*(rotation@minus)
    ratio = values/initial_qfi[:, None, None]
    error = float(np.max(abs(ratio-1)))
    assert norm_error < 1e-11 and derivative_error < 1e-7
    assert np.isfinite(values).all() and error < 1e-10
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for a, ax in enumerate(axes):
        for b, k in enumerate(KS):
            ax.plot(np.arange(N+1), ratio[a, b], lw=1.3, label=f'K={k:g}')
        ax.set(title=f'Preparation parameter θ = {PARAMETERS[a]:g} rad',
               xlabel='Complete kicks n', xlim=(0, N), ylim=(.95, 1.05))
        ax.grid(alpha=.3)
        ax.text(.5, .15, 'All eight K curves overlap at 1',
                ha='center', transform=ax.transAxes)
    axes[0].set_ylabel(r'$F_Q^{(\theta)}(n;\theta)/F_Q^{(\theta)}(0;\theta)$')
    axes[1].legend(ncol=2, fontsize=8)
    fig.suptitle('Initial Jy-rotation parameter, then common kicked-top evolution (j=10)')
    fig.tight_layout()
    fig.savefig(files[0], dpi=180)
    plt.close(fig)
    np.savez_compressed(files[1], qfi=values, ratio=ratio, theta=PARAMETERS,
                        Ks=KS, kicks=np.arange(N+1), initial_qfi=initial_qfi)
    report = {'family': 'U_K^n exp(-i theta Jy)|coherent>',
              'coherent_angles': [1.1, .7], 'j': J, 'dimension': 2*J+1,
              'theta': PARAMETERS, 'Ks': KS, 'kicks': [0, N],
              'initial_qfi': initial_qfi.tolist(), 'max_ratio_error': error,
              'max_norm_error': float(norm_error),
              'central_difference_step': h, 'max_tangent_error': float(derivative_error),
              'precision': 'complex128', 'sampling': 'deterministic, one coherent seed'}
    files[2].write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
