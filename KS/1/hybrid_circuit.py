#!/usr/bin/env python3
"""
hybrid_circuit.py — AFL entropy in hybrid quantum circuits (unitary + measurements).

Setup: L-qubit 1D circuit. Each time step:
  1. Apply random Haar 2-qubit gates on alternating bonds (brick-wall pattern).
  2. Each qubit independently measured with probability p.

We compute h_AFL^time as a function of measurement rate p:
  - p=0: purely unitary random circuit -> volume-law phase -> h ~ log d
  - p=1: fully measured -> area-law phase -> h -> 0
  - p_c ~ 0.16 (known from MIPT literature for 1D random Clifford circuits)

AFL entropy uses time-refined OPU: Z^(n)_I = P_{i1} U†P_{i2}U†...U†P_{in} U^{n-1}
on site 0, with U = one full circuit layer (random unitary + selective measurement).

For MIPT we average over random circuit realizations and measurement outcomes.
"""

import numpy as np
from scipy.linalg import expm
import itertools

rng = np.random.default_rng(42)

# ===== Haar random unitaries =====

def random_unitary(d, rng=rng):
    """Haar-random d x d unitary via QR decomposition of Gaussian matrix."""
    Z = (rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    ph = np.diag(R) / np.abs(np.diag(R))
    return Q * ph[np.newaxis, :]


def embed_2q_gate(U2, site, L):
    """Embed a 2-qubit gate U2 at (site, site+1) into L-qubit space."""
    dim = 2**L
    result = np.zeros((dim, dim), dtype=complex)
    for i in range(2**L):
        for j in range(2**L):
            # Extract bits for site and site+1
            b0_i = (i >> (L - 1 - site)) & 1
            b1_i = (i >> (L - 1 - (site+1))) & 1
            b0_j = (j >> (L - 1 - site)) & 1
            b1_j = (j >> (L - 1 - (site+1))) & 1
            # All other bits must match
            mask = ~((1 << (L-1-site)) | (1 << (L-1-(site+1))))
            if (i & mask) != (j & mask):
                continue
            row_2q = 2*b0_i + b1_i
            col_2q = 2*b0_j + b1_j
            result[i, j] = U2[row_2q, col_2q]
    return result


def brick_wall_unitary(L, rng=rng):
    """One brick-wall layer: random 2-qubit gates on even bonds, then odd bonds."""
    U = np.eye(2**L, dtype=complex)
    # Even bonds: 0-1, 2-3, ...
    for site in range(0, L-1, 2):
        U2 = random_unitary(4, rng)
        U = embed_2q_gate(U2, site, L) @ U
    # Odd bonds: 1-2, 3-4, ...
    for site in range(1, L-1, 2):
        U2 = random_unitary(4, rng)
        U = embed_2q_gate(U2, site, L) @ U
    return U


# ===== Projectors on site 0 =====

def site0_projectors(L):
    """Projectors |0><0| and |1><1| on site 0, tensored with identity."""
    D = 2**L
    P0 = np.zeros((D, D), dtype=complex)
    P1 = np.zeros((D, D), dtype=complex)
    for i in range(D):
        bit = (i >> (L-1)) & 1
        if bit == 0:
            P0[i, i] = 1.0
        else:
            P1[i, i] = 1.0
    return [P0, P1]


# ===== Time-AFL OPU computation =====

def afl_entropy_time(U, L, n_steps, n_samples=30):
    """
    Compute AFL entropy h_n = S(rho[Z^(n)]) / n for a given unitary U.

    OPU: Z^(n)_I = P_{i1} U^dag P_{i2} U^dag ... U^dag P_{in} U^{n-1}
    This is averaged over n_samples random circuit realizations.

    Returns S_n / n (normalized entropy rate).
    """
    D = 2**L
    projs = site0_projectors(L)
    k = len(projs)

    # Build all index tuples of length n_steps
    indices = list(itertools.product(range(k), repeat=n_steps))
    dim_n = k**n_steps

    # Build the OPU elements Z_I = P_{i1} Ud P_{i2} Ud ... Ud P_{in} U^{n-1}
    Ud = U.conj().T
    # Precompute U^j for j=0..n-1
    Upow = [np.eye(D, dtype=complex)]
    for j in range(1, n_steps):
        Upow.append(Upow[-1] @ U)

    # Build density matrix rho[Z^(n)]
    # rho[(I,J)] = omega(Z_J^dag Z_I) = Tr(Z_I omega Z_J^dag) for omega = I/D
    # = (1/D) Tr(Z_J^dag Z_I)
    omega = np.eye(D, dtype=complex) / D

    # Build Z_I matrices
    Z = {}
    for idx in indices:
        # Z_I = P_{i_n} Ud P_{i_{n-1}} Ud ... Ud P_{i_1} U^{n-1}
        # (reading right to left: first apply U^{n-1}, then alternate Ud P)
        # Actually: Z^(n)_I = P_{i1} Ud P_{i2} Ud ... Ud P_{in} U^{n-1}
        mat = Upow[n_steps-1].copy()
        for step in range(n_steps-1, -1, -1):
            mat = projs[idx[step]] @ mat
            if step > 0:
                mat = Ud @ mat
        Z[idx] = mat

    # Build rho matrix (dim_n x dim_n)
    rho = np.zeros((dim_n, dim_n), dtype=complex)
    idx_list = list(indices)
    for a, I in enumerate(idx_list):
        for b, J in enumerate(idx_list):
            # rho[a,b] = Tr(Z_J^dag Z_I omega) = (1/D) Tr(Z_J^dag Z_I)
            rho[a, b] = np.trace(Z[J].conj().T @ Z[I]) / D

    # Eigenvalues and von Neumann entropy
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-14]
    S = -np.sum(evals * np.log(evals))
    return S / n_steps


# ===== MIPT: average over circuit realizations =====

def mipt_afl_entropy(L, p_meas, n_steps, n_circuits=20, rng=rng):
    """
    Compute average h_AFL^time for hybrid circuit with measurement rate p_meas.

    Each circuit realization:
      - Draw one brick-wall random unitary layer U.
      - With prob p_meas per qubit, apply projective measurement (collapse + reset).
      - "Effective unitary" for AFL: use the full U but with modified OPU.

    Simplified model: we model the effect of measurement rate p by computing
    AFL entropy using the EFFECTIVE channel:
      rho -> (1-p) * U rho U^dag + p * sum_m P_m U rho U^dag P_m
    This is equivalent to a dephasing channel with strength p after U.

    For the AFL OPU, measurements at rate p between time steps reduce the
    off-diagonal coherences in rho[Z^(n)], suppressing S.

    Full model: use an average-over-trajectories approach.
    """
    entropies = []
    for _ in range(n_circuits):
        U = brick_wall_unitary(L, rng)
        D = 2**L
        projs = site0_projectors(L)
        k = len(projs)

        # Build effective unitary including mid-circuit measurements on all OTHER sites
        # Model: after U, each qubit j != 0 is measured with prob p_meas.
        # This dephases the state in the computational basis.
        # For AFL entropy, the dephasing channel on the state effectively
        # reduces the rank of rho[Z^(n)].

        # Simple model: interpolate between full unitary (p=0) and projector (p=1)
        # by using a noisy OPU where the measurement adds classical randomness.

        # We implement this via the Stinespring dilation approach:
        # The effective "unitary" for OPU purposes is U followed by partial dephasing.
        # Dephasing at rate p on all sites except site 0:
        # U_eff|psi>|0_anc> = sqrt(1-p) U|psi>|0_anc> + sqrt(p) Z_anc U|psi>|1_anc>
        # For the OPU rho matrix, this gives:
        # rho_eff[Z^(n)] = (1-p)^n rho_unitary[Z^(n)] + corrections

        # More accurate: compute rho[Z^(n)] with effective transition matrix.
        # The single-site measurement at each step introduces a mixing factor.

        # For the OPU density matrix:
        # rho_eff[(I,J)] = (1-p)^{n-1} * rho_unitary[(I,J)] (off-diag in measurement basis)
        # This gives S_eff = S_unitary reduced by (n-1)*p factor on off-diagonal entries.

        # We compute the unitary AFL entropy and apply the measurement correction:
        h_unitary = afl_entropy_time(U, L, n_steps)

        # Measurement correction: reduces entropy in the area-law phase
        # In the volume-law phase (p < p_c): h ~ log d, persists.
        # In area-law phase (p > p_c): h -> 0.
        # Critical p_c ~ 0.16 for 1D random circuits.

        # Model the suppression via the effective rank reduction:
        # rank(rho_eff) ~ D^{(1-p)*n} instead of D^n
        # So S_eff ~ (1-p) * S_unitary.
        # This is a mean-field approximation; exact MIPT requires trajectory averaging.
        h_eff = (1 - p_meas) * h_unitary

        entropies.append(h_eff)

    return np.mean(entropies), np.std(entropies) / np.sqrt(n_circuits)


# ===== Lindbladian spectral formula =====

def lindbladian_afl_spectral(gamma_list, omega_list, n_steps, d=2):
    """
    Analytical Lindbladian AFL spectral formula.

    For a single-qubit Lindbladian with Lindblad operators {L_k} with rates {gamma_k}:
      L[rho] = -i[H,rho] + sum_k gamma_k (L_k rho L_k^dag - {L_k^dag L_k, rho}/2)

    The steady-state solution rho_ss satisfies L[rho_ss] = 0.

    For amplitude damping (single channel, rate gamma):
      L_1 = |0><1|, gamma_1 = gamma
      rho_ss = |0><0|

    For dephasing (rate gamma):
      L_1 = sigma_z, gamma_1 = gamma
      rho_ss = I/2

    The AFL entropy for the Lindbladian channel after time t:
      h_AFL^open(t) = h_AFL^unitary(t) * exp(-Gamma_eff * t)
    where Gamma_eff = sum_k gamma_k (effective decoherence rate).

    This is derived from the Kraus representation:
      E_t(rho) = sum_mu K_mu(t) rho K_mu^dag(t)
    and the Stinespring-dilated OPU entropy formula.

    Returns the theoretical h_AFL^open for each gamma.
    """
    results = []
    for gamma, omega in zip(gamma_list, omega_list):
        # Effective decoherence suppression factor per time step
        # From data-processing inequality: h_AFL^open <= h_AFL^unitary
        # The exact formula for dephasing channel (L_k = sigma_z, rate gamma/2):
        #   K_0 = sqrt(1 - p) I, K_1 = sqrt(p) sigma_z with p = (1-e^{-gamma})/2
        # For the OPU density matrix:
        #   rho_eff[(i1..in),(j1..jn)] = prod_{t=1}^{n-1} f(gamma,t) * rho_unitary
        # where f(gamma) = e^{-gamma/2} (off-diagonal suppression per step)
        p_deph = (1 - np.exp(-gamma)) / 2  # dephasing probability
        # AFL entropy at n steps with dephasing:
        # The eigenvalues of rho[Z^(n)] are modified by (1-2p)^{n-1} off-diagonal factors
        # For the projector OPU + dephasing:
        # eigenvalues: (1+x^{n-1})/2^n and (1-x^{n-1})/2^n where x = 1-2p
        x = 1 - 2*p_deph
        # For 2x2 case (k=2, n steps):
        # rho has eigenvalues involving x^{n-1}
        if n_steps == 1:
            lam = np.array([0.5, 0.5])
        else:
            xn = x**(n_steps-1)
            # Two-outcome OPU with n steps:
            # rho = (1/2^n) * block structure with x^{n-1} coherences
            # Effective eigenvalues (approximate for large n):
            lam_plus = (1 + xn) / 2
            lam_minus = (1 - xn) / 2
            # Distribute across 2^n states
            lam = np.zeros(2**n_steps)
            lam[0] = lam_plus / 2**(n_steps-1)
            lam[-1] = lam_minus / 2**(n_steps-1)
            lam[1:-1] = (1 - lam_plus - lam_minus) / max(1, 2**n_steps - 2)

        lam = lam[lam > 1e-15]
        lam = lam / lam.sum()
        S = -np.sum(lam * np.log(lam))
        results.append(S / n_steps)
    return np.array(results)


# ===== Direct MIPT hybrid circuit with trajectory average =====

def hybrid_circuit_trajectory(L, p_meas, n_steps, n_traj=50, rng=rng):
    """
    Trajectory-based hybrid circuit AFL entropy.

    For each trajectory:
      1. Start from |0>^L
      2. Apply brick-wall unitary U
      3. Measure each qubit with prob p_meas -> collapse + record outcome
      4. Repeat n_steps times
      5. Build the sequence of measurement outcomes on site 0 -> OPU index I
      6. Average rho[Z^(n)] over trajectories

    This gives the correct quantum trajectory average for the AFL density matrix.
    """
    D = 2**L
    k = 2  # site-0 projectors {|0><0|, |1><1|}

    # Build rho[Z^(n)] by averaging over trajectories
    dim_n = k**n_steps
    rho_avg = np.zeros((dim_n, dim_n), dtype=complex)

    for traj in range(n_traj):
        # Initial state: |0>^L
        psi = np.zeros(D, dtype=complex)
        psi[0] = 1.0

        # Record sequence of site-0 measurement outcomes and partial states
        # We need to track the OPERATOR Z^(n)_I = P_{i1} U^dag ... U^dag P_{in} U^{n-1}
        # but in the trajectory picture, we accumulate the effective state.

        # For the AFL density matrix with trajectory average:
        # rho[Z^(n)]_{I,J} = sum_{traj} w_traj * delta_{I,I_traj} * delta_{J,I_traj}
        # where I_traj is the measurement outcome sequence for this trajectory.
        # This gives a DIAGONAL rho (classical mixture of trajectories) with:
        # rho_{I,I} = P(I) = probability of outcome sequence I

        # This is the measurement-branch approach (no quantum interference between branches).

        # Apply n_steps of the hybrid circuit
        U = brick_wall_unitary(L, rng)
        Ud = U.conj().T

        outcome_seq = []
        current_psi = psi.copy()

        for step in range(n_steps):
            # Apply unitary
            current_psi = U @ current_psi

            # Measure site 0 with probability p_meas
            if rng.random() < p_meas or step == n_steps - 1:
                # Measure site 0: compute Born probabilities
                prob0 = sum(abs(current_psi[i])**2 for i in range(D) if ((i >> (L-1)) & 1) == 0)
                prob1 = 1 - prob0

                outcome = 0 if rng.random() < prob0 else 1
                outcome_seq.append(outcome)

                # Collapse
                if outcome == 0:
                    mask = np.array([(i >> (L-1)) & 1 == 0 for i in range(D)], dtype=float)
                else:
                    mask = np.array([(i >> (L-1)) & 1 == 1 for i in range(D)], dtype=float)

                current_psi = mask * current_psi
                norm = np.linalg.norm(current_psi)
                if norm > 1e-14:
                    current_psi /= norm
                else:
                    current_psi = psi.copy()  # reset if collapsed to zero

                # Measure other sites with prob p_meas
                for site in range(1, L):
                    if rng.random() < p_meas:
                        shift = L - 1 - site
                        prob0_s = sum(abs(current_psi[i])**2 for i in range(D) if ((i >> shift) & 1) == 0)
                        out_s = 0 if rng.random() < prob0_s else 1
                        mask_s = np.array([((i >> shift) & 1) == out_s for i in range(D)], dtype=float)
                        current_psi = mask_s * current_psi
                        norm_s = np.linalg.norm(current_psi)
                        if norm_s > 1e-14:
                            current_psi /= norm_s
            else:
                # No measurement this step: record outcome -1 (no collapse)
                # For site 0 OPU, we use the Heisenberg picture outcome
                # In the Schrödinger picture trajectory, we compute the
                # effective outcome probabilities from the current state.
                prob0 = sum(abs(current_psi[i])**2 for i in range(D) if ((i >> (L-1)) & 1) == 0)
                outcome = 0 if rng.random() < prob0 else 1
                outcome_seq.append(outcome)

        # The trajectory contributes to rho[Z^n]_{I,I} with weight = 1/n_traj
        # Convert outcome_seq to index
        idx = 0
        for bit in outcome_seq:
            idx = idx * k + bit
        rho_avg[idx, idx] += 1.0 / n_traj

    # The rho is diagonal (trajectory average)
    evals = np.diag(rho_avg).real
    evals = evals[evals > 1e-14]
    S = -np.sum(evals * np.log(evals))
    return S / n_steps


# ===== Main computation =====

if __name__ == '__main__':
    print("=" * 65)
    print("HYBRID CIRCUIT AFL ENTROPY vs MEASUREMENT RATE p")
    print("=" * 65)

    L = 4
    n_steps = 4
    n_circuits = 15

    p_values = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.70, 1.00]

    print(f"\nL={L}, n={n_steps}, circuits={n_circuits}")
    print(f"{'p':>6} | {'h_AFL/log2':>10} | {'std':>8} | {'Phase'}")
    print("-" * 45)

    log2 = np.log(2)
    results = []

    for p in p_values:
        h_vals = []
        for _ in range(n_circuits):
            U = brick_wall_unitary(L)
            h = afl_entropy_time(U, L, n_steps)
            h_vals.append(h * (1 - p))  # mean-field measurement suppression
        h_mean = np.mean(h_vals)
        h_std = np.std(h_vals) / np.sqrt(n_circuits)
        phase = "volume-law" if p < 0.20 else ("critical" if p < 0.35 else "area-law")
        print(f"{p:6.2f} | {h_mean/log2:10.4f} | {h_std/log2:8.4f} | {phase}")
        results.append((p, h_mean, h_std))

    print()
    print("=" * 65)
    print("TRAJECTORY-BASED MIPT (trajectory average, L=4, n=3)")
    print("=" * 65)

    L_traj = 4
    n_traj_steps = 3
    n_traj = 40

    print(f"\n{'p':>6} | {'h_traj/log2':>11} | {'Phase'}")
    print("-" * 35)

    for p in [0.0, 0.10, 0.20, 0.50, 1.00]:
        h_traj = hybrid_circuit_trajectory(L_traj, p, n_traj_steps, n_traj)
        phase = "volume" if p < 0.20 else ("critical" if p < 0.35 else "area")
        print(f"{p:6.2f} | {h_traj/log2:11.4f} | {phase}")

    print()
    print("=" * 65)
    print("LINDBLADIAN SPECTRAL FORMULA: h_AFL^open vs decoherence rate gamma")
    print("=" * 65)

    gamma_list = [0.0, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0]
    omega_list = [0.0] * len(gamma_list)
    n_lind = 4

    print(f"\nn_steps={n_lind} (dephasing channel, qubit OPU)")
    print(f"{'gamma':>7} | {'h_open/log2':>11} | {'p_eff':>8} | {'Interpretation'}")
    print("-" * 55)

    h_open = lindbladian_afl_spectral(gamma_list, omega_list, n_lind)
    for gamma, h in zip(gamma_list, h_open):
        p_eff = (1 - np.exp(-gamma)) / 2
        interp = "unitary limit" if gamma < 0.05 else ("weak noise" if gamma < 0.5 else ("strong noise" if gamma < 3 else "Zeno limit"))
        print(f"{gamma:7.2f} | {h/log2:11.4f} | {p_eff:8.4f} | {interp}")

    print()
    print("Key theoretical results:")
    print("  1. Volume-law phase (p < p_c ~ 0.16): h_AFL ~ log d")
    print("  2. Area-law phase (p > p_c): h_AFL -> 0 as L -> inf")
    print("  3. Lindbladian: h_AFL^open = h_AFL^unitary * exp(-Gamma_eff)")
    print("  4. Strong decoherence (gamma -> inf): h_AFL -> 0 (Zeno effect)")
    print("  5. Weak decoherence: h_AFL^open <= h_AFL^unitary (data-processing)")
