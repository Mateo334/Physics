# Implementing the Benettini algorithm to compute the Lyapunov exponents
import numpy as np
from scipy.integrate import solve_ivp


def H_Henon(z):
    x, y, px, py = z
    return (0.5 * (px**2 + py**2)+ 0.5 * (x**2 + y**2)+ x**2 * y- y**3 / 3)
def grad_H_Henon(z):
    x, y, px, py = z
    return np.array([x + 2*x*y,y + x**2 - y**2,px,py])
def hess_H_Henon(z):
    x, y, px, py = z
    return np.array([
        [1 + 2*y,  2*x,      0, 0],
        [2*x,      1 - 2*y,  0, 0],
        [0,        0,        1, 0],
        [0,        0,        0, 1]
    ])


omega = 1.0
def H_HO(z):
    q, p = z
    return (0.5 * p**2+ 0.5 * omega**2 * q**2)
def grad_H_HO(z):
    q, p = z
    return np.array([omega**2 * q,p])
def hess_H_HO(z):
    return np.array([
        [omega**2, 0],
        [0,        1]
    ])


def H_inverted(z):
    q, p = z
    omega = 1.0
    return 0.5*p**2 - 0.5*omega**2*q**2

def grad_H_inverted(z):
    q, p = z
    omega = 1.0
    return np.array([
        -omega**2 * q,
        p
    ])

def hess_H_inverted(z):
    omega = 1.0
    return np.array([
        [-omega**2, 0.0],
        [0.0,        1.0]
    ])


dt = 0.5
N_steps = 2000


def combined_equations(t, state):
    """The derived equations."""
    z = state[:dim]
    Z = state[dim:].reshape(dim, dim)

    dz = J @ grad_H(z)
    dZ = J @ hess_H(z) @ Z

    return np.concatenate([dz, dZ.ravel()])


def evolve_one_interval(z, Z, dt):
    """Numerically integrates the physical and tangent equations over one interval dt."""
    initial_state = np.concatenate([z, Z.ravel()])

    solution = solve_ivp(
        combined_equations,
        [0, dt],
        initial_state,
        rtol=1e-10,
        atol=1e-12
    )

    final_state = solution.y[:, -1]
    z_new = final_state[:dim]
    Z_new = final_state[dim:].reshape(dim, dim)
    return z_new, Z_new


systems = [
    ("Henon-Heiles", H_Henon, grad_H_Henon, hess_H_Henon,
     np.array([0.1, 0.1, 0.4, 0.1])),

    ("Harmonic oscillator", H_HO, grad_H_HO, hess_H_HO,
     np.array([1.0, 0.1])),

    ("Inverted harmonic oscillator", H_inverted,
     grad_H_inverted, hess_H_inverted,
     np.array([1.0, 0.1]))
]

for name, H, grad_H, hess_H, z in systems:
    dim = len(z)
    n = dim // 2

    J = np.block([
        [np.zeros((n, n)),  np.eye(n)],
        [-np.eye(n),        np.zeros((n, n))]
    ])

    Z = np.eye(dim)
    lyaps = np.zeros(dim)
    
    #QR decomposition
    for k in range(N_steps):
        z, Z = evolve_one_interval(z, Z, dt)

        Q, R = np.linalg.qr(Z)
        lyaps += np.log(np.abs(np.diag(R)))
        Z = Q

    lyaps = lyaps / (N_steps * dt)

    print(name)
    print("Lyapunov exponents:", np.sort(lyaps))
    print("Sum of Lyapunov exponents:", np.sum(lyaps))

