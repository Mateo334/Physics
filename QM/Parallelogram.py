import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from pathlib import Path

d = 10
N = 100
n = d**2-1
dt = 10**-6
# r = np.random.rand(n)/1000
r = np.full(n, 1e-3)

def random_hermitian(d):
    X = []
    for i in range(n):
        A = np.random.normal(size = (d,d)) +1j*np.random.normal(size = (d,d))
        A = (A + A.conj().T)/2 #Hermitian
        X.append(A)
    return np.array(X)
def coefficients(A, X, G):
    comm = 1j * (A @ X - X @ A)
    F = np.einsum("lab,iba->li", X, comm).real
    return np.linalg.solve(G, F)
