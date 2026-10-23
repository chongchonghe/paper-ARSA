"""Linear wave analysis of the RHD equations with the (I)RSLA.

Implements the linearized system of Krumholz & Wibking, 'An Improved
Reduced Speed of Light Approximation', section 2.3: eigenvalues Omega of

    (A + i S) dW = Omega dW,   dW = (drho, dv, dT, dE, dF)

where A is the advection matrix and S = tau * Stilde the source matrix,
tau = L / khat the optical depth across a perturbation wavelength.

phi handling:
  - exact RHD: phi = 1
  - RSLA: phi = phi0 (constant)
  - IRSLA: phi = (phi0 + (tau/3)^4) / (1 + (tau/3)^4)  [eq. phi of paper]
"""

import numpy as np


def advection_matrix(M, phi, C, gamma):
    g1 = gamma - 1.0
    A = np.array([
        [M, 1, 0, 0, 0],
        [1, M, 1, 0, 0],
        [0, g1, M, 0, 0],
        [0, 0, 0, 0, phi * C],
        [0, 0, 0, phi * C / 3.0, 0],
    ], dtype=float)
    return A


def source_matrix(M, phi, C, P, gamma, tau):
    g1 = gamma - 1.0
    D = 3.0 * C**2 - M**2
    S = tau * np.array([
        [0, 0, 0, 0, 0],
        [0, -4.0 * C * P / D, -4.0 * P * M / C, -P * M / (3.0 * C), P],
        [0, 0, -4.0 * (C**2 - M**2) * P * g1 / C,
         (3.0 * C**2 + M**2) * P * g1 / (3.0 * C), -2.0 * P * M * g1],
        [0, 4.0 * phi * C * M / D, 4.0 * phi * C, -phi * C, phi * M],
        [0, 4.0 * phi * C**2 / D, 4.0 * phi * M, phi * M / 3.0, -phi * C],
    ], dtype=float)
    return S


def phi_irsla(tau, phi0):
    x = (tau / 3.0) ** 4
    return (phi0 + x) / (1.0 + x)


def eigensystem(tau, M, C, P, gamma, phi):
    """Eigenvalues Omega and right eigenvectors of (A + iS)."""
    A = advection_matrix(M, phi, C, gamma)
    S = source_matrix(M, phi, C, P, gamma, tau)
    lam, vec = np.linalg.eig(A + 1j * S)
    return lam, vec


def classify_modes(lam, vec):
    """Sort the 5 eigenvalues into (entropy, acoustic-, acoustic+, rad-, rad+)
    by tracking dominant eigenvector components: radiation modes dominated by
    (dE, dF), entropy mode has smallest |Re| among matter modes."""
    # matter weight vs radiation weight
    idx = np.arange(5)
    radw = (np.abs(vec[3]) + np.abs(vec[4])) / np.abs(vec).sum(axis=0)
    rad_idx = idx[np.argsort(-radw)][:2]
    mat_idx = np.array([i for i in idx if i not in rad_idx])
    # among matter modes: entropy = smallest |Re(Omega) - mean|; simpler:
    # acoustic pair are the two with largest |Re| spread around advection speed
    re = lam[mat_idx].real
    order = np.argsort(re)
    mat_sorted = mat_idx[order]  # ascending real part: ac-, entropy, ac+
    rad_sorted = rad_idx[np.argsort(lam[rad_idx].real)]  # rad-, rad+
    return dict(entropy=lam[mat_sorted[1]],
                acoustic=(lam[mat_sorted[0]], lam[mat_sorted[2]]),
                radiation=(lam[rad_sorted[0]], lam[rad_sorted[1]]))
