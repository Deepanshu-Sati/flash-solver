"""
Peng-Robinson Equation of State
================================
Implements the PR-EOS for vapor-liquid equilibrium.

Key equations:
    P = RT/(V-b) - a(T)/[V(V+b) + b(V-b)]

    a_i  = 0.45724 * R^2 * Tc_i^2 / Pc_i * alpha_i(T)
    b_i  = 0.07780 * R * Tc_i / Pc_i
    kappa_i = 0.37464 + 1.54226*omega_i - 0.26992*omega_i^2
    alpha_i = [1 + kappa_i*(1 - sqrt(Tr_i))]^2

Mixing rules (van der Waals, kij=0):
    a_mix = sum_i sum_j yi*yj*sqrt(ai*aj)
    b_mix = sum_i yi*bi

Reference: Peng & Robinson, Ind. Eng. Chem. Fundam., 1976.

Author: Deepanshu Sati
"""

import numpy as np


# Universal gas constant [bar·L/(mol·K)]
R = 0.08314462


def pr_pure_params(Tc: np.ndarray, Pc: np.ndarray, omega: np.ndarray,
                   T_kelvin: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute PR-EOS pure component a_i and b_i parameters.

    Parameters
    ----------
    Tc : np.ndarray
        Critical temperatures [K], shape (n,).
    Pc : np.ndarray
        Critical pressures [bar], shape (n,).
    omega : np.ndarray
        Acentric factors [-], shape (n,).
    T_kelvin : float
        System temperature [K].

    Returns
    -------
    a : np.ndarray
        Attractive parameter [bar·L²/mol²], shape (n,).
    b : np.ndarray
        Repulsive parameter [L/mol], shape (n,).
    """
    Tr = T_kelvin / Tc                                          # Reduced temperature

    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega ** 2   # Temperature dependency factor
    alpha = (1.0 + kappa * (1.0 - np.sqrt(Tr))) ** 2           # Temperature-dependent alpha

    a = 0.45724 * (R ** 2) * (Tc ** 2) / Pc * alpha            # [bar·L²/mol²]
    b = 0.07780 * R * Tc / Pc                                   # [L/mol]

    return a, b


def pr_mixing_rules(a_pure: np.ndarray, b_pure: np.ndarray,
                    composition: np.ndarray,
                    kij: np.ndarray = None) -> tuple[float, float]:
    """
    Apply van der Waals mixing rules.

    Parameters
    ----------
    a_pure : np.ndarray
        Pure component a values [bar·L²/mol²].
    b_pure : np.ndarray
        Pure component b values [L/mol].
    composition : np.ndarray
        Mole fractions (x or y), must sum to 1.
    kij : np.ndarray, optional
        Binary interaction parameters, shape (n,n). Default zeros.

    Returns
    -------
    a_mix : float
        Mixture attractive parameter.
    b_mix : float
        Mixture repulsive parameter.
    """
    n = len(composition)
    if kij is None:
        kij = np.zeros((n, n))

    # Combine pure parameters: a_ij = sqrt(ai * aj) * (1 - kij)
    a_ij = np.sqrt(np.outer(a_pure, a_pure)) * (1.0 - kij)

    # a_mix = sum_i sum_j yi*yj*a_ij
    a_mix = float(composition @ a_ij @ composition)

    # b_mix = sum_i yi*bi
    b_mix = float(composition @ b_pure)

    return a_mix, b_mix


def solve_Z_factor(A: float, B: float, phase: str = "vapor") -> float:
    """
    Solve cubic PR-EOS for compressibility factor Z.

    Cubic form:
        Z³ - (1-B)Z² + (A-3B²-2B)Z - (AB-B²-B³) = 0

    Parameters
    ----------
    A : float
        Dimensionless parameter A = a_mix*P/(R*T)^2
    B : float
        Dimensionless parameter B = b_mix*P/(R*T)
    phase : str
        'vapor'  → returns largest real root
        'liquid' → returns smallest positive real root

    Returns
    -------
    float
        Compressibility factor Z.
    """
    coeffs = [
        1.0,
        -(1.0 - B),
        (A - 3.0 * B ** 2 - 2.0 * B),
        -(A * B - B ** 2 - B ** 3)
    ]

    roots = np.roots(coeffs)

    # Keep only real positive roots (physically meaningful)
    real_roots = roots[
        (np.abs(roots.imag) < 1e-8) & (roots.real > B)
    ].real

    if len(real_roots) == 0:
        raise ValueError(f"No valid Z-factor found. A={A:.4f}, B={B:.4f}")

    if phase == "vapor":
        return float(np.max(real_roots))   # Gas: largest root
    else:
        return float(np.min(real_roots))   # Liquid: smallest root


def compute_A_B(a_mix: float, b_mix: float,
                T_kelvin: float, P_bar: float) -> tuple[float, float]:
    """
    Compute dimensionless PR-EOS parameters A and B.

    A = a_mix * P / (R*T)^2
    B = b_mix * P / (R*T)
    """
    RT = R * T_kelvin
    A = a_mix * P_bar / (RT ** 2)
    B = b_mix * P_bar / RT
    return A, B