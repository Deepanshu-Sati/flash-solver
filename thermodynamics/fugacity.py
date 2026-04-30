"""
Fugacity Coefficient Module (PR-EOS)
=====================================
Computes fugacity coefficients for each component in a mixture
using the Peng-Robinson EOS analytical expression.

ln(phi_i) = (bi/bm)*(Z-1) - ln(Z-B)
            - A/(2√2 * B) * (2*sum_j[yj*aij]/am - bi/bm)
            * ln[(Z + (1+√2)*B) / (Z + (1-√2)*B)]

Real K-values:
    Ki = phi_i_liquid / phi_i_vapor

Author: Deepanshu Sati
"""

import numpy as np
from thermodynamics.peng_robinson import (
    pr_pure_params, pr_mixing_rules, compute_A_B, solve_Z_factor, R
)


def fugacity_coefficients(composition: np.ndarray,
                           a_pure: np.ndarray,
                           b_pure: np.ndarray,
                           T_kelvin: float,
                           P_bar: float,
                           phase: str,
                           kij: np.ndarray = None) -> np.ndarray:
    """
    Compute fugacity coefficients for all components in a phase.

    Parameters
    ----------
    composition : np.ndarray
        Mole fractions [x_i or y_i], shape (n,).
    a_pure : np.ndarray
        Pure component PR-EOS a parameters.
    b_pure : np.ndarray
        Pure component PR-EOS b parameters.
    T_kelvin : float
        Temperature [K].
    P_bar : float
        Pressure [bar].
    phase : str
        'vapor' or 'liquid' — determines Z root selection.
    kij : np.ndarray, optional
        Binary interaction parameters.

    Returns
    -------
    np.ndarray
        Fugacity coefficients phi_i, shape (n,).
    """
    n = len(composition)
    if kij is None:
        kij = np.zeros((n, n))

    # --- Mixing rules ---
    a_mix, b_mix = pr_mixing_rules(a_pure, b_pure, composition, kij)

    # --- Dimensionless parameters ---
    A, B = compute_A_B(a_mix, b_mix, T_kelvin, P_bar)

    # --- Z-factor ---
    Z = solve_Z_factor(A, B, phase=phase)

    # --- Cross parameters a_ij ---
    a_ij = np.sqrt(np.outer(a_pure, a_pure)) * (1.0 - kij)

    # --- sum_j(yj * aij) for each i ---
    sum_aij = a_ij @ composition   # shape (n,)

    # --- Fugacity coefficient formula ---
    sqrt2 = np.sqrt(2.0)
    ln_phi = (
        (b_pure / b_mix) * (Z - 1.0)
        - np.log(Z - B)
        - (A / (2.0 * sqrt2 * B))
        * (2.0 * sum_aij / a_mix - b_pure / b_mix)
        * np.log((Z + (1.0 + sqrt2) * B) / (Z + (1.0 - sqrt2) * B))
    )

    return np.exp(ln_phi)


def compute_real_K_values(z: np.ndarray,
                           x: np.ndarray,
                           y: np.ndarray,
                           Tc: np.ndarray,
                           Pc: np.ndarray,
                           omega: np.ndarray,
                           T_kelvin: float,
                           P_bar: float,
                           kij: np.ndarray = None) -> np.ndarray:
    """
    Compute K-values from PR-EOS fugacity coefficients.
        Ki = phi_i_liquid / phi_i_vapor

    Parameters
    ----------
    z : np.ndarray
        Feed composition (used for initialization only).
    x : np.ndarray
        Liquid phase mole fractions.
    y : np.ndarray
        Vapor phase mole fractions.
    Tc, Pc, omega : np.ndarray
        Critical properties and acentric factor.
    T_kelvin : float
        Temperature [K].
    P_bar : float
        Pressure [bar].
    kij : np.ndarray, optional
        Binary interaction parameters.

    Returns
    -------
    np.ndarray
        K-values for each component.
    """
    a_pure, b_pure = pr_pure_params(Tc, Pc, omega, T_kelvin)

    phi_liq = fugacity_coefficients(x, a_pure, b_pure, T_kelvin, P_bar,
                                     phase="liquid", kij=kij)
    phi_vap = fugacity_coefficients(y, a_pure, b_pure, T_kelvin, P_bar,
                                     phase="vapor", kij=kij)

    # Avoid division by zero
    K = phi_liq / np.where(phi_vap > 1e-15, phi_vap, 1e-15)
    return K