"""
Rachford-Rice Equation Solver
==============================
Solves the Rachford-Rice equation for vapor fraction ψ (psi):

    f(ψ) = Σ [ zi*(Ki - 1) / (1 + ψ*(Ki - 1)) ] = 0

Analytical bounds for the root bracket (Whitson & Brulé, 2000):
    ψ_min = 1 / (1 - K_max)
    ψ_max = 1 / (1 - K_min)

The function is guaranteed monotonically decreasing → Brent's method converges.

Author: Deepanshu Sati
"""

import numpy as np
from scipy.optimize import brentq


def rachford_rice(psi: float, z: np.ndarray, K: np.ndarray) -> float:
    """
    Evaluate the Rachford-Rice residual function.

    Parameters
    ----------
    psi : float
        Vapor fraction (0 ≤ ψ ≤ 1).
    z : np.ndarray
        Feed mole fractions, shape (n,).
    K : np.ndarray
        Equilibrium K-values, shape (n,).

    Returns
    -------
    float
        Residual f(ψ). Zero when ψ is the correct vapor fraction.
    """
    numerator = z * (K - 1.0)
    denominator = 1.0 + psi * (K - 1.0)

    # Guard against denominator singularities
    denominator = np.where(np.abs(denominator) < 1e-12,
                           np.sign(denominator) * 1e-12,
                           denominator)
    return np.sum(numerator / denominator)


def rachford_rice_derivative(psi: float, z: np.ndarray, K: np.ndarray) -> float:
    """
    Analytical derivative df/dψ (useful for Newton-Raphson convergence check).

    df/dψ = -Σ [ zi*(Ki-1)² / (1 + ψ*(Ki-1))² ]
    """
    num = z * (K - 1.0) ** 2
    denom = (1.0 + psi * (K - 1.0)) ** 2
    return -np.sum(num / denom)


def compute_psi_bounds(K: np.ndarray) -> tuple[float, float]:
    """
    Compute valid bounds for ψ to avoid singularities.

    Parameters
    ----------
    K : np.ndarray
        K-values for all components.

    Returns
    -------
    psi_min, psi_max : float
        Valid bracket for Brent's method.
    """
    K_min = np.min(K)
    K_max = np.max(K)

    # Theoretical bounds (strict)
    psi_lower = 1.0 / (1.0 - K_max)
    psi_upper = 1.0 / (1.0 - K_min)

    # Add small tolerance to stay inside the bracket
    eps = 1e-8
    psi_min = psi_lower + eps
    psi_max = psi_upper - eps

    return psi_min, psi_max


def solve_rachford_rice(z: np.ndarray, K: np.ndarray,
                         tol: float = 1e-10) -> tuple[float, str]:
    """
    Solve Rachford-Rice equation for vapor fraction ψ.

    First checks if mixture is single-phase (all liquid or all vapor),
    then applies Brent's method within the valid bracket.

    Parameters
    ----------
    z : np.ndarray
        Feed mole fractions (must sum to 1).
    K : np.ndarray
        Equilibrium K-values.
    tol : float
        Convergence tolerance for Brent's method.

    Returns
    -------
    psi : float
        Vapor fraction (0 ≤ ψ ≤ 1).
    phase : str
        'two-phase', 'all-vapor', or 'all-liquid'.
    """
    # --- Single-phase checks ---
    # If f(0) < 0 → all liquid (ψ = 0, bubble point not reached)
    if rachford_rice(0.0, z, K) < 0:
        return 0.0, "all-liquid"

    # If f(1) > 0 → all vapor (ψ = 1, dew point exceeded)
    if rachford_rice(1.0, z, K) > 0:
        return 1.0, "all-vapor"

    # --- Two-phase: solve with Brent's method ---
    psi_min, psi_max = compute_psi_bounds(K)

    # Clamp to [0, 1] to handle edge cases
    psi_min = max(psi_min, 1e-8)
    psi_max = min(psi_max, 1.0 - 1e-8)

    psi_solution = brentq(
        rachford_rice,
        psi_min,
        psi_max,
        args=(z, K),
        xtol=tol,
        rtol=tol,
        maxiter=500
    )

    return psi_solution, "two-phase"