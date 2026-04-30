"""
Flash Calculation Orchestrator
================================
Implements the complete isothermal flash algorithm (T, P flash).

Algorithm flow:
    1. Load component properties
    2. Compute initial K-values (Wilson correlation as initial guess)
    3. Solve Rachford-Rice → vapor fraction ψ
    4. Compute phase compositions x_i and y_i
    5. (Real model) Update K-values from PR-EOS fugacity coefficients
    6. Repeat until convergence
    7. Return FlashResult

Convergence criterion:
    max|Ki_new/Ki_old - 1| < tolerance

Author: Deepanshu Sati
"""

import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Optional

from thermodynamics.raoults_law import compute_ideal_K_values
from thermodynamics.peng_robinson import pr_pure_params
from thermodynamics.fugacity import compute_real_K_values
from flash.rachford_rice import solve_rachford_rice


@dataclass
class FlashResult:
    """Container for flash calculation outputs."""
    psi: float                    # Vapor fraction [-]
    x: np.ndarray                 # Liquid mole fractions
    y: np.ndarray                 # Vapor mole fractions
    K: np.ndarray                 # Final K-values
    phase: str                    # 'two-phase', 'all-liquid', 'all-vapor'
    T_celsius: float              # Operating temperature [°C]
    P_bar: float                  # Operating pressure [bar]
    components: list              # Component names
    iterations: int               # PR-EOS iterations (real model)
    converged: bool               # Convergence flag


def wilson_K_estimate(Tc: np.ndarray, Pc: np.ndarray,
                       omega: np.ndarray,
                       T_kelvin: float, P_bar: float) -> np.ndarray:
    """
    Wilson correlation for initial K-value estimate.
    Ki = (Pci/P) * exp(5.373*(1+ωi)*(1 - Tci/T))

    Essential to start PR-EOS iterations in the correct region.
    """
    return (Pc / P_bar) * np.exp(5.373 * (1.0 + omega) * (1.0 - Tc / T_kelvin))


def compute_phase_compositions(psi: float, z: np.ndarray,
                                K: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute liquid and vapor mole fractions from ψ and K-values.

    x_i = z_i / (1 + ψ*(Ki - 1))
    y_i = Ki * x_i

    Compositions are normalized to handle numerical drift.
    """
    x = z / (1.0 + psi * (K - 1.0))
    y = K * x

    # Normalize (should already sum to 1, but guard against drift)
    x = x / x.sum()
    y = y / y.sum()

    return x, y


def isothermal_flash_ideal(z: np.ndarray,
                            T_celsius: float,
                            P_bar: float,
                            params_df: pd.DataFrame) -> FlashResult:
    """
    Isothermal flash using ideal thermodynamics (Raoult's Law).

    Parameters
    ----------
    z : np.ndarray
        Feed mole fractions.
    T_celsius : float
        Temperature [°C].
    P_bar : float
        Pressure [bar].
    params_df : pd.DataFrame
        Component database with Antoine constants.

    Returns
    -------
    FlashResult
    """
    z = np.asarray(z, dtype=float)
    z = z / z.sum()  # Normalize feed

    # --- Step 1: Compute ideal K-values ---
    K = compute_ideal_K_values(T_celsius, P_bar, params_df)

    # --- Step 2: Solve Rachford-Rice ---
    psi, phase = solve_rachford_rice(z, K)

    # --- Step 3: Phase compositions ---
    if phase == "all-liquid":
        x, y = z.copy(), z.copy()
    elif phase == "all-vapor":
        x, y = z.copy(), z.copy()
    else:
        x, y = compute_phase_compositions(psi, z, K)

    return FlashResult(
        psi=psi, x=x, y=y, K=K, phase=phase,
        T_celsius=T_celsius, P_bar=P_bar,
        components=list(params_df.index),
        iterations=1, converged=True
    )


def isothermal_flash_pr(z: np.ndarray,
                         T_celsius: float,
                         P_bar: float,
                         params_df: pd.DataFrame,
                         kij: Optional[np.ndarray] = None,
                         max_iter: int = 100,
                         tol: float = 1e-8) -> FlashResult:
    """
    Isothermal flash using Peng-Robinson EOS (real thermodynamics).

    Successive substitution algorithm:
        1. Wilson K-values as initial estimate
        2. Solve Rachford-Rice → ψ, x, y
        3. Compute fugacity coefficients for liquid and vapor
        4. Update K = phi_L / phi_V
        5. Check convergence: max|ΔK/K| < tol
        6. Repeat until converged

    Parameters
    ----------
    z : np.ndarray
        Feed mole fractions.
    T_celsius : float
        Temperature [°C].
    P_bar : float
        Pressure [bar].
    params_df : pd.DataFrame
        Component database with Tc, Pc, omega.
    kij : np.ndarray, optional
        Binary interaction parameters (n×n matrix).
    max_iter : int
        Maximum successive substitution iterations.
    tol : float
        Convergence tolerance.

    Returns
    -------
    FlashResult
    """
    z = np.asarray(z, dtype=float)
    z = z / z.sum()

    T_kelvin = T_celsius + 273.15

    Tc = params_df["Tc"].values
    Pc = params_df["Pc"].values
    omega = params_df["omega"].values

    # --- Step 1: Wilson K-value initial estimate ---
    K = wilson_K_estimate(Tc, Pc, omega, T_kelvin, P_bar)

    converged = False
    iterations = 0
    x, y = z.copy(), z.copy()

    for iteration in range(max_iter):
        iterations = iteration + 1

        # --- Step 2: Solve Rachford-Rice ---
        psi, phase = solve_rachford_rice(z, K)

        if phase == "all-liquid":
            x, y = z.copy(), z.copy()
            converged = True
            break
        elif phase == "all-vapor":
            x, y = z.copy(), z.copy()
            converged = True
            break

        x, y = compute_phase_compositions(psi, z, K)

        # --- Step 3: Update K from PR-EOS fugacity coefficients ---
        K_new = compute_real_K_values(z, x, y, Tc, Pc, omega,
                                       T_kelvin, P_bar, kij)

        # --- Step 4: Check convergence ---
        K_ratio = K_new / np.where(K > 1e-15, K, 1e-15)
        error = np.max(np.abs(K_ratio - 1.0))

        K = K_new

        if error < tol:
            converged = True
            break

    return FlashResult(
        psi=psi, x=x, y=y, K=K, phase=phase,
        T_celsius=T_celsius, P_bar=P_bar,
        components=list(params_df.index),
        iterations=iterations, converged=converged
    )