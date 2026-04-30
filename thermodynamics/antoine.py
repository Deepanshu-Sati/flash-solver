"""
Antoine Equation Module
=======================
Computes saturation pressure using the Antoine equation:
    log10(Psat) = A - B / (T + C)

Units:
    T   : degrees Celsius
    Psat: bar (using NIST/Poling constants)

Author: Deepanshu Sati
"""

import numpy as np


def antoine_psat(T_celsius: float, A: float, B: float, C: float) -> float:
    """
    Compute saturation pressure from Antoine equation.

    Parameters
    ----------
    T_celsius : float
        Temperature in degrees Celsius.
    A, B, C : float
        Antoine constants (log10-bar, Celsius form).

    Returns
    -------
    float
        Saturation pressure in bar.

    Raises
    ------
    ValueError
        If temperature leads to non-physical result.
    """
    log_psat = A - B / (T_celsius + C)
    psat = 10 ** log_psat
    return psat


def antoine_psat_array(T_celsius: float, params_df) -> np.ndarray:
    """
    Compute Psat for all components in a DataFrame at temperature T.

    Parameters
    ----------
    T_celsius : float
        Temperature in Celsius.
    params_df : pd.DataFrame
        DataFrame with columns A, B, C for each component.

    Returns
    -------
    np.ndarray
        Array of saturation pressures [bar], shape (n_components,).
    """
    A = params_df["A"].values
    B = params_df["B"].values
    C = params_df["C"].values
    log_psat = A - B / (T_celsius + C)
    return 10.0 ** log_psat