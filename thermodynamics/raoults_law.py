"""
Raoult's Law Module
===================
Computes ideal vapor-liquid equilibrium K-values:
    Ki = Psat_i / P

Valid assumption for:
    - Low to moderate pressures
    - Chemically similar components
    - Ideal or near-ideal liquid mixtures

Author: Deepanshu Sati
"""

import numpy as np
from thermodynamics.antoine import antoine_psat_array


def compute_ideal_K_values(T_celsius: float, P_bar: float, params_df) -> np.ndarray:
    """
    Compute ideal K-values using Raoult's Law.

    Parameters
    ----------
    T_celsius : float
        Temperature in Celsius.
    P_bar : float
        System pressure in bar.
    params_df : pd.DataFrame
        Component properties with Antoine constants.

    Returns
    -------
    np.ndarray
        K-values for each component, shape (n_components,).
    """
    psat = antoine_psat_array(T_celsius, params_df)
    K = psat / P_bar
    return K