"""
Phase Diagram Visualization Module
====================================
Generates T-x-y, P-x-y, and vapor fraction vs pressure diagrams.

Author: Deepanshu Sati
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from flash.flash_calculation import isothermal_flash_ideal, isothermal_flash_pr


def plot_Pxy_binary(component_1: str, component_2: str,
                     params_df: pd.DataFrame,
                     T_celsius: float,
                     P_range: tuple = (0.1, 10.0),
                     n_points: int = 50,
                     model: str = "ideal") -> plt.Figure:
    """
    Generate P-x-y diagram for a binary mixture at fixed temperature.

    Parameters
    ----------
    component_1, component_2 : str
        Component names (must exist in params_df).
    params_df : pd.DataFrame
        Component database.
    T_celsius : float
        Temperature [°C].
    P_range : tuple
        (P_min, P_max) in bar.
    n_points : int
        Number of pressure points.
    model : str
        'ideal' or 'pr'.
    """
    comps = params_df.loc[[component_1, component_2]]
    pressures = np.linspace(P_range[0], P_range[1], n_points)

    x1_list, y1_list = [], []

    # Sweep over different feed compositions at each pressure
    z1_values = np.linspace(0.01, 0.99, 20)

    bubble_P, dew_P, bubble_x, dew_y = [], [], [], []

    # For each z1, find bubble and dew points (simplified approach)
    z1_range = np.linspace(0.01, 0.99, 40)

    for z1 in z1_range:
        z = np.array([z1, 1.0 - z1])
        for P in pressures:
            if model == "ideal":
                result = isothermal_flash_ideal(z, T_celsius, P, comps)
            else:
                result = isothermal_flash_pr(z, T_celsius, P, comps)

            if result.phase == "two-phase":
                bubble_x.append(result.x[0])
                bubble_P.append(P)
                dew_y.append(result.y[0])
                dew_P.append(P)
                break  # Found two-phase, take bubble/dew at this z

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(bubble_x, bubble_P, c="navy", s=10, label="Bubble curve (liquid)", zorder=3)
    ax.scatter(dew_y, dew_P, c="crimson", s=10, label="Dew curve (vapor)", zorder=3)

    ax.set_xlabel(f"Mole fraction {component_1}", fontsize=13)
    ax.set_ylabel("Pressure [bar]", fontsize=13)
    ax.set_title(f"P-x-y Diagram: {component_1}/{component_2} at T = {T_celsius:.1f}°C\n"
                 f"Model: {'Ideal (Raoult)' if model=='ideal' else 'Peng-Robinson'}",
                 fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)

    return fig


def plot_vapor_fraction_vs_pressure(z: np.ndarray,
                                     params_df: pd.DataFrame,
                                     T_celsius: float,
                                     P_range: tuple = (0.5, 20.0),
                                     n_points: int = 100,
                                     model: str = "ideal") -> plt.Figure:
    """
    Plot vapor fraction ψ vs pressure for a fixed feed at fixed T.

    Parameters
    ----------
    z : np.ndarray
        Feed mole fractions.
    params_df : pd.DataFrame
        Component database.
    T_celsius : float
        Temperature [°C].
    P_range : tuple
        (P_min, P_max) in bar.
    n_points : int
        Number of pressure points.
    model : str
        'ideal' or 'pr'.
    """
    pressures = np.linspace(P_range[0], P_range[1], n_points)
    psi_values = []
    phase_colors = []

    for P in pressures:
        if model == "ideal":
            result = isothermal_flash_ideal(z, T_celsius, P, params_df)
        else:
            result = isothermal_flash_pr(z, T_celsius, P, params_df)

        psi_values.append(result.psi)
        phase_colors.append(result.phase)

    psi_values = np.array(psi_values)

    fig, ax = plt.subplots(figsize=(9, 5))

    # Color-code by phase region
    two_phase_mask = np.array(phase_colors) == "two-phase"
    ax.plot(pressures[two_phase_mask], psi_values[two_phase_mask],
            "b-", lw=2.5, label="Two-phase region")
    ax.axhline(y=0, color="navy", ls="--", lw=1.5, alpha=0.6, label="All liquid (ψ=0)")
    ax.axhline(y=1, color="crimson", ls="--", lw=1.5, alpha=0.6, label="All vapor (ψ=1)")

    # Mark transition pressures
    if any(two_phase_mask):
        P_bubble = pressures[two_phase_mask][-1]  # Highest P in two-phase = bubble point
        P_dew = pressures[two_phase_mask][0]       # Lowest P in two-phase = dew point
        ax.axvline(P_bubble, color="navy", ls=":", alpha=0.8,
                   label=f"Bubble P ≈ {P_bubble:.2f} bar")
        ax.axvline(P_dew, color="crimson", ls=":", alpha=0.8,
                   label=f"Dew P ≈ {P_dew:.2f} bar")

    comp_names = list(params_df.index)
    z_str = ", ".join([f"{comp_names[i]}: {z[i]:.2f}" for i in range(len(z))])

    ax.set_xlabel("Pressure [bar]", fontsize=13)
    ax.set_ylabel("Vapor Fraction ψ [-]", fontsize=13)
    ax.set_title(f"Vapor Fraction vs Pressure at T = {T_celsius:.1f}°C\n"
                 f"Feed: {z_str}", fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.05, 1.05)

    return fig


def plot_composition_profile(result: object) -> plt.Figure:
    """
    Bar chart comparing feed, liquid, and vapor compositions.

    Parameters
    ----------
    result : FlashResult
        Output from isothermal_flash_ideal or isothermal_flash_pr.
    """
    n = len(result.components)
    x_pos = np.arange(n)
    width = 0.25

    fig, ax = plt.subplots(figsize=(max(8, n * 1.5), 5))

    # Feed composition (not stored in FlashResult, reconstruct as z = ψy + (1-ψ)x)
    z_reconstructed = result.psi * result.y + (1 - result.psi) * result.x

    ax.bar(x_pos - width, z_reconstructed, width, label="Feed z_i",
           color="slategray", alpha=0.85)
    ax.bar(x_pos, result.x, width, label="Liquid x_i",
           color="steelblue", alpha=0.85)
    ax.bar(x_pos + width, result.y, width, label="Vapor y_i",
           color="tomato", alpha=0.85)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(result.components, fontsize=11)
    ax.set_ylabel("Mole Fraction [-]", fontsize=12)
    ax.set_title(
        f"Flash Composition Results | T={result.T_celsius:.1f}°C, P={result.P_bar:.2f} bar\n"
        f"Phase: {result.phase} | ψ = {result.psi:.4f}",
        fontsize=12
    )
    ax.legend(fontsize=11)
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_ylim(0, 1.05)

    return fig