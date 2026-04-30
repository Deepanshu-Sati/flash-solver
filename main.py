"""
Multicomponent Flash Solver — CLI Entry Point
=============================================
Run a flash calculation from the command line.

Usage:
    python main.py --components methane ethane propane
                   --feed 0.5 0.3 0.2
                   --temperature 25
                   --pressure 10
                   --model pr

Author: Deepanshu Sati
"""

import argparse
import numpy as np
import pandas as pd
import sys

from flash.flash_calculation import isothermal_flash_ideal, isothermal_flash_pr
from visualization.phase_diagrams import plot_composition_profile
import matplotlib.pyplot as plt


def load_database(path: str = "data/component_database.csv") -> pd.DataFrame:
    """Load component property database."""
    df = pd.read_csv(path, index_col="component")
    return df


def print_results(result) -> None:
    """Pretty-print flash results to terminal."""
    print("\n" + "=" * 60)
    print("  FLASH CALCULATION RESULTS")
    print("=" * 60)
    print(f"  Temperature : {result.T_celsius:.2f} °C")
    print(f"  Pressure    : {result.P_bar:.4f} bar")
    print(f"  Phase state : {result.phase.upper()}")
    print(f"  Vapor frac. : ψ = {result.psi:.6f}")
    if hasattr(result, 'iterations'):
        print(f"  Iterations  : {result.iterations}")
        print(f"  Converged   : {result.converged}")
    print("-" * 60)
    print(f"  {'Component':<15} {'z_i':>8}  {'x_i':>8}  {'y_i':>8}  {'K_i':>8}")
    print("-" * 60)
    z = result.psi * result.y + (1 - result.psi) * result.x
    for i, comp in enumerate(result.components):
        print(f"  {comp:<15} {z[i]:>8.4f}  {result.x[i]:>8.4f}  "
              f"{result.y[i]:>8.4f}  {result.K[i]:>8.4f}")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Multicomponent Flash Solver — Ideal & PR-EOS"
    )
    parser.add_argument("--components", nargs="+", required=True,
                        help="Component names from database")
    parser.add_argument("--feed", nargs="+", type=float, required=True,
                        help="Feed mole fractions (need not sum to 1, auto-normalized)")
    parser.add_argument("--temperature", type=float, required=True,
                        help="Temperature [°C]")
    parser.add_argument("--pressure", type=float, required=True,
                        help="Pressure [bar]")
    parser.add_argument("--model", choices=["ideal", "pr"], default="ideal",
                        help="Thermodynamic model: 'ideal' or 'pr'")
    parser.add_argument("--plot", action="store_true",
                        help="Show composition bar chart")

    args = parser.parse_args()

    # --- Load database ---
    try:
        db = load_database()
    except FileNotFoundError:
        print("ERROR: data/component_database.csv not found.")
        sys.exit(1)

    # --- Validate components ---
    for comp in args.components:
        if comp not in db.index:
            print(f"ERROR: Component '{comp}' not in database.")
            print(f"Available: {list(db.index)}")
            sys.exit(1)

    params_df = db.loc[args.components]
    z = np.array(args.feed, dtype=float)

    if len(z) != len(args.components):
        print("ERROR: Feed composition length must match number of components.")
        sys.exit(1)

    # --- Run flash ---
    if args.model == "ideal":
        result = isothermal_flash_ideal(z, args.temperature, args.pressure, params_df)
    else:
        result = isothermal_flash_pr(z, args.temperature, args.pressure, params_df)

    # --- Print results ---
    print_results(result)

    # --- Optional plot ---
    if args.plot:
        fig = plot_composition_profile(result)
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()