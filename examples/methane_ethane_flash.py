"""
Example: Methane/Ethane Binary Flash
======================================
Demonstrates a simple two-component flash at various pressures.

Expected behavior:
    - Methane (lighter) preferentially partitions to vapor phase
    - K_methane > 1, K_ethane < 1 at these conditions
"""

import numpy as np
import pandas as pd
import sys
sys.path.insert(0, "..")

import matplotlib.pyplot as plt
from flash.flash_calculation import isothermal_flash_ideal, isothermal_flash_pr
from visualization.phase_diagrams import plot_vapor_fraction_vs_pressure

# Load database
db = pd.read_csv("data/component_database.csv", index_col="component")
comps = db.loc[["methane", "ethane"]]

# Feed: 60% methane, 40% ethane
z = np.array([0.60, 0.40])
T_celsius = -20.0   # °C — below methane normal boiling point

print("=" * 55)
print("  Methane / Ethane Flash | T = -20°C")
print("=" * 55)

for P_bar in [5.0, 10.0, 20.0, 40.0]:
    result_ideal = isothermal_flash_ideal(z, T_celsius, P_bar, comps)
    result_pr = isothermal_flash_pr(z, T_celsius, P_bar, comps)

    print(f"\n  P = {P_bar:5.1f} bar")
    print(f"  {'Model':<10} {'ψ':>8}  {'x_CH4':>8}  {'y_CH4':>8}  {'K_CH4':>8}")
    print(f"  {'Ideal':<10} {result_ideal.psi:>8.4f}  "
          f"{result_ideal.x[0]:>8.4f}  {result_ideal.y[0]:>8.4f}  {result_ideal.K[0]:>8.4f}")
    print(f"  {'PR-EOS':<10} {result_pr.psi:>8.4f}  "
          f"{result_pr.x[0]:>8.4f}  {result_pr.y[0]:>8.4f}  {result_pr.K[0]:>8.4f}")

# --- Visualization ---
fig = plot_vapor_fraction_vs_pressure(
    z, comps, T_celsius, P_range=(2.0, 50.0), model="pr"
)
plt.tight_layout()
plt.savefig("methane_ethane_psi_vs_P.png", dpi=150)
plt.show()
print("\nPlot saved to methane_ethane_psi_vs_P.png")