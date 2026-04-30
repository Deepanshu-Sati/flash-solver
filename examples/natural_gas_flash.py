"""
Example: Natural Gas Mixture Flash (5 components)
===================================================
Simulates a typical wellhead natural gas flash separation.

Feed composition (approximate Gulf Coast gas):
    methane:  70%
    ethane:   15%
    propane:   8%
    n-butane:  4%
    n-pentane: 3%

Operating conditions: T = 40°C, P = 70 bar (separator)
"""

import numpy as np
import pandas as pd
import sys
sys.path.insert(0, "..")

import matplotlib.pyplot as plt
from flash.flash_calculation import isothermal_flash_pr
from visualization.phase_diagrams import plot_composition_profile

db = pd.read_csv("data/component_database.csv", index_col="component")

components = ["methane", "ethane", "propane", "n-butane", "n-pentane"]
comps = db.loc[components]

z = np.array([0.70, 0.15, 0.08, 0.04, 0.03])
T_celsius = 40.0
P_bar = 70.0

result = isothermal_flash_pr(z, T_celsius, P_bar, comps)

print("=" * 60)
print("  Natural Gas Flash Separator — PR-EOS")
print(f"  T = {T_celsius}°C, P = {P_bar} bar")
print(f"  Phase: {result.phase} | ψ = {result.psi:.4f}")
print("=" * 60)
print(f"  {'Component':<12} {'Feed z':>8}  {'Liquid x':>8}  {'Vapor y':>8}  {'K':>8}")
for i, comp in enumerate(components):
    print(f"  {comp:<12} {z[i]:>8.4f}  {result.x[i]:>8.4f}  "
          f"{result.y[i]:>8.4f}  {result.K[i]:>8.4f}")

fig = plot_composition_profile(result)
plt.tight_layout()
plt.savefig("natural_gas_flash_compositions.png", dpi=150)
plt.show()