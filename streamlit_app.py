"""
Multicomponent Flash Solver — Streamlit Web Interface
======================================================
Interactive web app for running flash calculations and
visualizing phase behavior.

Run with:
    streamlit run app/streamlit_app.py

Author: Deepanshu Sati
"""

import sys
sys.path.insert(0, "..")

import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from flash.flash_calculation import isothermal_flash_ideal, isothermal_flash_pr
from visualization.phase_diagrams import (
    plot_composition_profile,
    plot_vapor_fraction_vs_pressure,
    plot_Pxy_binary
)


# ─── Page Configuration ────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Flash Solver Pro",
    page_icon="⚗️",
    layout="wide"
)

# ─── Load Database ─────────────────────────────────────────────────────────────
@st.cache_data
def load_database():
    return pd.read_csv("data/component_database.csv", index_col="component")

db = load_database()
available_components = list(db.index)


# ─── Sidebar: User Inputs ──────────────────────────────────────────────────────
st.sidebar.title("⚗️ Flash Solver")
st.sidebar.markdown("### Operating Conditions")

T_celsius = st.sidebar.slider(
    "Temperature (°C)", min_value=-100.0, max_value=300.0, value=25.0, step=1.0
)
P_bar = st.sidebar.slider(
    "Pressure (bar)", min_value=0.1, max_value=200.0, value=5.0, step=0.1
)

st.sidebar.markdown("### Components & Feed Composition")
selected_components = st.sidebar.multiselect(
    "Select components",
    available_components,
    default=["methane", "ethane", "propane"]
)

model = st.sidebar.radio("Thermodynamic model", ["Ideal (Raoult's Law)", "Peng-Robinson EOS"])

# ─── Main Panel ────────────────────────────────────────────────────────────────
st.title("⚗️ Multicomponent Flash Solver")
st.markdown(
    "A **professional isothermal flash calculator** implementing both "
    "Raoult's Law and the Peng–Robinson Equation of State."
)

if len(selected_components) < 2:
    st.warning("Please select at least 2 components from the sidebar.")
    st.stop()

# Feed composition inputs
st.markdown("### Feed Composition (mole fractions)")
cols = st.columns(len(selected_components))
z_values = []

for i, (col, comp) in enumerate(zip(cols, selected_components)):
    with col:
        default_z = round(1.0 / len(selected_components), 3)
        z_i = st.number_input(
            comp, min_value=0.001, max_value=0.999,
            value=default_z, step=0.01, key=f"z_{i}"
        )
        z_values.append(z_i)

z = np.array(z_values)
z = z / z.sum()  # Normalize

st.info(f"Normalized feed: " +
        ", ".join([f"{selected_components[i]}: {z[i]:.3f}" for i in range(len(z))]))

# ─── Run Flash ─────────────────────────────────────────────────────────────────
if st.button("🚀 Run Flash Calculation", type="primary"):
    params_df = db.loc[selected_components]
    use_pr = "Peng" in model

    with st.spinner("Solving flash equations..."):
        try:
            if use_pr:
                result = isothermal_flash_pr(z, T_celsius, P_bar, params_df)
            else:
                result = isothermal_flash_ideal(z, T_celsius, P_bar, params_df)
        except Exception as e:
            st.error(f"Flash calculation failed: {e}")
            st.stop()

    # ─── Results Display ───────────────────────────────────────────────────────
    col1, col2, col3 = st.columns(3)
    col1.metric("Phase State", result.phase.replace("-", " ").title())
    col2.metric("Vapor Fraction ψ", f"{result.psi:.4f}")
    if use_pr:
        col3.metric("Iterations", result.iterations)

    # Composition table
    st.markdown("### Phase Compositions")
    z_recon = result.psi * result.y + (1 - result.psi) * result.x
    comp_data = {
        "Component": result.components,
        "Feed z_i": z_recon.round(4),
        "Liquid x_i": result.x.round(4),
        "Vapor y_i": result.y.round(4),
        "K-value": result.K.round(4),
    }
    st.dataframe(pd.DataFrame(comp_data), use_container_width=True)

    # Composition bar chart
    st.markdown("### Composition Distribution")
    fig = plot_composition_profile(result)
    st.pyplot(fig)
    plt.close()

    # Vapor fraction sensitivity plot
    st.markdown("### Vapor Fraction vs Pressure")
    P_low = max(0.1, P_bar * 0.1)
    P_high = P_bar * 5.0
    fig2 = plot_vapor_fraction_vs_pressure(
        z, params_df, T_celsius,
        P_range=(P_low, P_high),
        model="pr" if use_pr else "ideal"
    )
    st.pyplot(fig2)
    plt.close()

# ─── Educational Footer ────────────────────────────────────────────────────────
with st.expander("📚 Theory — Rachford-Rice Equation"):
    st.latex(r"f(\psi) = \sum_{i=1}^{n} \frac{z_i(K_i - 1)}{1 + \psi(K_i - 1)} = 0")
    st.markdown("""
    **Where:**
    - ψ (psi) = vapor fraction
    - z_i = feed mole fraction of component i
    - K_i = equilibrium ratio (y_i/x_i)
    
    Solved numerically with **Brent's method** over the guaranteed bracket.
    """)

with st.expander("📚 Theory — Peng-Robinson EOS"):
    st.latex(r"P = \frac{RT}{V-b} - \frac{a(T)}{V(V+b)+b(V-b)}")
    st.markdown("""
    K-values computed from fugacity coefficient ratio:
    """)
    st.latex(r"K_i = \frac{\hat{\phi}_i^L}{\hat{\phi}_i^V}")