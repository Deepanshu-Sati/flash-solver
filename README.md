# Multicomponent Flash Solver

A **thermodynamic flash calculation simulator** implemented in Python that predicts vapor–liquid equilibrium (VLE) of multicomponent mixtures using both **ideal models** and the **Peng–Robinson equation of state**.

This project demonstrates the core algorithms used in chemical process simulators such as **Aspen Plus** and **HYSYS**.

---

# Overview

In many chemical processes, mixtures of components split into **vapor and liquid phases** depending on temperature and pressure.

This project computes:

- Vapor fraction of the mixture
- Composition of liquid phase
- Composition of vapor phase

for a given **temperature, pressure, and feed composition**.

The solver implements the **Rachford–Rice equation** together with thermodynamic models to determine phase equilibrium.

---

# Features

- Multicomponent **isothermal flash calculations**
- Ideal vapor-liquid equilibrium using **Raoult's Law**
- Real thermodynamics using **Peng–Robinson EOS**
- Fugacity coefficient calculations
- Numerical solution of the **Rachford–Rice equation**
- Visualization of phase behavior
- Command line interface
- Interactive **Streamlit web application**

---

# Project Structure

```

flash-solver/
│
├── flash/
│   ├── flash_calculation.py
│   └── rachford_rice.py
│
├── thermodynamics/
│   ├── antoine.py
│   ├── raoults_law.py
│   ├── peng_robinson.py
│   └── fugacity.py
│
├── visualization/
│   └── phase_diagrams.py
│
├── utils/
│   └── numerical_methods.py
│
├── examples/
│   ├── methane_ethane_flash.py
│   ├── natural_gas_flash.py
│   └── propane_butane_flash.py
│
├── data/
│   └── component_database.csv
│
├── streamlit_app.py
├── main.py
├── requirements.txt
└── README.md

```

---

# Theory

## Phase Equilibrium Condition

At equilibrium:

```

f_i^L = f_i^V

```

which leads to the equilibrium ratio:

```

K_i = y_i / x_i

```

where

- `x_i` = liquid mole fraction  
- `y_i` = vapor mole fraction

---

## Rachford–Rice Equation

The vapor fraction `ψ` is obtained by solving:

```

f(ψ) = Σ [ z_i (K_i − 1) / (1 + ψ (K_i − 1)) ] = 0

```

where

- `z_i` = feed composition
- `K_i` = equilibrium constant

---

## Thermodynamic Models

### Ideal Model

Raoult's Law:

```

K_i = P_sat,i / P

```

where saturation pressure is computed using the **Antoine equation**.

---

### Real Model

The project implements the **Peng–Robinson Equation of State**:

```

P = RT/(V − b) − a(T)/(V(V + b) + b(V − b))

````

Fugacity coefficients are used to compute real K-values.

---

# Installation

Clone the repository:

```bash
git clone https://github.com/yourusername/flash-solver.git
cd flash-solver
````

Install dependencies:

```bash
pip install -r requirements.txt
```

---

# Running the Flash Solver

## Command Line Interface

Example flash calculation:

```bash
python main.py --components methane ethane propane \
--feed 0.5 0.3 0.2 \
--temperature 25 \
--pressure 10 \
--model pr \
--plot
```

The output includes:

* phase state
* vapor fraction
* liquid composition
* vapor composition

---

# Running the Streamlit App

Launch the interactive interface:

```bash
streamlit run streamlit_app.py
```

The web interface allows users to:

* choose components
* adjust temperature and pressure
* run flash calculations
* visualize phase compositions and phase behavior

---

# Example Simulations

Run the example scripts:

```bash
python -m examples.natural_gas_flash
```

or

```bash
python -m examples.methane_ethane_flash
```

These simulate realistic mixtures such as **natural gas separation**.

---

# Applications

Flash calculations are widely used in:

* Natural gas processing
* Petroleum refining
* Distillation column design
* Separation process simulation
* Process simulators (Aspen Plus, HYSYS)

---

# Future Improvements

Possible extensions include:

* Phase envelope calculation
* PT flash (energy balance)
* Distillation column simulation
* Support for additional equations of state (SRK)
* Thermodynamic parameter estimation

---

# License

This project is licensed under the **MIT License**.

---

# Author

Developed as a chemical engineering simulation project.

````
