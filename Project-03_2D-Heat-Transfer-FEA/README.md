# Project 03 — 2D Heat Transfer FEA Solver

## Overview

This project implements a **2D finite element solver for steady-state and transient heat transfer** problems. It solves the heat conduction equation over 2D domains with support for anisotropic conductivity, convective boundary conditions, and internal heat generation.

The solver supports:
- **Steady-state** heat conduction
- **Transient** heat conduction (time-dependent)
- **Linear and non-linear** thermal conductivity (temperature-dependent)
- CST (3-node triangle) and Q4 (4-node quadrilateral) elements
- Dirichlet (prescribed temperature), Neumann (heat flux), and Robin (convection) boundary conditions

---

## Governing Equation

The 2D steady-state heat conduction equation is:

```
-d/dx [ kx(x,y,T) * dT/dx ] - d/dy [ ky(x,y,T) * dT/dy ] + Q(x,y,T) = 0
```

For transient analysis:

```
rho * cp * dT/dt - d/dx [ kx * dT/dx ] - d/dy [ ky * dT/dy ] + Q = 0
```

where:
- `kx`, `ky` = thermal conductivity in x and y directions [W/m·K]
- `Q` = volumetric heat generation rate [W/m³]
- `rho` = density [kg/m³]
- `cp` = specific heat capacity [J/kg·K]
- `T` = temperature [K or °C]

---

## Boundary Condition Types

| Type | Description | Math Form |
|---|---|---|
| Dirichlet | Prescribed temperature | T = T_prescribed |
| Neumann | Prescribed heat flux | -k * dT/dn = q_n |
| Robin | Convection (Newton's law) | -k * dT/dn = h*(T - T_inf) |

---

## Features

| Feature | Description |
|---|---|
| Element Types | CST (3-node triangle), Q4 (4-node quad) |
| Analysis Type | Steady-state and transient |
| Conductivity | Isotropic or anisotropic; constant or T-dependent |
| Solvers | Direct, Picard iteration, Newton-Raphson |
| Time Integration | Backward Euler (implicit), Crank-Nicolson |
| Output | Temperature field, heat flux, gradient maps |

---

## Input File

The solver reads problem parameters from `input-file.yaml`. Key sections include:

- **MESH**: Node coordinates and element connectivity
- **MATERIAL**: Thermal conductivity `k`, density `rho`, specific heat `cp`
- **BOUNDARY_CONDITIONS**: Prescribed temperatures, heat fluxes, convection coefficients
- **SOURCE**: Internal heat generation `Q`
- **SOLVER**: Steady-state or transient, linear or nonlinear, time-stepping parameters

---

## Usage

```bash
python solver.py input-file.yaml
```

---

## Output

- Nodal temperature distribution `T(x,y)`
- Element heat flux vectors `{q_x, q_y}`
- Temperature gradient field
- Convergence history (nonlinear/transient analysis)
- Time history at selected nodes (transient only)

---

## References

- Lewis, R.W., Nithiarasu, P., Seetharamu, K.N., *Fundamentals of the Finite Element Method for Heat and Fluid Flow*, Wiley
- Reddy, J.N., Gartling, D.K., *The Finite Element Method in Heat Transfer and Fluid Dynamics*, CRC Press
- Bathe, K.J., *Finite Element Procedures*, Prentice Hall
