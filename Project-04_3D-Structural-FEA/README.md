# Project 04 — 3D Structural FEA Solver

## Overview

This project implements a **3D finite element solver for structural mechanics**, capable of analyzing solid bodies under static and dynamic loading. It solves the full 3D elasticity equations with support for linear and non-linear material behavior, large deformations, and multiple element types.

The solver supports:
- **Linear static** analysis
- **Non-linear static** analysis (geometric and/or material nonlinearity)
- **Modal analysis** (natural frequencies and mode shapes)
- **Transient dynamic** analysis (implicit and explicit time integration)
- Tetrahedral (TET4, TET10) and hexahedral (HEX8, HEX20) elements
- Dirichlet (prescribed displacement), Neumann (traction/force), and Robin (elastic foundation) boundary conditions

---

## Governing Equations

The 3D equilibrium equation in vector form is:

```
div(sigma) + b = rho * u_tt
```

For static analysis (`u_tt = 0`):

```
div(sigma) + b = 0
```

The generalized constitutive relation is:

```
{sigma} = [D] {epsilon}
```

where `[D]` is the 6x6 isotropic elasticity matrix:

```
         E*(1-nu)                 nu                 nu
[D] = -------------- * |  1      nu/(1-nu)   nu/(1-nu)   0   0   0  |
      (1+nu)*(1-2nu)   | nu/(1-nu)   1      nu/(1-nu)   0   0   0  |
                       | nu/(1-nu)  nu/(1-nu)   1        0   0   0  |
                       |  0          0          0   (1-2nu)/(2(1-nu)) 0  0 |
                       |  0          0          0    0  (1-2nu)/(2(1-nu)) 0 |
                       |  0          0          0    0   0  (1-2nu)/(2(1-nu))|
```

The strain-displacement relation:

```
{epsilon} = [B] {u}    =>    {epsilon_xx, epsilon_yy, epsilon_zz, gamma_xy, gamma_yz, gamma_xz}
```

---

## Features

| Feature | Description |
|---|---|
| Element Types | TET4, TET10, HEX8, HEX20 |
| Analysis Types | Linear static, nonlinear static, modal, transient dynamic |
| Material Models | Linear elastic, elastoplastic (von Mises, Drucker-Prager) |
| Solvers | Direct (LU), iterative (PCG), Newton-Raphson, Arc-Length |
| Time Integration | Newmark-beta, HHT-alpha (implicit); Central Difference (explicit) |
| BCs | Dirichlet, Neumann (traction/point force), Pressure, Robin |
| Output | Displacements, stresses, strains, reaction forces, mode shapes |

---

## Boundary Condition Types

| Type | Description | Math Form |
|---|---|---|
| Dirichlet | Prescribed displacement | u = u_prescribed |
| Neumann | Applied traction or point force | sigma*n = t |
| Pressure | Normal surface pressure | sigma*n = -p*n |
| Robin | Elastic foundation (Winkler) | sigma*n = -k_f * u |

---

## Input File

The solver reads problem parameters from `input-file.yaml`. Key sections include:

- **MESH**: 3D node coordinates (x, y, z) and element connectivity
- **MATERIAL**: Elastic constants `E`, `nu`, density `rho`, and optional plasticity parameters
- **BOUNDARY_CONDITIONS**: Prescribed displacements, tractions, pressures, and body forces
- **SOLVER**: Analysis type, linear/nonlinear mode, iteration parameters
- **DYNAMICS**: Mass matrix type, time integration scheme, time-stepping parameters (for transient/modal)

---

## Usage

```bash
python solver.py input-file.yaml
```

---

## Output

- Nodal displacement vector `{u_x, u_y, u_z}`
- Element stress tensor `{sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_xz}`
- von Mises stress field
- Principal stresses and directions
- Reaction forces at constrained nodes
- Natural frequencies and mode shapes (modal analysis)
- Time history of displacements/stresses (transient analysis)

---

## References

- Bathe, K.J., *Finite Element Procedures*, Prentice Hall
- Hughes, T.J.R., *The Finite Element Method*, Dover Publications
- Zienkiewicz, O.C., Taylor, R.L., *The Finite Element Method for Solid and Structural Mechanics*, Elsevier
- Belytschko, T., Liu, W.K., Moran, B., *Nonlinear Finite Elements for Continua and Structures*, Wiley
