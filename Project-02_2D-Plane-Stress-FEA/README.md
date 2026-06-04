# Project 02 — 2D Plane Stress Non-Linear FEA Solver

## Overview

This project implements a **2D plane stress finite element analysis (FEA) solver** capable of handling both linear and non-linear material behavior. It solves the governing equations of 2D elasticity under the plane stress assumption (thin structures where out-of-plane stress is zero).

The solver supports:
- Linear elastic analysis
- Non-linear analysis using **Picard (successive substitution)** and **Newton-Raphson** iteration schemes
- CST (Constant Strain Triangle) and/or quadrilateral elements
- Multiple boundary condition types: Dirichlet, Neumann, and Robin

---

## Governing Equations

Under plane stress, the equilibrium equations are:

```
∂σ_xx/∂x + ∂σ_xy/∂y + b_x = 0
∂σ_xy/∂x + ∂σ_yy/∂y + b_y = 0
```

where `σ` is the stress tensor and `b` is the body force vector.

The constitutive relation (for isotropic linear elastic material) is:

```
{σ} = [D] {ε}
```

where `[D]` is the plane stress elasticity matrix:

```
[D] = E/(1-ν²) * | 1   ν   0         |
                  | ν   1   0         |
                  | 0   0   (1-ν)/2   |
```

---

## Features

| Feature | Description |
|---|---|
| Element Types | CST (3-node triangle), Q4 (4-node quad) |
| Material Models | Linear elastic, nonlinear (user-defined) |
| Solvers | Direct, Picard iteration, Newton-Raphson |
| BCs | Dirichlet (displacement), Neumann (traction), Robin |
| Output | Displacement field, stress/strain distributions |

---

## Input File

The solver reads problem parameters from `input-file.yaml`. Key sections include:

- **MESH**: Node coordinates and element connectivity
- **MATERIAL**: Young's modulus `E`, Poisson's ratio `ν`, and optional nonlinear parameters
- **BOUNDARY_CONDITIONS**: Prescribed displacements, tractions, and body forces
- **SOLVER**: Linear or nonlinear mode, iteration parameters

---

## Usage

```bash
python solver.py input-file.yaml
```

---

## Output

- Nodal displacement vector `{u}`
- Element stress and strain tensors
- Reaction forces at constrained nodes
- Convergence history (for nonlinear analysis)

---

## References

- Hughes, T.J.R., *The Finite Element Method*, Dover Publications
- Bathe, K.J., *Finite Element Procedures*, Prentice Hall
- Reddy, J.N., *An Introduction to Nonlinear Finite Element Analysis*, Oxford University Press
