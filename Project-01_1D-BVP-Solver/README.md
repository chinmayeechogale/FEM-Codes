# Project 01 — 1D Boundary Value Problem (BVP) / Non-Linear FEM Solver

## Overview

This project implements a general-purpose **1D finite element solver** for second-order boundary value problems (BVPs) of the form:

```
-d/dx [ A(x,u,u_x) * du/dx ] + B(x,u,u_x) * du/dx + C(x,u,u_x) * u = F(x)
```

The solver supports **linear**, **Picard (successive substitution)**, and **Newton-Raphson** nonlinear solution schemes, and handles all three standard boundary condition types.

---

## Files

| File | Description |
|------|-------------|
| `1D-beam-solver.m` | Main MATLAB solver script |
| `input-file.yaml` | YAML input configuration file |

---

## Governing PDE

The weak form of the general 1D BVP solved here is:

```
A(x, u, u_x) = AX0 + AX1*x + AU1*u + AU2*u^2 + AUX1*u_x + AUX2*u_x^2
B(x, u, u_x) = BX0 + BX1*x + BU1*u + BU2*u^2 + BUX1*u_x + BUX2*u_x^2
C(x, u, u_x) = CX0 + CX1*x + CU1*u + CU2*u^2 + CUX1*u_x + CUX2*u_x^2
F(x)         = FX0 + FX1*x + FX2*x^2
```

All coefficients are specified in the YAML input file.

---

## Input File Parameters

| Parameter | Description |
|-----------|-------------|
| `X0` | Left boundary x-coordinate |
| `L` | Domain length |
| `A`, `B`, `C` | PDE coefficient arrays `[X0, X1, U1, U2, UX1, UX2]` |
| `F` | Source term coefficients `[FX0, FX1, FX2]` |
| `P` | Element polynomial order |
| `NEM` | Number of elements |
| `NGP` | Number of Gauss quadrature points |
| `NSPV` | Number of Dirichlet (essential) BCs |
| `NSSV` | Number of Neumann (natural) BCs |
| `NSMB` | Number of Robin (mixed) BCs |
| `NONLIN` | Solver mode: `0`=Linear, `1`=Picard, `2`=Newton |
| `ITERMAX` | Maximum nonlinear iterations |
| `EPSILON` | Convergence tolerance |
| `GAMMA` | Relaxation parameter (Picard) |

---

## Features

- Arbitrary polynomial order Lagrange elements (P=1, 2, 3, ...)
- Automatic differentiation via dual numbers for shape functions
- Three nonlinear solution strategies
- Robin (mixed), Neumann, and Dirichlet boundary conditions
- Gauss quadrature integration from external database

---

## Usage

```matlab
% Set input file in script header:
INPFILE = 'input-file.yaml';

% Run:
1D-beam-solver
```

### Dependencies
- YAML-MATLAB v0.4.3 (add to path)
- `Gauss.g` quadrature database

---

## Example Problem (from input file)

The default input file solves a nonlinear ODE with:
- Domain: [0, 1]
- `A = -1` (diffusion coefficient)
- `B = -2u` (nonlinear convection)
- Neumann BC at left end, Robin BC at right end
- Picard iteration with 100 max iterations, tolerance 1e-3

---

## Algorithm

1. Read YAML input and build mesh via `MESH1D`
2. Precompute Lagrange shape functions and derivatives at Gauss points
3. Loop over nonlinear iterations:
   - Assemble element stiffness matrix `ELK` and force vector `ELF`
   - Assemble global `GLK`, `GLF`
   - Apply boundary conditions
   - Solve `GLK * U = GLF`
   - Check convergence
4. Output nodal solution `GCU`
