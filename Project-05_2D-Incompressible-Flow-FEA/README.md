# Project 05 — 2D Incompressible Viscous Flow FEA Solver

## Overview

This project implements a **2D finite element solver for incompressible viscous fluid flow**, based on the steady-state and transient **Navier-Stokes equations**. It handles laminar flow problems in complex geometries using mixed velocity-pressure formulations.

The solver supports:
- **Steady-state** Stokes and Navier-Stokes flow
- **Transient** Navier-Stokes flow (time-dependent)
- **Linear Stokes** flow (low Reynolds number)
- **Non-linear** convective effects via Picard and Newton-Raphson iterations
- Triangular (T3, T6) and quadrilateral (Q4, Q9) elements with stable mixed formulations
- Dirichlet (no-slip / inlet), Neumann (outflow / traction), and periodic boundary conditions

---

## Governing Equations

The 2D incompressible Navier-Stokes equations are:

**Momentum (x and y):**
```
rho * (du/dt + u*du/dx + v*du/dy) = -dp/dx + mu*(d2u/dx2 + d2u/dy2) + bx
rho * (dv/dt + u*dv/dx + v*dv/dy) = -dp/dy + mu*(d2v/dx2 + d2v/dy2) + by
```

**Continuity (incompressibility constraint):**
```
du/dx + dv/dy = 0
```

where:
- `u`, `v` = velocity components in x and y [m/s]
- `p` = pressure [Pa]
- `rho` = fluid density [kg/m³]
- `mu` = dynamic viscosity [Pa·s]
- `bx`, `by` = body force per unit volume [N/m³]
- `Re = rho * U * L / mu` = Reynolds number

---

## Mixed Formulation (Velocity-Pressure)

The weak form yields the coupled system:

```
[K + N(u)   G ] {u}   {f}
[G^T        0 ] {p} = {0}
```

where:
- `[K]` = viscous stiffness matrix
- `[N(u)]` = convective (nonlinear) matrix
- `[G]` = pressure-velocity coupling (divergence) matrix
- `{f}` = external force / boundary traction vector

**Stable element pairs used (satisfying LBB condition):**
| Velocity Element | Pressure Element | Name |
|---|---|---|
| T6 (6-node quad. tri.) | T3 (3-node tri.) | Taylor-Hood T6/T3 |
| Q9 (9-node quad.) | Q4 (4-node quad.) | Taylor-Hood Q9/Q4 |
| Q4 + stabilization | Q4 | SUPG/PSPG stabilized |

---

## Features

| Feature | Description |
|---|---|
| Element Types | T3, T6, Q4, Q9 (mixed velocity-pressure) |
| Flow Regimes | Stokes (Re≪1), laminar Navier-Stokes |
| Analysis Type | Steady-state and transient |
| Stabilization | SUPG, PSPG, GLS for equal-order elements |
| Nonlinear Solvers | Picard iteration, Newton-Raphson |
| Time Integration | Backward Euler, Crank-Nicolson, BDF2 |
| BCs | No-slip (Dirichlet), inlet velocity, outflow, traction |
| Output | Velocity field, pressure field, streamlines, vorticity |

---

## Benchmark Problems

- **Lid-driven cavity flow** (Re = 100, 400, 1000)
- **Poiseuille flow** (analytical solution validation)
- **Flow over a backward-facing step**
- **Flow around a cylinder** (drag/lift coefficient)

---

## Input File

The solver reads problem parameters from `input-file.yaml`. Key sections include:

- **MESH**: Node coordinates and element connectivity
- **FLUID**: Density `rho`, dynamic viscosity `mu`
- **BOUNDARY_CONDITIONS**: Inlet velocity profile, no-slip walls, outflow/traction conditions
- **SOLVER**: Stokes or Navier-Stokes mode, linear/nonlinear, iteration parameters
- **TRANSIENT**: Time-stepping parameters (for unsteady flow)

---

## Usage

```bash
python solver.py input-file.yaml
```

---

## Output

- Nodal velocity field `{u(x,y), v(x,y)}`
- Nodal pressure field `p(x,y)`
- Vorticity field `omega = dv/dx - du/dy`
- Streamfunction `psi(x,y)`
- Drag and lift coefficients on immersed bodies
- Convergence history (nonlinear iterations)
- Time history at monitor points (transient only)

---

## References

- Donea, J., Huerta, A., *Finite Element Methods for Flow Problems*, Wiley
- Zienkiewicz, O.C., Taylor, R.L., Nithiarasu, P., *The Finite Element Method for Fluid Dynamics*, Elsevier
- Gresho, P.M., Sani, R.L., *Incompressible Flow and the Finite Element Method*, Wiley
- Brooks, A.N., Hughes, T.J.R., *Streamline Upwind/Petrov-Galerkin Formulations*, CMAME, 1982
