# Project 06 — Fluid-Structure Interaction (FSI) FEA Solver

## Overview

This project implements a **coupled Fluid-Structure Interaction (FSI) finite element solver** that simulates the two-way coupling between an incompressible viscous fluid and a deformable elastic structure. The solver handles problems where fluid forces cause structural deformation, and structural motion in turn modifies the fluid domain.

The solver supports:
- **Monolithic** (fully coupled, simultaneous) and **partitioned** (staggered) FSI solution strategies
- **Arbitrary Lagrangian-Eulerian (ALE)** formulation for the fluid domain
- **Linear and geometrically nonlinear** structural response
- Steady-state and transient FSI problems
- 2D and quasi-3D configurations
- Dynamic mesh update via **mesh motion** solvers (Laplacian smoothing, radial basis functions)

---

## Governing Equations

### Fluid Domain (ALE Navier-Stokes)

In the ALE frame, the incompressible Navier-Stokes equations become:

```
rho_f * [du/dt + (u - w)*grad(u)] = -grad(p) + mu*lap(u) + b_f
div(u) = 0
```

where `w` is the mesh (ALE frame) velocity and `u - w` is the convective velocity relative to the moving mesh.

### Structure Domain (Lagrangian Elastodynamics)

```
rho_s * d2d/dt2 - div(sigma_s) = b_s
```

where `d` is the structural displacement, `sigma_s` is the Cauchy stress tensor, and `b_s` is the body force.

### Interface Coupling Conditions

At the fluid-structure interface `Gamma_FSI`:

```
  Kinematic: u_f = dd/dt        (velocity continuity)
  Dynamic:   sigma_f * n = sigma_s * n   (traction equilibrium)
  Geometric: x_f = x_s          (position continuity)
```

---

## Solution Strategies

| Strategy | Description | Stability | Cost |
|---|---|---|---|
| Monolithic | Solve fluid + structure simultaneously | Unconditionally stable | High (large system) |
| Partitioned (strong) | Iterate fluid-structure until convergence | Stable (with subiterations) | Moderate |
| Partitioned (weak) | One fluid solve + one structure solve per step | Conditionally stable | Low |

---

## Features

| Feature | Description |
|---|---|
| FSI Formulation | ALE-based monolithic and partitioned |
| Fluid Elements | Q4/Q9, T3/T6 (stabilized or Taylor-Hood) |
| Structure Elements | Q4, T3, beam elements (Euler-Bernoulli, Timoshenko) |
| Mesh Motion | Laplacian smoothing, RBF interpolation |
| Nonlinear Solvers | Newton-Raphson (monolithic), Aitken acceleration (partitioned) |
| Time Integration | Generalized-alpha, Newmark-beta, BDF2 |
| BCs | Inlet/outlet flow, no-slip, fixed/free structure ends, FSI interface |
| Output | Fluid velocity/pressure, structural displacement/stress, interface forces |

---

## Benchmark Problems

- **Flexible beam in channel flow** (Turek & Hron FSI benchmark)
- **Flow-induced vibration of a cylinder with elastic tail**
- **Pressure-driven inflation of an elastic tube**
- **Flag-in-wind flutter simulation**

---

## Input File

The solver reads problem parameters from `input-file.yaml`. Key sections include:

- **FLUID_MESH**: Fluid domain node coordinates and element connectivity
- **STRUCTURE_MESH**: Structural domain node coordinates and element connectivity
- **FSI_INTERFACE**: Node/element pairs defining the fluid-structure coupling boundary
- **FLUID**: Density, viscosity, flow model parameters
- **STRUCTURE**: Material properties (E, nu, rho), element type, nonlinear flag
- **SOLVER**: FSI coupling strategy (monolithic/partitioned), nonlinear iteration parameters
- **MESH_MOTION**: Mesh update method and parameters
- **TRANSIENT**: Time integration scheme, time step, total simulation time

---

## Usage

```bash
python fsi_solver.py input-file.yaml
```

---

## Output

- Fluid velocity field `{u(x,y,t), v(x,y,t)}`
- Fluid pressure field `p(x,y,t)`
- Structural displacement field `{d_x(x,y,t), d_y(x,y,t)}`
- Structural stress and strain fields
- Interface traction forces
- Drag and lift coefficients on the structure
- Tip displacement time history (for beam/flag problems)
- Convergence history (FSI subiterations)

---

## References

- Donea, J., Huerta, A., *Finite Element Methods for Flow Problems*, Wiley
- Bazilevs, Y., Takizawa, K., Tezduyar, T.E., *Computational Fluid-Structure Interaction*, Wiley
- Hou, G., Wang, J., Layton, A., *Numerical Methods for Fluid-Structure Interaction*, Commun. Comput. Phys., 2012
- Turek, S., Hron, J., *Proposal for Numerical Benchmarking of Fluid-Structure Interaction*, Springer, 2006
- Farhat, C., Lesoinne, M., *Two Efficient Staggered Algorithms for Serial and Parallel Solution of FSI Problems*, CMAME, 2000
