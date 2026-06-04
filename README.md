# FEM Codes — Finite Element Method Solvers in MATLAB

A collection of MATLAB-based Finite Element Method (FEM) solvers for 1D and 2D boundary value problems, beam bending, fluid flow, and structural mechanics. Each project is self-contained with its own source code, input file, and documentation.

---

## Repository Structure

```
FEM-Codes/
├── README.md                          ← You are here
├── Project-01_1D-BVP-Solver/
│   ├── README.md
│   ├── 1D-beam-solver.m
│   └── input-file.yaml
├── Project-02_2D-FEM-Solver/
│   ├── README.md
│   ├── FEM2D.m
│   └── input-file.yaml
├── Project-03_1D-Beam-Bending/
│   ├── README.md
│   ├── 1D-beam-bending.m
│   └── input-file.yaml
├── Project-04_2D-Plane-Stress/
│   ├── README.md
│   ├── FEM2D-PlaneStress.m
│   └── input-file.yaml
├── Project-05_2D-Plate-Bending/
│   ├── README.md
│   ├── FEM2D-PlateBending.m
│   └── input-file.yaml
├── Project-06_2D-Fluid-Flow/
│   ├── README.md
│   ├── FEM2D-FluidFlow.m
│   └── input-file.yaml
└── Project-07_2D-Geometric-Nonlinearity/
    ├── README.md
    ├── FEM2D-GeomNonlin.m
    └── input-file.yaml
```

---

## Projects Overview

| # | Project | Dimension | Analysis Type | PDE / Model |
|---|---------|-----------|---------------|-------------|
| 1 | [1D BVP / Non-Linear FEM Solver](./Project-01_1D-BVP-Solver/) | 1D | Linear / Picard / Newton | General 2nd-order BVP |
| 2 | [2D FEM Solver](./Project-02_2D-FEM-Solver/) | 2D | Linear / Picard / Newton | 2D diffusion / convection |
| 3 | [1D Beam Bending](./Project-03_1D-Beam-Bending/) | 1D | Linear / Nonlinear | EBT & TBT beam models |
| 4 | [2D Plane Stress](./Project-04_2D-Plane-Stress/) | 2D | Nonlinear | Plane stress elasticity |
| 5 | [2D Plate Bending](./Project-05_2D-Plate-Bending/) | 2D | Nonlinear | Kirchhoff / Mindlin plate |
| 6 | [2D Fluid Flow](./Project-06_2D-Fluid-Flow/) | 2D | Nonlinear | Viscous fluid flow (Navier-Stokes) |
| 7 | [2D Geometric Non-Linearity](./Project-07_2D-Geometric-Nonlinearity/) | 2D | Nonlinear | Geometric nonlinearity |

---

## Dependencies

- **MATLAB** (R2020a or later recommended)
- **YAML-MATLAB** library (v0.4.3) — [YAMLMatlab on GitHub](https://github.com/ewiger/yamlmatlab)
  - Add to path: `addpath(genpath('YAMLMatlab_0.4.3'))`
- **Gauss quadrature database** (`Gauss.g`) — included in each project folder

---

## How to Run

1. Clone the repository:
   ```bash
   git clone https://github.com/chinmayeechogale/FEM-Codes.git
   ```
2. Open MATLAB and navigate to the desired project folder.
3. Ensure the YAML-MATLAB library is on the path.
4. Edit the `input-file.yaml` to set problem parameters.
5. Run the main `.m` script.

---

## Solver Architecture

All solvers follow a common modular architecture:

```
 Input File (YAML)
       ↓
 Mesh Generation  →  Node coords + Connectivity
       ↓
 Shape Function Precomputation (Gauss quadrature)
       ↓
 ┌─────────────────────────────────┐
 │   Nonlinear Iteration Loop      │
 │  (Linear / Picard / Newton)     │
 │  ┌──────────────────────────┐   │
 │  │  Element-level assembly  │   │
 │  │  ELK, ELF computation    │   │
 │  └──────────────────────────┘   │
 │  Global assembly → Solve KU=F   │
 │  Apply BCs (Dirichlet/Neumann/  │
 │  Robin) → Convergence check     │
 └─────────────────────────────────┘
       ↓
   Solution U (nodal values)
```

### Boundary Condition Types Supported
| Type | Description |
|------|-------------|
| Dirichlet (Essential) | Prescribed primary variable values |
| Neumann (Natural) | Prescribed secondary variable (flux/force) |
| Robin (Mixed) | Linear combination: β₀·u + βᵤ·u = q_ref |

### Nonlinear Solvers
| Mode | Method | Notes |
|------|--------|-------|
| 0 | Linear | Single solve |
| 1 | Picard (successive substitution) | Relaxation parameter γ |
| 2 | Newton–Raphson | Tangent stiffness matrix |

---

## Author

**Chinmayee Chogale**  
Graduate Researcher — Computational Mechanics  
[GitHub](https://github.com/chinmayeechogale)

---

## License

This repository is for academic and educational use.
