# Project 04 — 2D Transient FEM Solver

## Overview

This project implements a **2D finite element solver for transient (time-dependent) scalar field problems**. It solves a general second-order parabolic/hyperbolic PDE using the Galerkin weighted residual method, covering problems such as:
- **Transient heat conduction** (parabolic: first-order in time)
- **Wave propagation** (hyperbolic: second-order in time)
- **Coupled diffusion-reaction** problems

The solver supports:
- **First-order** (parabolic) and **second-order** (hyperbolic) time dependence, selected via `ITEM`
- **α-family time integration scheme** (α=0: Forward Euler; α=0.5: Crank-Nicolson; α=1: Backward Euler)
- **Newmark-β method** for second-order problems, parameterized by `NEW` (β)
- **Linear, Picard, and Newton-Raphson** nonlinear solvers
- **Dirichlet, Neumann, and Robin** boundary conditions
- **Q4/Q8/Q9** Lagrangian quadrilateral elements

---

## Files

| File | Description |
|------|-------------|
| `2D-Transient.m` | Main MATLAB solver script |
| `input-file.yaml` | YAML input file (sample: 2D transient diffusion on a unit square) |
| `GAUSS.g` | Gauss quadrature point/weight database |
| `ReadYaml.m` | YAML reader (YAML-MATLAB v0.4.3 dependency) |

---

## Governing Equation

The general 2D transient PDE solved is:

```
A0*u + C1(x,y)*∂u/∂t + C2(x,y)*∂²u/∂t²
   - ∂/∂x[A11(x,y,u)*∂u/∂x]
   - ∂/∂y[A22(x,y,u)*∂u/∂y] = F(x,y)
```

where:
- `u` = primary field variable (e.g., temperature, pressure)
- `A11`, `A22` = diffusion/conductivity coefficients (can be nonlinear functions of x, y, u, ∇u)
- `A0` = reaction coefficient
- `C1` = damping/capacitance coefficient (first-order time term)
- `C2` = inertia coefficient (second-order time term; set to 0 for parabolic problems)
- `F` = source/body force term

### PDE Coefficients (general polynomial form)

```
A11 = A10 + A1X*x + A1Y*y + A1U*u + A1UX*du/dx + A1UY*du/dy
A22 = A20 + A2X*x + A2Y*y + A2U*u + A2UX*du/dx + A2UY*du/dy
C1  = C10 + C1X1*x + C1Y1*y
C2  = C20 + C2X1*x + C2Y1*y
F   = FX0 + FX1*x + FY1*y
```

### Time Integration

**For first-order problems (`ITEM = 1`, parabolic):**

The α-family scheme:
```
[C/Δt + α*K] * U^(n+1) = [C/Δt - (1-α)*K] * U^n + F
```

The effective stiffness and force are assembled using:
```
DATA.A1 = α * Δt
DATA.A2 = (1 - α) * Δt
```

**For second-order problems (`ITEM = 1` with Newmark-β activated):**

Newmark-β parameters:
```
DATA.A3 = 2/(β*Δt²)    (acceleration coefficient)
DATA.A4 = DATA.A3 * Δt
DATA.A5 = 1/β - 1
DATA.A6 = 2*α*β/Δt
DATA.A7 = 2*α*β - 1
DATA.A8 = Δt*α*(β - 1)
```

Update relations each time step:
```
a^(n+1) = A3*(U^(n+1) - U^n) - A4*v^n - A5*a^n
v^(n+1) = v^n + A1*a^(n+1) + A2*a^n
```

---

## Element & Mesh

| Feature | Description |
|---------|-------------|
| Element family | 2D Lagrangian quadrilateral |
| Supported NPE | 4 (Q4), 8 (Q8), 9 (Q9) |
| Default (sample input) | 9-node biquadratic Q9 |
| Mesh generation | Built-in `MESH2DR`: structured rectangular grid |
| Shape functions | `INTERPLN2D` — Lagrangian, mapped from reference [-1,1]² |
| Integration | Gauss quadrature; `NGP` points per direction |
| DOFs per node | 1 (scalar field: temperature, pressure, etc.) |

---

## Nonlinear Solvers

| `NONLIN` | Mode | Description |
|----------|------|-------------|
| `0` | Linear | Direct solve at each time step |
| `1` | Picard (direct iteration) | Iterative substitution within each time step |
| `2` | Newton-Raphson | Tangent-based correction; uses `TAN` matrix assembled in `ELEMAT2D` |

---

## Boundary Conditions

| Type | Parameters | Description |
|------|-----------|-------------|
| Dirichlet (essential) | `NSPV`, `ISPV`, `VSPV` | Prescribed nodal field values |
| Neumann (natural) | `NSSV`, `ISSV`, `VSSV` | Prescribed normal flux / heat flux |
| Robin (mixed) | `NSMB`, `ISMB`, `UREF`, `BETA0`, `BETAU` | Convective / radiation boundary |

---

## Input File (`input-file.yaml`)

```yaml
# DOMAIN DATA
X0: 0   # X-coordinate of bottom-left corner
Y0: 0   # Y-coordinate of bottom-left corner

# PDE COEFFICIENTS
# A11 = [A10, A1X, A1Y, A1U, A1UX, A1UY]
A11: [1,0,0,0,0,0]   # Diffusion in X (isotropic, constant k=1)
A22: [1,0,0,0,0,0]   # Diffusion in Y
A0:  0                # Reaction coefficient
C1:  [1,0,0]          # Capacitance/damping [C10, C1X1, C1Y1]
C2:  [0,0,0]          # Inertia [C20, C2X1, C2Y1] (0 = parabolic/heat)
F:   [1,0,0]          # Source term [FX0, FX1, FY1]

# MESH DATA
NPE: 9                        # Nodes per element (4, 8, or 9)
DX:  [0.25,0.25,0.25,0.25]   # Element lengths along X (4 elements)
DY:  [0.25,0.25,0.25,0.25]   # Element lengths along Y (4 elements)
NGP: 5                        # Gauss points per direction

# BOUNDARY CONDITIONS
NSPV: 17   # Number of Dirichlet BCs
NSSV: 0    # Number of Neumann BCs
NSMB: 0    # Number of Robin BCs

# Dirichlet: prescribe u=0 on right edge and top edge
ISPV:
  NODES: [9,18,27,36,45,54,63,72,81,73,74,75,76,77,78,79,80]
  DOFS:  [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
VSPV: [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

# NONLINEAR SOLVER
NONLIN:  0      # 0 = Linear; 1 = Picard; 2 = Newton
ITERMAX: 25     # Max iterations per time step
EPSILON: 1e-3   # Convergence tolerance
GAMMA:   1.0    # Relaxation parameter

# TIME INTEGRATION PARAMETERS
ITEM:     1       # 1 = transient (first-order); 0 = steady-state
TIMESTEP: 0.05   # Time step size (Δt)
NTIME:    25      # Number of time steps

# STABILITY PARAMETERS (α-family / Newmark-β)
ALPHA: 1            # α parameter (0=Forward Euler, 0.5=Crank-Nicolson, 1=Backward Euler)
NEW:   0.333333333  # β parameter for Newmark method

# POST-PROCESSING
POINT: 1   # Node index to print solution history
```

---

## Usage

1. **Install dependencies:** Add `YAML-MATLAB v0.4.3` to the MATLAB path.
2. **Prepare the input file:** Edit `input-file.yaml` with domain geometry, PDE coefficients, time parameters, BCs, and solver settings.
3. **Set `ITEM`:** `1` for transient, `0` for steady-state.
4. **Set `ALPHA` and `NEW`:** Control time integration accuracy and stability.
5. **Ensure support files:** `GAUSS.g` must be in the working directory.
6. **Run the solver in MATLAB:**
   ```matlab
   run('2D-Transient')
   ```
7. **Review output:** The solver prints the solution at node `POINT` at every time step. The full solution history is stored in `TIMEHIST`.

---

## Key Functions

| Function | Description |
|----------|-------------|
| `ELEMAT2D` | Computes element stiffness (K), mass (M), damping (C) matrices and force vector; assembles effective system `ELKEFF`, `ELFEFF` for the chosen time integration scheme |
| `MESH2DR` | Generates structured 2D mesh: node coordinates and element connectivity |
| `INTERPLN2D` | Evaluates Lagrangian shape functions and derivatives for Q4/Q8/Q9 elements |
| `PRECOMPUTE_SF2D` | Pre-computes shape functions at all Gauss points |
| `GAUSS` | Reads Gauss quadrature data from `GAUSS.g` |

---

## Sample Problem

The default input file solves **2D transient heat conduction** on a unit square [0,1]×[0,1]:
- Isotropic conductivity `k = 1` (A11 = A22 = 1)
- Uniform heat source `F = 1`
- Capacitance `C1 = 1` (first-order time dependence)
- Zero temperature on the right edge (x=1) and top edge (y=1)
- **Backward Euler** time integration (α=1), `Δt = 0.05`, 25 time steps
- 4×4 mesh of Q9 elements (9-node biquadratic), `NGP=5`
- Linear solver (`NONLIN=0`)

---

## References

1. Reddy, J. N. (2006). *An Introduction to the Finite Element Method* (3rd ed.). McGraw-Hill.
2. Reddy, J. N., & Gartling, D. K. (2010). *The Finite Element Method in Heat Transfer and Fluid Dynamics* (3rd ed.). CRC Press.
3. Newmark, N. M. (1959). A Method of Computation for Structural Dynamics. *ASCE Journal of the Engineering Mechanics Division*, 85(3), 67–94.
