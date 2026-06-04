# Project 07 — 2D Geometric Non-Linearity FEM Solver

## Overview

This project implements a **2D finite element solver for geometrically nonlinear structural analysis** under large deformations. Unlike classical linear elasticity, this solver accounts for finite strains and rotations, making it suitable for problems where deformations are large enough to significantly alter the structure's geometry.

The solver supports:
- **Linear** static analysis
- **Picard (successive substitution)** iteration for nonlinear problems
- **Newton-Raphson** iteration for fast quadratic convergence
- **Incremental load-stepping** to trace the nonlinear equilibrium path
- **Plane stress** and **plane strain** assumptions
- Post-processing of **Cauchy**, **First Piola-Kirchhoff (PK1)**, and **Second Piola-Kirchhoff (PK2)** stress tensors

---

## Files

| File | Description |
|------|-------------|
| `2D Geometric Non-Linearity` | Main MATLAB script (input, solver, post-processing) |
| `input-file.yaml` | YAML input file controlling all problem parameters |
| `GAUSS.g` | Gauss quadrature point/weight database |
| `ReadYaml.m` | YAML reader (YAML-MATLAB v0.4.3 dependency) |

---

## Governing Equations

The solver is based on the **Updated Lagrangian formulation** for 2D solid mechanics. The weak form of the equilibrium equation in the current configuration is:

```
∫_Ω [ σ : δε + σ : (∇δu · ∇u) ] dΩ = ∫_Ω f · δu dΩ + ∫_Γ t · δu dΓ
```

where:
- **σ** = Cauchy stress tensor
- **ε** = Euler-Almansi strain tensor
- **δu** = virtual displacement
- **f** = body force vector
- **t** = surface traction

### Euler-Almansi Strains

The nonlinear strain-displacement relations used in the element routine (`ELEMAT2D`) are:

```
e_xx = u_x - 0.5*(u_x² + v_x²)
e_yy = v_y - 0.5*(u_y² + v_y²)
γ_xy = u_y + v_x - u_x*u_y - v_x*v_y
```

where `u_x`, `u_y`, `v_x`, `v_y` are spatial displacement gradients.

### Cauchy Stress (via Constitutive Law)

**Plane Stress:**
```
C11 = E/(1-ν²),  C22 = C11,  C12 = ν*C11,  C66 = E/(2*(1+ν))
```

**Plane Strain:**
```
C11 = E*(1-ν)/((1+ν)*(1-2ν)),  C12 = E*ν/((1+ν)*(1-2ν)),  C66 = E/(2*(1+ν))
```

```
σ_xx = C11*e_xx + C12*e_yy
σ_yy = C12*e_xx + C22*e_yy
σ_xy = C66*γ_xy
```

### Tangent Stiffness Matrix

The element stiffness matrix (`ELK`) is decomposed into:
- **Material stiffness** (`K¹`): derived from the constitutive law
- **Geometric stiffness** (`KG`): accounts for current Cauchy stresses

```
KG = ∫ [ σ_xx*(∂N/∂x)(∂N/∂x)ᵀ + σ_xy*((∂N/∂y)(∂N/∂x)ᵀ + (∂N/∂x)(∂N/∂y)ᵀ) + σ_yy*(∂N/∂y)(∂N/∂y)ᵀ ] t dΩ
```

---

## Element & Mesh

| Feature | Description |
|---------|-------------|
| Element family | 2D Lagrangian quadrilateral |
| Supported NPE | 4 (bilinear Q4), 8 (serendipity Q8), 9 (biquadratic Q9) |
| Default (sample input) | 9-node biquadratic (Q9) |
| Mesh generation | Built-in `MESH2DR` function (structured rectangular grid) |
| Shape functions | `INTERPLN2D` — Lagrangian, mapped from reference [-1,1]² |
| Integration | Gauss quadrature; `NGP` points per direction (default: 3×3 = 9 pts) |
| Stress recovery | Separate 2×2 Gauss rule used in `POST_STRESS_2D` |

---

## Nonlinear Solvers

| `NONLIN` | Mode | Update rule |
|----------|------|-------------|
| `0` | Linear | Direct solve — single iteration |
| `1` | Picard (direct iteration) | `GCU ← GCU + SOLU`, convergence: ‖ΔU‖/‖U‖ < ε |
| `2` | Newton-Raphson | `GCU ← GCU + dU` (tangent-based increment), convergence: ‖dU‖/‖U‖ < ε |

**Incremental loading** is handled via `LSTEP` load steps with increments defined in the `DP` array. At each load step, the iterative solver runs until convergence (or `ITERMAX` is reached).

A **relaxation parameter** `GAMMA` blends the current and previous solution during element assembly:
```
ELU = (1 - GAMMA)*GPU + GAMMA*GCU
```

---

## Post-Processing

The `POST_STRESS_2D` function computes, at each Gauss point:

| Quantity | Symbol | Description |
|----------|--------|-------------|
| Euler-Almansi strains | e_xx, e_yy, γ_xy | Strains in current (spatial) config |
| Green-Lagrange strains | E_XX, E_YY, Γ_XY | Strains in reference config |
| Cauchy stress | σ_xx, σ_yy, σ_xy | True stress in current config |
| 1st Piola-Kirchhoff | P = F·S | Nominal stress (force/reference area) |
| 2nd Piola-Kirchhoff | S = C:E_GL | Reference stress — work conjugate to E_GL |
| Jacobian | J = det(F) | Volume ratio (current / reference) |

---

## Boundary Conditions

| Type | Parameter | Description |
|------|-----------|-------------|
| Dirichlet (essential) | `NSPV`, `ISPV`, `VSPV` | Prescribed nodal displacements |
| Neumann (natural) | `NSSV`, `ISSV`, `VSSV` | Prescribed nodal forces / tractions |
| Robin (mixed) | `NSMB`, `ISMB`, `UREF`, `BETA0`, `BETAU` | Nonlinear spring / foundation |
| Edge traction | `TOP_EDGES`, `BOTTOM_EDGES`, `Q0`, `LOAD_DOF` | Distributed pressure applied to element edges |

For **Newton-Raphson** (`NONLIN = 2`), Dirichlet BCs are enforced as zero-increment corrections; for Linear/Picard modes, total values are prescribed directly.

---

## Input File (`input-file.yaml`)

```yaml
# DOMAIN DATA
X0: 0        # X-coordinate of bottom-left corner
Y0: 0        # Y-coordinate of bottom-left corner

# PDE COEFFICIENTS (not active in geometric nonlinear mode)
A1: [0.2,0,0,0.0004,0,0]   # A10, A1X, A1Y, A1U, A1UX, A1UY
A2: [0.2,0,0,0.0004,0,0]   # A20, A2X, A2Y, A2U, A2UX, A2UY
F:  [0,0,0]                 # FX0, FX1, FY1

# GEOMETRIC & MATERIAL DATA
E:     1.2e7   # Young's modulus
NU:    0.3     # Poisson's ratio
THICK: 0.1    # Plate thickness

# MESH DATA
NPE: 9                    # Nodes per element (4, 8, or 9)
DX: [2,2,2,2,2]           # Element lengths along X
DY: [1]                   # Element lengths along Y
NGP: 3                    # Gauss points per direction

# BOUNDARY CONDITIONS
NSPV: 6                   # Number of Dirichlet BCs
NSSV: 0                   # Number of Neumann BCs
NSMB: 0                   # Number of Robin BCs
ISPV:
  NODES: [1,12,18,1,12,18]
  DOFS:  [1,1,1,2,2,2]
VSPV: [0,0,0,0,0,0]

# EDGE TRACTION LOADING
LOADTYPE: EDGE_TRACTION
Q0: 500                   # Traction magnitude (psi per increment)
TOP_EDGES:
  - [18,19,20]
  - [20,21,22]
  - [22,23,24]
  - [24,25,26]
  - [26,27,28]
BOTTOM_EDGES:
  - [1,2,3]
  - [3,4,5]
  - [5,6,7]
  - [7,8,9]
  - [9,10,11]
LOAD_DOF: 2               # DOF direction for traction (2 = Y)

# NONLINEAR ANALYSIS PARAMETERS
NONLIN:  1       # 0 = Linear; 1 = Picard; 2 = Newton
ITERMAX: 25     # Max iterations per load step
EPSILON: 1e-3   # Convergence tolerance
GAMMA:   1.0    # Relaxation parameter
LSTEP:   18     # Number of load steps
DP: [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]  # Load increments per step
PLANE:   1      # 0 = Plane strain; 1 = Plane stress
```

---

## Usage

1. **Install dependencies:** Add `YAML-MATLAB v0.4.3` to the MATLAB path.
2. **Prepare the input file:** Edit `input-file.yaml` (or `SAMPLE_INPUT_FILE2D.yaml`) with your problem geometry, material properties, boundary conditions, and solver settings.
3. **Ensure support files are present:** `GAUSS.g` must be in the working directory.
4. **Run the solver:** Execute the main script in MATLAB:
   ```matlab
   run('2D Geometric Non-Linearity')
   ```
5. **Review output:** The solver prints per-load-step diagnostics:
   - Node displacements (u, v) at a monitor node
   - Current coordinates of Gauss points
   - Cauchy, PK1, and PK2 stress components at element Gauss points

---

## Key Functions

| Function | Description |
|----------|-------------|
| `ELEMAT2D` | Computes element stiffness (K¹ + KG) and internal/external force vectors using Euler-Almansi strains |
| `POST_STRESS_2D` | Post-processes displacements to compute Cauchy, PK1, PK2 stresses and Green-Lagrange strains |
| `MESH2DR` | Generates structured 2D mesh: node coordinates (GLXY) and connectivity (NOD) |
| `INTERPLN2D` | Evaluates Lagrangian shape functions and their derivatives for Q4/Q8/Q9 elements |
| `PRECOMPUTE_SF2D` | Pre-computes shape functions at all Gauss integration points |
| `GAUSS` | Reads Gauss quadrature points and weights from `GAUSS.g` database |

---

## References

1. Reddy, J. N. (2004). *An Introduction to Nonlinear Finite Element Analysis*. Oxford University Press.
2. Reddy, J. N. (2006). *An Introduction to the Finite Element Method* (3rd ed.). McGraw-Hill.
3. Bonet, J., & Wood, R. D. (2008). *Nonlinear Continuum Mechanics for Finite Element Analysis* (2nd ed.). Cambridge University Press.
4. Belytschko, T., Liu, W. K., & Moran, B. (2000). *Nonlinear Finite Elements for Continua and Structures*. Wiley.
