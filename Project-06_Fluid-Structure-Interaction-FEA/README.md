# Project 06 — 2D Plate Bending FEA Solver (FSDT)

## Overview

This project implements a **2D finite element solver for plate bending** based on **First-Order Shear Deformation Theory (FSDT)**, also known as **Reissner-Mindlin plate theory**. Unlike classical Kirchhoff thin-plate theory, FSDT includes transverse shear deformation, making it suitable for moderately thick and thick plates.

The solver supports:
- **Linear** static analysis
- **Picard (successive substitution)** iteration
- **Newton-Raphson** iteration with **tangent stiffness** for fast quadratic convergence
- **Geometrically nonlinear** behavior via **Von Kármán strain-displacement relations** (large transverse deflections)
- **Incremental load-stepping** for nonlinear path-dependent problems
- **Orthotropic or isotropic** material properties
- **Selective/reduced integration** to avoid shear locking (separate Gauss rules for bending and shear terms)
- **Dirichlet, Neumann, and Robin** (Winkler foundation) boundary conditions
- **Stress recovery** at top/bottom surfaces and transverse shear stresses

---

## Files

| File | Description |
|------|-------------|
| `2D_Bending.m` | Main MATLAB solver script |
| `input-file.yaml` | YAML input file (sample: square plate, clamped edges + Winkler foundation) |
| `GAUSS.g` | Gauss quadrature point/weight database |
| `ReadYaml.m` | YAML reader (YAML-MATLAB v0.4.3 dependency) |

---

## Governing Equations (FSDT)

### Kinematic Assumptions

FSDT assumes:
- **Straight normals remain straight** but **not necessarily perpendicular** to the mid-plane after deformation.
- Independent rotation variables `φ_x`, `φ_y` represent the rotations of normals about the x and y axes.

### Displacement Field

```
u(x, y, z) = u₀(x, y) + z * φ_x(x, y)
v(x, y, z) = v₀(x, y) + z * φ_y(x, y)
w(x, y, z) = w₀(x, y)
```

- `u₀`, `v₀` = in-plane displacements at the mid-plane
- `w₀` = transverse deflection
- `φ_x`, `φ_y` = rotations of the normal about x and y axes

**DOFs per node: 5** — `[u, v, w, φ_x, φ_y]`

### Strains

**Membrane strains** (mid-plane, z=0):
```
ε_xx₀ = ∂u₀/∂x + 0.5*(∂w₀/∂x)²  (Von Kármán)
ε_yy₀ = ∂v₀/∂y + 0.5*(∂w₀/∂y)²
γ_xy₀ = ∂u₀/∂y + ∂v₀/∂x + (∂w₀/∂x)*(∂w₀/∂y)
```

**Curvatures**:
```
κ_xx = ∂φ_x/∂x
κ_yy = ∂φ_y/∂y
κ_xy = ∂φ_x/∂y + ∂φ_y/∂x
```

**Transverse shear strains**:
```
γ_xz = φ_x + ∂w₀/∂x
γ_yz = φ_y + ∂w₀/∂y
```

### Constitutive Relations (Orthotropic)

**Membrane forces** (`A_ij` = membrane stiffness):
```
N_xx = A11*ε_xx₀ + A12*ε_yy₀
N_yy = A12*ε_xx₀ + A22*ε_yy₀
N_xy = A66*γ_xy₀
```

**Bending moments** (`D_ij` = flexural stiffness):
```
M_xx = D11*κ_xx + D12*κ_yy
M_yy = D12*κ_xx + D22*κ_yy
M_xy = D66*κ_xy
```

**Transverse shear forces** (`A_44`, `A_55` = shear stiffness with shear correction factor `Ks`):
```
Q_xz = A55*γ_xz
Q_yz = A44*γ_yz
```

For an **isotropic** plate:
```
A11 = A22 = E*h/(1-ν²)
A12 = ν*A11
A66 = E*h/(2*(1+ν))
D11 = D22 = E*h³/(12*(1-ν²))
D12 = ν*D11
D66 = E*h³/(24*(1+ν))
A55 = A44 = Ks * G * h
```

### Selective/Reduced Integration

To prevent **shear locking**, the solver uses:
- **Full integration** (`NGP` Gauss points) for membrane and bending terms
- **Reduced integration** (`LGP` Gauss points) for transverse shear terms

---

## Element & Mesh

| Feature | Description |
|---------|-------------|
| Element family | 2D Lagrangian quadrilateral (plate element) |
| Supported NPE | 4 (Q4), 8 (Q8), 9 (Q9) |
| Default (sample input) | 4-node bilinear Q4 |
| Mesh generation | Built-in `MESH2DR`: structured rectangular grid |
| Shape functions | `INTERPLN2D` — Lagrangian |
| Full integration | `NGP` points per direction (membrane/bending; default: 3) |
| Reduced integration | `LGP` points per direction (shear; default: 2) |
| DOFs per node | 5 (u, v, w, φ_x, φ_y) |

---

## Nonlinear Solvers

| `NONLIN` | Mode | Description |
|----------|------|-------------|
| `0` | Linear | Single direct solve (small deflections) |
| `1` | Picard (direct iteration) | Iterative substitution; convergence: ‖ΔU‖/‖U‖ < ε |
| `2` | Newton-Raphson | Tangent stiffness (`TAN` matrix); fast quadratic convergence: ‖dU‖ < ε |

**Incremental loading:** `LSTEP` load steps with per-step multipliers in array `DP`. The load is ramped up incrementally to capture the nonlinear equilibrium path.

Relaxation parameter `GAMMA` blends previous and current solution:
```
ELU = (1 - GAMMA)*GPU + GAMMA*GCU
```

---

## Material & Section Parameters (Orthotropic)

| Parameter | Symbol | Description |
|-----------|--------|-------------|
| `E1` | Young's modulus (direction 1) | Longitudinal modulus |
| `E2` | Young's modulus (direction 2) | Transverse modulus |
| `G12` | Shear modulus (1-2 plane) | In-plane shear |
| `G13` | Shear modulus (1-3 plane) | Transverse shear |
| `G23` | Shear modulus (2-3 plane) | Transverse shear |
| `Nu12` | Poisson's ratio (ν_12) | Major Poisson ratio |
| `H` | Plate thickness | Constant thickness |
| `Ks` | Shear correction factor | Typically 5/6 ≈ 0.8333 for rectangular cross-sections |

**Calculated:**
```
Nu21 = Nu12 * E2 / E1
A11 = E1*H / (1 - Nu12*Nu21)
A12 = Nu21 * A11
A22 = A11 * E2 / E1
D11 = A11 * H² / 12
D12 = Nu21 * D11
D22 = D11 * E2 / E1
A66 = G12 * H
A44 = Ks * G23 * H
A55 = Ks * G13 * H
D66 = G12 * H³ / 12
```

---

## Boundary Conditions

| Type | Parameters | Description |
|------|-----------|-------------|
| Dirichlet (essential) | `NSPV`, `ISPV`, `VSPV` | Prescribed displacements / rotations |
| Neumann (natural) | `NSSV`, `ISSV`, `VSSV` | Prescribed nodal forces / moments |
| Robin (Winkler foundation) | `NSMB`, `ISMB`, `UREF`, `BETA0`, `BETAU` | Elastic foundation: `q = BETA0*(w - UREF) + BETAU*(w - UREF)²` |

---

## Input File (`input-file.yaml`)

```yaml
# DOMAIN DATA
X0: 0
Y0: 0

# DISTRIBUTED TRANSVERSE LOAD
Q: [1,0,0]   # q(x,y) = QX0 + QZX1*x + QZY1*y
FX: [0,0,0]  # In-plane body forces (inactive)
FY: [0,0,0]

# GEOMETRIC & MATERIAL DATA (Orthotropic)
E1:   7.8e6    # Young's modulus (direction 1)
E2:   2.6e6    # Young's modulus (direction 2)
G12:  1.3e6    # In-plane shear modulus
G13:  1.3e6    # Transverse shear modulus
G23:  1.3e6    # Transverse shear modulus
Nu12: 0.25     # Poisson's ratio
H:    1        # Plate thickness
Ks:   0.83333333  # Shear correction factor

# MESH DATA
NPE: 4                                      # Nodes per element (4, 8, 9)
DX: [0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625]  # 8 elements in X
DY: [0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625]  # 8 elements in Y

# QUADRATURE
NGP: 3    # Full integration (membrane/bending)
LGP: 2    # Reduced integration (shear)

# BOUNDARY CONDITIONS
NSPV: 90   # Number of Dirichlet BCs (clamped edges)
NSSV: 0
NSMB: 0

# Dirichlet BCs: clamped on all four edges (u=v=w=φ_x=φ_y=0)
ISPV:
  NODES: [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,10,10,19,19,28,28,37,37,46,46,55,55,64,64,
          9,9,9,9,18,18,18,18,27,27,27,27,36,36,36,36,45,45,45,45,54,54,54,54,63,63,63,63,72,72,72,72,
          73,73,73,74,74,74,75,75,75,76,76,76,77,77,77,78,78,78,79,79,79,80,80,80,81,81,81,81]
  DOFS:  [2,5,2,5,2,5,2,5,2,5,2,5,2,5,2,5,1,4,1,4,1,4,1,4,1,4,1,4,1,4,
          1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,
          1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,4]
VSPV: [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

# Robin BC: Winkler elastic foundation on top edge (nodes 73-81, DOF 3 = w)
ISMB:
  NODES: [73,74,75,76,77,78,79,80,81]
  DOFS:  [3,3,3,3,3,3,3,3,3]
UREF:  [0,0,0,0,0,0,0,0,0]
BETA0: [2.2128e7,2.2128e7,2.2128e7,2.2128e7,2.2128e7,2.2128e7,2.2128e7,2.2128e7,2.218e7]
BETAU: [0,0,0,0,0,0,0,0,0]

# NONLINEAR ANALYSIS
NONLIN:  2       # 0 = Linear; 1 = Picard; 2 = Newton-Raphson
ITERMAX: 100     # Max iterations per load step
EPSILON: 1e-3    # Convergence tolerance
GAMMA:   1.0     # Relaxation parameter
LSTEP:   32      # Number of load steps
DP: [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]

# OUTPUT
RECDOF: 3   # Global DOF to print (transverse displacement at a monitored node)
```

---

## Usage

1. **Install dependencies:** Add `YAML-MATLAB v0.4.3` to the MATLAB path.
2. **Prepare the input file:** Edit `input-file.yaml` with plate geometry, material properties, BCs, and solver settings.
3. **Ensure support files:** `GAUSS.g` must be in the working directory.
4. **Run the solver in MATLAB:**
   ```matlab
   run('2D_Bending')
   ```
5. **Review output:**
   - Monitored DOF (`RECDOF`) printed at each load step
   - Stress tensors computed at top/bottom surfaces and Gauss points via `POSTPRO_FSDT`

---

## Key Functions

| Function | Description |
|----------|-------------|
| `ELEMAT2D_FDST` | Computes element stiffness matrix and force vector for FSDT plate element; includes Von Kármán geometric nonlinearity and tangent matrix (`TAN`) for Newton method |
| `POSTPRO_FSDT` | Post-processes nodal displacements to compute stresses at top/bottom surfaces and transverse shear stresses |
| `MESH2DR` | Generates structured 2D rectangular mesh |
| `INTERPLN2D` | Evaluates Lagrangian shape functions and derivatives |
| `PRECOMPUTE_SF2D` | Pre-computes shape functions at all Gauss points |
| `GAUSS` | Reads Gauss quadrature data |

---

## Sample Problem

The default input file models a **square orthotropic plate** (5 units × 5 units):
- Clamped boundary conditions on all four edges (u=v=w=φ_x=φ_y=0)
- Winkler elastic foundation on the top edge (nodes 73–81)
- Uniform distributed transverse load `q = 1` unit/area
- 8×8 mesh of Q4 elements
- **Newton-Raphson** nonlinear solver over 32 load steps with `DP=5` per step
- Material: E1=7.8e6, E2=2.6e6, G12=1.3e6, Nu12=0.25, H=1, Ks=0.8333

---

## References

1. Reddy, J. N. (2006). *Theory and Analysis of Elastic Plates and Shells* (2nd ed.). CRC Press.
2. Reddy, J. N. (2006). *An Introduction to the Finite Element Method* (3rd ed.). McGraw-Hill.
3. Zienkiewicz, O. C., & Taylor, R. L. (2005). *The Finite Element Method for Solid and Structural Mechanics* (6th ed.). Elsevier.
4. Bathe, K. J. (1996). *Finite Element Procedures*. Prentice Hall.
