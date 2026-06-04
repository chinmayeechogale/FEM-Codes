# Project 03 — 1D Beam Bending Solver (EBT / TBT)

## Overview

This project implements a **1D finite element solver for beam bending** based on two classical beam theories:
- **Euler-Bernoulli Beam Theory (EBT):** Assumes plane sections remain plane and perpendicular to the neutral axis; no shear deformation.
- **Timoshenko Beam Theory (TBT):** Includes transverse shear deformation, suitable for thick beams. Uses a shear correction factor `Ks`.

The solver supports:
- **Linear** static analysis
- **Picard (successive substitution)** iteration
- **Newton-Raphson** iteration
- **Incremental load-stepping** for nonlinear problems
- **Selective/reduced integration** to avoid shear locking (separate Gauss rules for bending and shear terms)
- **Dirichlet, Neumann, and Robin** boundary conditions

---

## Files

| File | Description |
|------|-------------|
| `1D Beam Bending Solver` | Main MATLAB solver script |
| `input-file.yaml` | YAML input file (sample: simply-supported beam, L=50) |
| `GAUSS.g` | Gauss quadrature point/weight database |
| `ReadYaml.m` | YAML reader (YAML-MATLAB v0.4.3 dependency) |

---

## Governing Equations

### Euler-Bernoulli Beam Theory (MODEL = 1)

The EBT beam bending equation (strong form):

```
d²/dx² [ EI * d²w/dx² ] = q(x)
```

DOFs per node: **3** — axial displacement `u`, transverse deflection `w`, rotation `θ = dw/dx`

Element stiffness contributions:
```
K_bending = ∫ EI * (d²N/dx²)(d²N/dx²)ᵀ dx
K_axial   = ∫ EA * (dN/dx)(dN/dx)ᵀ dx
```

### Timoshenko Beam Theory (MODEL = 2)

The TBT accounts for shear deformation with independent rotation `φ_x`:

```
d/dx [ EI * dφ_x/dx ] - Ks*G*A*(dw/dx - φ_x) = 0
d/dx [ Ks*G*A*(dw/dx - φ_x) ] + q(x) = 0
```

DOFs per node: **3** — axial `u`, transverse `w`, rotation `φ_x`

Shear stiffness:
```
G_shear = Ks * E * CS_A / (2*(1 + Nu))
```

**Selective/reduced integration:** Bending terms integrated with `NGP` (full) Gauss points; shear terms integrated with `LGP` (reduced) Gauss points to prevent shear locking.

### Distributed Load

```
q(x) = QX0 + QX1*x + QX2*x²
```

---

## Element & Mesh

| Feature | Description |
|---------|-------------|
| Element type | 1D Lagrangian beam element |
| Element order | `P` (default: P=1, linear, NPE=2) |
| Mesh generation | Built-in `MESH1D`: uniform mesh of `NEM` elements over `[X0, X0+L]` |
| Full integration | `NGP` Gauss points (bending/axial terms; default: 2) |
| Reduced integration | `LGP` Gauss points (shear terms, TBT only; default: 1) |
| DOFs per node | 3 (u, w, φ_x or θ) |

---

## Material & Section Parameters

| Parameter | Symbol | Description |
|-----------|--------|-------------|
| `E` | Young's modulus | Elastic modulus of beam material |
| `CS_A` | A | Cross-sectional area |
| `I` | I | Second moment of area (moment of inertia) |
| `Nu` | ν | Poisson's ratio (used to compute shear modulus) |
| `Ks` | K_s | Shear correction coefficient (TBT only; e.g., 5/6 ≈ 0.8333) |

---

## Nonlinear Solvers

| `NONLIN` | Mode | Description |
|----------|------|-------------|
| `0` | Linear | Single direct solve |
| `1` | Picard (direct iteration) | Iterative substitution; convergence: ‖ΔU‖/‖U‖ < ε |
| `2` | Newton-Raphson | Tangent-based increment; convergence: ‖dU‖ < ε |

**Incremental loading:** `LSTEP` load steps with per-step multipliers in array `DP`. At each step the nonlinear loop iterates up to `ITERMAX` times.

Relaxation parameter `GAMMA` blends previous and current solution during assembly:
```
ELU = (1 - GAMMA)*GPU + GAMMA*GCU
```

---

## Boundary Conditions

| Type | Parameters | Description |
|------|-----------|-------------|
| Dirichlet (essential) | `NSPV`, `ISPV`, `VSPV` | Prescribed nodal displacements / rotations |
| Neumann (natural) | `NSSV`, `ISSV`, `VSSV` | Prescribed nodal forces / moments |
| Robin (mixed/spring) | `NSMB`, `ISMB`, `UREF`, `BETA0`, `BETAU` | Elastic foundation or nonlinear spring support |

---

## Input File (`input-file.yaml`)

```yaml
# DOMAIN DATA
X0: 0      # Left end x-coordinate
L:  50     # Length of the beam

# DISTRIBUTED TRANSVERSE LOAD: q(x) = QX0 + QX1*x + QX2*x^2
Q: [1,0,0]   # QX0, QX1, QX2

# BEAM MATERIAL & SECTION PROPERTIES
E:    30000000    # Young's modulus
CS_A: 1           # Cross-sectional area
Ks:   0.8333      # Shear correction coefficient
Nu:   0.3         # Poisson's ratio
I:    0.08333     # Moment of inertia

# MESH DATA
P:   1    # Element polynomial order (NPE = P+1)
NEM: 8    # Number of elements

# BEAM MODEL
MODEL: 1   # 1 = Euler-Bernoulli (EBT); 2 = Timoshenko (TBT)

# QUADRATURE
NGP: 2    # Full integration Gauss points (bending/axial)
LGP: 1    # Reduced integration Gauss points (shear, TBT only)

# BOUNDARY CONDITIONS
NSPV: 3   # Number of Dirichlet BCs
NSSV: 0   # Number of Neumann BCs
NSMB: 0   # Number of Robin BCs

# DIRICHLET BCs: pin at left end (node 1: w=0), pin at right end (node 9: u=0, phi=0)
ISPV:
  NODES: [1, 9, 9]
  DOFS:  [2, 1, 3]
VSPV: [0, 0, 0]

# ROBIN BC: elastic spring at mid-span (node 5, DOF 2 = transverse)
ISMB:
  NODES: [5]
  DOFS:  [2]
UREF:  [0]
BETA0: [3703.54]
BETAU: [0]

# NONLINEAR ANALYSIS PARAMETERS
NONLIN:  1       # 0 = Linear; 1 = Picard; 2 = Newton-Raphson
ITERMAX: 25      # Max iterations per load step
EPSILON: 1e-3    # Convergence tolerance
GAMMA:   1.0     # Relaxation parameter
LSTEP:   10      # Number of load steps
DP: [1,1,1,1,1,1,1,1,1,1]   # Load multiplier per step

# OUTPUT
RECDOF: 26   # Global DOF printed each load step (mid-span transverse displacement)
```

---

## Usage

1. **Install dependencies:** Add `YAML-MATLAB v0.4.3` to the MATLAB path.
2. **Prepare the input file:** Edit `input-file.yaml` with beam geometry, material properties, BCs, and solver settings.
3. **Set the beam model:** `MODEL: 1` for EBT, `MODEL: 2` for TBT.
4. **Ensure support files:** `GAUSS.g` must be in the working directory.
5. **Run the solver in MATLAB:**
   ```matlab
   run('1D Beam Bending Solver')
   ```
6. **Review output:** The solver prints the value of `RECDOF` (the monitored DOF) at the end of each converged load step.

---

## Key Functions

| Function | Description |
|----------|-------------|
| `ELEMBEAM_EBT` | Euler-Bernoulli element stiffness matrix and force vector |
| `ELEMBEAM_TBT` | Timoshenko element stiffness and force vectors (with selective integration) |
| `MESH1D` | Generates 1D mesh: node coordinates and element connectivity |
| `SHAPE1D` | Evaluates 1D Lagrangian shape functions and derivatives |
| `GAUSS` | Reads Gauss quadrature data from `GAUSS.g` |

---

## Sample Problem

The default input file models a **simply-supported beam** of length 50 units:
- Fixed in transverse direction at both ends (nodes 1 and 9, DOF 2 = 0)
- Fixed axial and rotational DOFs at the right end (node 9, DOFs 1 and 3 = 0)
- Uniform distributed transverse load `q = 1` unit/length
- Elastic spring foundation at the mid-span node (node 5)
- Solved using **Picard iteration** over 10 equal load steps

---

## References

1. Reddy, J. N. (2006). *An Introduction to the Finite Element Method* (3rd ed.). McGraw-Hill.
2. Reddy, J. N. (1997). *On Locking-Free Shear Deformable Beam Finite Elements*. Computer Methods in Applied Mechanics and Engineering.
3. Timoshenko, S. P., & Goodier, J. N. (1951). *Theory of Elasticity*. McGraw-Hill.
