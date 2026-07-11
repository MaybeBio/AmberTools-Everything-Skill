# Chapter 33: QM/MM -- Detailed Reference

## Core Commands & Syntax

QM/MM calculations in sander partition the system into a QM region and an MM region. The implementation supports semi-empirical NDDO, DFTB, and ab initio methods via external interfaces.

### Basic QM/MM MD with Semi-empirical Methods

```bash
sander -O -i mdin -o mdout -p prmtop -c inpcrd
```

### QM/MM with External QM Program (e.g., Gaussian)

```bash
sander -O -i mdin -o mdout -p prmtop -c inpcrd
```

### QM/MM with QUICK (ab initio, GPU-accelerated)

```bash
# Using QUICK library linked to sander (recommended)
sander -O -i mdin -o mdout -p prmtop -c inpcrd

# Standalone QUICK
quick water.in

# QUICK with CUDA (multi-GPU)
mpirun -np 2 quick.cuda.MPI water.in
```

---

## QM/MM Theory

The effective Hamiltonian is:

```
H_eff * Psi(x_e, x_QM, x_MM) = E_eff(x_QM, x_MM) * Psi(x_e, x_QM, x_MM)
```

```
E_eff = <Psi| H_QM + H_QM/MM |Psi> + E_MM
```

### Electrostatic Embedding (default)

The QM-MM interaction Hamiltonian (no covalent bonds across boundary):

```
H_QM/MM = sum_q sum_m [ Q_m * h_electron(x_e, x_MM) - Q_m * Z_q * h_core(x_QM, x_MM) + A/r_qm^12 - B/r_qm^6 ]
```

The MM point charges polarize the QM electron density.

### Mechanical Embedding (qmmm_int=5)

```
H_QM/MM = sum_q sum_m [ Q_m*Q_q/r_qm + A/r_qm^12 - B/r_qm^6 ]
```

QM-MM interactions treated classically; no polarization of QM density by MM.

---

## Key Namelists / Input Files

### `&cntrl` Namelist: QM/MM Activation

| Variable | Default | Description |
|----------|---------|-------------|
| `ifqnt` | 0 | Enable QM/MM: 0=off, 1=on |

### `&qmmm` Namelist: General QM/MM Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `iqmatoms` | (none) | Comma-separated atom numbers for QM region |
| `qmmask` | (none) | Mask specifying QM atoms (e.g., `':753'`, `':1-2'`) |
| `qmcut` | same as `cut` | QM/MM electrostatic cutoff (Angstrom) |
| `qm_theory` | `'PM3'` | QM Hamiltonian (see table below) |
| `qmcharge` | 0 | Net charge of QM region |
| `spin` | 1 | Spin multiplicity of QM region |
| `qmmm_int` | 1 | QM/MM interaction method (see table below) |
| `qm_ewald` | 1 (periodic) | Long-range electrostatics: 0=cutoff, 1=PME/Ewald, 2=fixed Mulliken charges |
| `qm_pme` | 1 | PME for QM-MM: 0=regular Ewald, 1=PME (default) |
| `qmshake` | 1 | SHAKE on QM H atoms: 0=no, 1=yes |
| `qmgb` | 2 | GB treatment: 2=polarize QM density, 3=debug (gas-phase charges) |
| `qmmm_switch` | 0 | Switching function at cutoff: 0=no, 1=yes |
| `r_switch_hi` | `qmcut` | Upper switching boundary (Angstrom) |
| `r_switch_lo` | `r_switch_hi - 2` | Lower switching boundary (Angstrom) |
| `kmaxqx,y,z` | 8 | Max k-space vectors for QM-MM Ewald |
| `ksqmaxq` | 100 | Max K^2 for spherical cutoff in reciprocal space |
| `printdipole` | 0 | Dipole moment output: 0=no, 1=QM region, 2=QM+MM |
| `writepdb` | 0 | Write QM region PDB: 0=no, 1=yes |
| `verbosity` | 0 | SQM verbosity level |
| `scfconv` | 1e-8 | SCF convergence criterion |
| `tight_p_conv` | 0 | Tight convergence: 0=no, 1=yes |
| `diag_routine` | 0 | Diagonalization routine |
| `density_predict` | 0 | Density prediction method |
| `ndiis_matrices` | (default) | Number of DIIS matrices |
| `ndiis_iters` | (default) | Number of DIIS iterations |
| `print_charges` | 0 | Print Mulliken charges: 0=no, 1=yes |
| `peptide_corr` | 0 | Peptide bond correction: 0=no, 1=yes |
| `qmqmdx` | 0 | QM-QM distance cutoff |
| `qmqm_anal` | 0 | Analytical QM-QM interactions |
| `qm_mm_anal` | 0 | Analytical QM-MM interactions |
| `itrmax` | (default) | Maximum SCF iterations |
| `parameter_file` | (none) | External parameter file |
| `qxd` | `.false.` | Charge-dependent exchange-dispersion corrections |

### QM Theory Options (`qm_theory`)

| Value | Hamiltonian | Description |
|-------|-------------|-------------|
| `'PM3'` | PM3 | Default semi-empirical NDDO |
| `'AM1'` | AM1 | Austin Model 1 |
| `'MNDO'` | MNDO | Modified Neglect of Diatomic Overlap |
| `'PM6'` | PM6 | PM6 semi-empirical |
| `'PM7'` | PM7 | PM7 semi-empirical |
| `'DFTB'` | SCC-DFTB | Self-consistent charge DFTB |
| `'DFTB3'` | DFTB3 | Third-order DFTB |
| `'EXTERN'` | External | External QM program via interface |
| `'QUICK'` | QUICK | QUICK ab initio (HF/DFT) |

### QM/MM Interaction Methods (`qmmm_int`)

| Value | Name | Description |
|-------|------|-------------|
| 0 | No QM-MM | Electrostatic QM-MM interactions turned off |
| 1 | Default | Standard QM-MM (electrostatic embedding); MM charges in 1-electron Hamiltonian |
| 2 | CHARMM-style | As 1 plus Gaussian core-core terms for AM1/PM3 QM-MM interactions |
| 3 | PM3/MM* | Reformulated QM core-MM charge potential (QM limited to H,C,N,O) |
| 4 | PM3/MMX2 | Extended PM3/MM* with rho_mm parameter (H,C,N,O,S) |
| 5 | Mechanical embedding | Classical point-charge QM-MM interactions; no QM polarization |

### `&qmmm` Namelist: Link Atom Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `lnk_dis` | 1.09 | Distance from QM atom to link atom (Angstrom); negative = place on MM link pair atom |
| `lnk_method` | 1 | Valence terms at boundary: 1=include MM terms crossing boundary, 2=exclude link-pair-containing terms |
| `lnk_atomic_no` | 1 | Atomic number of link atom (default: Hydrogen) |
| `adjust_q` | 2 | Charge conservation: 0=none, 1=nearest MM atoms, 2=all MM atoms (default) |

### DFTB-specific Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `dftb_disper` | 0 | Dispersion correction for DFTB |
| `dftb_3rd_order` | 0 | Third-order DFTB |
| `dftb_chg` | 0 | DFTB charge parameter |
| `dftb_telec` | 0 | Electronic temperature for DFTB |
| `dftb_maxiter` | (default) | Max SCF iterations for DFTB |

### External QM Program Namelists

To use external QM programs, set `qm_theory='EXTERN'` and include one of:

- `&adf` -- ADF (Amsterdam Density Functional)
- `&gms` -- GAMESS-US
- `&nw` -- NWChem
- `&gau` -- Gaussian
- `&orc` -- ORCA
- `&qc` -- Q-Chem
- `&tc` -- TeraChem
- `&mrcc` -- MRCC
- `&quick` -- QUICK (also available as linked library)

Common external namelist variable: `use_template` (0=use namelist parameters, 1=use template input file).

---

## Link Atom Placement

Link atoms are placed along the bond vector between the QM and MM atom of the QM-MM covalent pair:

```
d_L-QM = 1.09 Angstrom (default, equals C-H equilibrium distance)
```

Link atom type: Hydrogen by default (`lnk_atomic_no=1`).

### Valence Term Treatment

**`lnk_method=1` (default):**
- Bonds: MM-MM, MM-MML, MML-QM
- Angles: MM-MM-MM, MM-MM-MML, MM-MML-QM, MML-QM-QM
- Dihedrals: MM-MM-MM-MM, MM-MM-MM-MML, MM-MM-MML-QM, MM-MML-QM-QM, MML-QM-QM-QM

**`lnk_method=2`:**
- Bonds: MM-MM, MM-MML
- Angles: MM-MM-MM, MM-MM-MML, MM-MML-QM
- Dihedrals: MM-MM-MM-MM, MM-MM-MM-MML, MM-MM-MML-QM, MM-MML-QM-QM

### Electrostatic Interactions Around Link Atoms

- Link atoms interact with the full MM electrostatic field (within cutoff)
- MM link pair atoms (MM atoms bound to QM atoms) are excluded from link atom electrostatics
- VDW interactions are NOT calculated for link atoms
- VDW interactions between real QM atoms and all MM atoms (including MM link pair) are calculated

### Link Atom Charge Conservation

When `adjust_q=2` (default), the difference between `qmcharge` and the sum of prmtop charges of QM + link-pair atoms is distributed equally over all MM atoms (excluding link pair atoms).

---

## Ewald and PME for QM/MM

### Inside cutoff (r < qmcut)

```
E_QM/MM = -q_m * (mu_a nu_a, s_m s_m) + Z_a * q_m * (s_a s_a, s_m s_m) * (1 + scale)
```

### Outside cutoff (r > qmcut)

```
E_QM/MM = q_m * (Z_a - sum c_mu_mu) / r
```

### Switching Function (qmmm_switch=1)

```
E_QM/MM = E_QM/MM(r<cutoff) * s(r) + E_QM/MM(r>cutoff) * (1 - s(r))
```

Smoothly connects the two potentials to avoid energy drift in NVE simulations.

---

## QUICK: Ab Initio QM

### Supported Features (QUICK 24.03)

- HF and DFT (LDA, GGA, Hybrid-GGA) energy calculations
- Gradient and geometry optimization
- Restricted closed-shell and unrestricted open-shell
- Grimme dispersion corrections
- Mulliken charge analysis
- Basis functions up to d (f on CUDA)
- MPI parallelization for CPU
- Single-GPU CUDA/HIP and multi-GPU MPI+CUDA/MPI+HIP

### QUICK Input Example

```
B3LYP BASIS=cc-pVDZ CHARGE=0 MULT=1 GRADIENT

O      -0.06756756  -0.31531531   0.00000000
H       0.89243244  -0.31531531   0.00000000
H      -0.38802215   0.58962052   0.00000000
```

---

## Generalized Born for QM/MM

When `qmgb=2` (default), the GB polarization terms are included in the Fock matrix. Mulliken charges from the converged QM calculation are used for GB energy. The SCF procedure solves:

```
dE_eff/dc_ij = 0
```

with a Fock matrix modified by both MM charges and GB polarization terms.

---

## SQM Variables (from Chapter 9)

These are also available in the `&qmmm` namelist when using built-in semi-empirical methods:

| Variable | Description |
|----------|-------------|
| `scfconv` | SCF convergence criterion |
| `tight_p_conv` | Tight convergence flag |
| `diag_routine` | Diagonalization routine selection |
| `pseudo_diag` | Pseudo-diagonalization |
| `pseudo_diag_criteria` | Pseudo-diagonalization criteria |
| `printcharges` | Print Mulliken charges |
| `density_predict` | Density prediction method |
| `ndiis_matrices` | Number of DIIS matrices |
| `ndiis_iters` | DIIS iterations |
| `verbosity` | Output verbosity level |

---

## Common Workflows

### Workflow 1: Basic Semi-empirical QM/MM MD

```
&cntrl
    imin=0, nstlim=10000,
    dt=0.002,
    ntt=1, tempi=0.1, temp0=300.0,
    ntb=1,
    ntf=2, ntc=2,
    cut=8.0,
    ifqnt=1
/
&qmmm
    qmmask=':753',
    qmcharge=-2,
    qm_theory='PM3',
    qmcut=8.0
/
```

### Workflow 2: QM/MM with Link Atoms

```
&qmmm
    qmmask=':1-50',
    qmcharge=0,
    qm_theory='DFTB',
    lnk_dis=1.09,
    lnk_method=1,
    lnk_atomic_no=1,
    adjust_q=2,
/
```

### Workflow 3: External QM Program (Gaussian)

```
&cntrl
    ifqnt=1,
/
&qmmm
    qmmask='@*',
    qmcharge=0,
    spin=1,
    qm_theory='EXTERN',
/
&gau
    method='HF',
    basis='6-31G*',
    scf='tight',
/
```

### Workflow 4: QM/MM with PM3/MM* Reformulated Potential

```
&qmmm
    qmmask=':1-100',
    qmcharge=0,
    qm_theory='PM3',
    qmmm_int=3,
/
```

### Workflow 5: Mechanical Embedding

```
&qmmm
    qmmask=':1-100',
    qmcharge=0,
    qm_theory='PM3',
    qmmm_int=5,
/
```

---

## Hints for Successful QM/MM Calculations

1. **Prmtop requirements**: Mass, charges, vdW, and GB radii must be present. Bond/angle/dihedral parameters involving QM atoms are neglected (except SHAKE).
2. **Choosing QM region**: 80-100 atoms max for semi-empirical. Cut non-polar bonds (C-C). Avoid cutting unsaturated or polar bonds.
3. **Link atom restrictions**: One link atom per MM link pair atom. Cannot cut C-H bonds.
4. **Electrostatic cutoff**: For non-periodic GB, 15-20 Angstrom. For periodic, 8-9 Angstrom.
5. **Parallel scaling**: Good to ~8 CPU cores for semi-empirical. Matrix diagonalization dominates for large QM regions.
6. **SHAKE caution**: QM region SHAKE constrains to MM equilibrium bond lengths. Use `noshmask` to remove SHAKE from perturbed regions.

---

## Worked Example: QM/MM MD with PM3 and Link Atoms

Simulate an enzyme active site with residue 753 treated quantum mechanically:

```
&cntrl
    imin=0, nstlim=10000,
    dt=0.002,
    ntt=1, tempi=0.1, temp0=300.0,
    ntb=1,
    ntf=2, ntc=2,
    cut=8.0,
    ifqnt=1
/
&qmmm
    qmmask=':753',
    qmcharge=-2,
    qm_theory='PM3',
    qmcut=8.0,
    lnk_dis=1.09,
    lnk_method=1,
    adjust_q=2,
/
```

This treats residue 753 with PM3, adds hydrogen link atoms where QM-MM bonds are cut, and distributes excess charge.

---

## Key Takeaways

1. **Electrostatic embedding is default**: QM electron density is polarized by MM charges. Use `qmmm_int=5` for mechanical embedding.
2. **Link atoms are automatic**: Hydrogen link atoms are placed at 1.09 Angstrom from QM atoms along QM-MM bond vectors. No manual specification needed.
3. **External QM programs require template or namelist**: Set `qm_theory='EXTERN'` and include the program-specific namelist (`&gau`, `&orc`, `&tc`, etc.).
4. **QUICK provides GPU-accelerated ab initio**: HF and DFT on GPUs via CUDA/HIP. Use the linked library version for best performance.
5. **PME is not supported for external QM**: Long-range QM-MM electrostatics are truncated at `qmcut` for external programs. Use non-periodic simulations with large cutoff.

## Connects To

- Chapter 9: SQM (semi-empirical quantum mechanics)
- Chapter 10: QUICK (ab initio quantum chemistry)
- Chapter 23: sander (MD engine)
- Chapter 4: Generalized Born
- Chapter 34: Free energy calculations (QM/MM TI)