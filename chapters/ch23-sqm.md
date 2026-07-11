# Chapter 23: sqm -- Semi-empirical Quantum Chemistry

## Core Commands & Syntax

```bash
sqm [-O] -i <input-file> -o <output-file>
```

Default input file name is `mdin`, default output is `mdout`. The `-O` flag overwrites existing output files.

The input file format uses a `&qmmm` namelist followed by atom coordinates:

```
Run semi-empirical minimization
 &qmmm
   qm_theory='AM1',    qmcharge=0,
 /
   6    CG        -1.9590         0.1020                       0.7950
   6   CD1        -1.2490         0.6020                      -0.3030
   6   CD2        -2.0710         0.8650                       1.9630
   ...
```

One line per atom: atomic number, atom name, and Cartesian coordinates (free format).

For external point charges, add `#EXCHARGES` block between `#EXCHARGES` and `#END` markers:

```
#EXCHARGES
6 C -4.7106131    0.0413373    2.1738637   -0.03140
...
#END
```

This enables simplified QM/MM calculations (electrostatic embedding only, no vdW) for single-point energy evaluations. Use `qmmm_int` to set coupling method.

## Key Namelists / Input Files

### `&qmmm` Namelist Variables

**Hamiltonian Selection:**
| Variable | Description | Default |
|----------|-------------|---------|
| `qm_theory` | Hamiltonian: `AM1`, `RM1`, `MNDO`, `PM3`, `PDDG/PM3`, `PDDG/MNDO`, `PM3-CARB1`, `PM3-MAIS`, `MNDO/d` (or `MNDOD`), `AM1/d` (or `AM1D`), `PM6`, `DFTB2` (or `DFTB`), `DFTB3` | `PM3` |
| `qmcharge` | Charge on QM system (integer) | 0 |
| `spin` | Multiplicity; currently only singlet=1 is supported | 1 |

**DFTB-Specific:**
| Variable | Description | Default |
|----------|-------------|---------|
| `dftb_slko_path` | Path to Slater-Koster parameter files | `$AMBERHOME/dat/slko/mio-1-1/` (DFTB2), `$AMBERHOME/dat/slko/3ob-3-1/` (DFTB3) |
| `dftb_disper` | Dispersion correction for DFTB2 (mio-1-1); 0=off, 1=on | 0 |
| `dftb_3rd_order` | Third-order diagonal corrections: `''` (none), `'PA'`, `'PR'`, `'READ'`, or filename | `''` |
| `dftb_chg` | Charge type: 0=Mulliken, 2=CM3 (DFTB2+mio-1-1 only) | 0 |
| `dftb_telec` | Electronic temperature (K) for SCC convergence | 0.0 |
| `dftb_maxiter` | Max SCC iterations before resetting Broyden | 70 |

**SCF Convergence:**
| Variable | Description | Default |
|----------|-------------|---------|
| `scfconv` | SCF convergence criterion (kcal/mol) | 1.0e-8 |
| `tight_p_conv` | 0=loose density convergence, 1=tight (same as scfconv) | 0 |
| `itrmax` | Maximum SCF iterations | 1000 |
| `errconv` | Max absolute error matrix element (hartree); NDDO only | 1.0e-1 |
| `vshift` | Level shift (eV) for virtual orbitals; NDDO only | 0.0 |

**Diagonalization Control (NDDO only):**
| Variable | Description | Default |
|----------|-------------|---------|
| `diag_routine` | 0=auto-select, 1=internal, 2=dspev, 3=dspevd, 4=dspevx, 5=dsyev, 6=dsyevd, 7=dsyevr | 0 |
| `pseudo_diag` | 0=full diagonalization, 1=pseudo diagonalizations when possible | 1 |
| `pseudo_diag_criteria` | Density matrix difference threshold for pseudo diag | 0.05 |

**DIIS Extrapolation (NDDO only):**
| Variable | Description | Default |
|----------|-------------|---------|
| `ndiis_attempts` | Number of DIIS extrapolation attempts (after iteration 100) | 0 |
| `ndiis_matrices` | Number of matrices in DIIS extrapolation | 6 |

**Output Control:**
| Variable | Description | Default |
|----------|-------------|---------|
| `verbosity` | 0=minimal, 1=SCF per step, 2=+memory/pairs, 3=+SCF convergence, 4=+forces, 5=+kJ/mol | 0 |
| `printcharges` | 0=no charges, 1=print Mulliken charges every ntpr | 0 |
| `print_eigenvalues` | 0=off, 1=end of calc, 2=every SCF cycle, 3=each SCF step | 1 |

**Minimization Control:**
| Variable | Description | Default |
|----------|-------------|---------|
| `maxcyc` | Max minimization cycles (TNCG method); 0=single-point | 9999 |
| `ntpr` | Print frequency for minimization progress | 10 |
| `grms_tol` | Gradient RMS tolerance for termination | 0.02 |

**Additional:**
| Variable | Description | Default |
|----------|-------------|---------|
| `qmqmdx` | 1=analytical QM derivatives (default), 2=numerical | 1 |
| `qxd` | Charge-dependent exchange-dispersion vdW corrections | `.false.` |
| `peptide_corr` | 0=off, 1=apply MM correction to peptide linkages | 0 |
| `parameter_file` | Read user-defined parameters from file (NDDO methods) | `''` |
| `qmmm_int` | QM/MM coupling for external point charges | 0 |

## Available Hamiltonians

### MNDO-type (NDDO):
- **PM3**: H, Be, C, N, O, F, Mg, Al, Si, P, S, Cl, Zn, Ga, Ge, As, Se, Br, Cd, In, Sn, Sb, Te, I, Hg, Tl, Pb, Bi
- **AM1**: H, C, N, O, F, Al, Si, P, S, Cl, Zn, Ge, Br, I, Hg
- **RM1**: H, C, N, O, P, S, F, Cl, Br, I
- **MNDO**: H, Li, Be, B, C, N, O, F, Al, Si, P, S, Cl, Zn, Ge, Br, Cd, Sn, I, Hg, Pb
- **PDDG/PM3**: H, C, N, O, F, Si, P, S, Cl, Br, I
- **PDDG/MNDO**: H, C, N, O, F, Cl, Br, I
- **PM3CARB1**: H, C, O
- **PM3-MAIS**: H, O, Cl
- **MNDO/d**: H, Li, Be, B, C, N, O, F, Na, Mg, Al, Si, P, S, Cl, Zn, Ge, Br, Sn, I, Hg, Pb
- **AM1/d**: H, C, N, O, F, Mg, Al, Si, P, S, Cl, Zn, Ge, Br, I, Hg
- **PM6**: H, He, Li, Be, ..., Bi (most elements through atomic number 83)

### Dispersion/H-bond corrections (AM1, PM6):
- `AM1-D*` / `PM6-D`: Dispersion correction only
- `AM1-DH+` / `PM6-DH+`: Dispersion + hydrogen bond correction

### DFTB-type:
- **DFTB2** (SCC-DFTB): self-consistent charge, 2nd order Taylor expansion
- **DFTB3**: 3rd order expansion
- Elements: any atoms for which `.skf` parameters are available from www.dftb.org
- Default parameters: mio-1-1 (DFTB2), 3ob-3-1 (DFTB3)

## Common Workflows

### Single-point energy calculation
```
&qmmm
  qm_theory='PM3', qmcharge=0, maxcyc=0,
/
[atom coordinates]
```

### Geometry optimization
```
&qmmm
  qm_theory='AM1', qmcharge=0, maxcyc=9999,
  scfconv=1.0d-10, tight_p_conv=1,
/
[atom coordinates]
```

### DFTB3 calculation with custom parameters
```
&qmmm
  qm_theory='DFTB3', qmcharge=-1,
  dftb_slko_path='/path/to/3ob-3-1/',
  dftb_maxiter=100, dftb_telec=100.0,
  scfconv=1.0d-9,
/
[atom coordinates]
```

### Single-point with external point charges
```
&qmmm
  qm_theory='PM3', qmcharge=0, maxcyc=0,
  qmmm_int=1,
/
[QM atom coordinates]
#EXCHARGES
6 C ... charge
#END
```

## Key Takeaways

1. **sqm** is a standalone pure-QM program (no MM region); use `sqm -i input -o output` with `&qmmm` namelist.
2. **Hamiltonian selection** via `qm_theory`: use `PM6-DH+` for dispersion+H-bond corrections, `DFTB3` for tight-binding with 3rd-order terms.
3. **SCF convergence** is controlled by `scfconv` (kcal/mol), `tight_p_conv`, and `itrmax`; tighten `scfconv` to 1.0e-10 or better for dynamics.
4. **DFTB parameters** (.skf files) must be obtained from www.dftb.org; use `dftb_slko_path` to specify location.
5. **External point charges** can be included via `#EXCHARGES...#END` block for simplified electrostatic embedding in single-point calculations.

## Connects To

- Chapter 11: QM/MM calculations with sander (full QM/MM using the same `&qmmm` namelist)
- Chapter 10: QUICK ab initio quantum chemistry
- Chapter 24: paramfit (fits parameters to QM energies)
- Chapter 31: External library interface (external QM programs)