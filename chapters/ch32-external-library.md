# Chapter 31: External Library Interface -- QM Programs

## Core Commands & Syntax

The external library interface enables sander to call external QM programs (ADF, GAMESS, Gaussian, NWChem, ORCA, Q-Chem, MRCC, TeraChem) for ab initio and DFT calculations in QM/MM simulations.

```bash
sander -i mdin -o mdout -p prmtop -c inpcrd
sander.MPI -i mdin -o mdout -p prmtop -c inpcrd
```

### Enabling the External Interface

In the `&cntrl` namelist:
```
ifqnt = 1,    ! Enable QM/MM
```

In the `&qmmm` namelist:
```
qm_theory = 'EXTERN',    ! Use external QM program
qmmask = '@*',           ! QM atom selection (here: all atoms)
qmcharge = 0,            ! QM region charge
spin = 1,                ! Spin multiplicity
```

For QM/MM with electronic embedding (default), include all MM point charges as external field:
```
qmcut = <value larger than system size>
```

## Key Namelists / Input Files

### External Program Namelist Selection

| Program | Namelist | Required |
|---------|----------|----------|
| ADF | `&adf` | ADF in PATH |
| GAMESS | `&gms` | GAMESS compiled; `$GMS_PATH` set |
| NWChem | `&nw` | NWChem installed |
| Gaussian | `&gau` | g16/g09/g03 in PATH |
| ORCA | `&orc` | ORCA in PATH |
| Q-Chem | `&qc` | Q-Chem installed |
| MRCC | `&mrcc` | MRCC installed |
| TeraChem | `&tc` | TeraChem installed |

### Common Namelist Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `use_template` | Use user-provided template file (0=no, 1=yes) | 0 |
| `ntpr` | Dipole moment printing frequency (defaults to &cntrl ntpr) | -- |
| `dipole` | Toggle dipole moment writing (0=off, 1=on) | 0 |

### ADF Namelist (&adf)

| Variable | Description | Default |
|----------|-------------|---------|
| `basis` | Basis set: SZ, DZ, DZP, TZP, TZ2P, TZ2P+, ZORA/QZ4P | DZP |
| `core` | Frozen core: None, Small, Medium, Large | None |
| `zlmfit` | Density fit quality for ZLM method | good |
| `fit_type` | Old pair fit basis set type | "" |
| `xc` | XC functional: LDA VWN, GGA BLYP, GGA PBE, HYBRID B3LYP, HYBRID PBE0 | GGA BLYP |
| `scf_iter` | Max SCF cycles | 50 |
| `scf_conv` | SCF convergence threshold | 1.0e-6 |
| `beckegrid` | Becke integration grid: Normal, Good, VeryGood | Good |
| `integration` | Numerical integration accuracy (-1=use Becke) | -1.0 |
| `num_threads` | CPU threads for ADF (0=all cores) | 0 |
| `use_dftb` | Use ADF's DFTB (0=no, 1=yes) | 0 |
| `exactdensity` | Exact electron density for XC (0=no, 1=yes) | 0 |

**ADF template file** (`adf_job.tpl`): Must contain `BASIS ... END` and `SAVE TAPE21`. Do NOT include: UNITS, FRAGMENTS ... END, RESTART, GRADIENT, ATOMS ... END.

### GAMESS Namelist (&gms)

| Variable | Description | Default |
|----------|-------------|---------|
| `basis` | Basis set: STO-3G, 6-31G, 6-31G*, 6-31G**, 6-31+G*, 6-31++G*, 6-311G, 6-311G*, 6-311G**, KTZV, KTZVP, KTZVPP, CCn (n=D,T,Q,5,6), ACCn | 6-31G* |
| `method` | Method: HF, MP2, or DFT functional name | BP86 |
| `nrad` | Radial quadrature points | 96 |
| `nleb` | Angular Lebedev grid points (GAMESS default 302 is too loose) | 590 |
| `scf_conv` | SCF convergence threshold | 1.0e-6 |
| `maxit` | Max SCF iterations | 50 |
| `gms_version` | GAMESS version number | 00 |
| `num_threads` | CPU threads (requires rungms setup) | 1 |
| `mwords` | Memory in millions of 64-bit words | 50 |
| `chelpg` | Calculate CHELPG charges (0=no, 1=yes) | 0 |

**GAMESS template file** (`gms_job.tpl`): `$CONTRL` must contain `RUNTYP=GRADIENT`, `UNIT=ANGS`, `COORD=UNIQUE`. Do NOT include `$DATA`.

### Gaussian Namelist (&gau)

| Variable | Description | Default |
|----------|-------------|---------|
| `basis` | Basis set: STO-3G, 3-21G, 6-31G, 6-311G, with +/++ and */** | 6-31G* |
| `method` | Method: RHF, MP2, BLYP, PBE, B3LYP, etc. | BLYP |
| `scf_conv` | SCF convergence as 10^-N | 8 |
| `num_threads` | CPU threads for Gaussian | 1 |
| `executable` | Gaussian executable name (tries g16, g09, g03) | auto |
| `mem` | Memory allocation string | 256MB |

**Gaussian template file** (`gau_job.tpl`): Route section only, e.g., `#P B3LYP/6-31G* SCF=(Conver=8)`. Do NOT include coordinates, point charges, or Link 0 Commands.

### ORCA Namelist (&orc)

Namelist parameters correspond to ORCA keywords. See ORCA manual for details.

### Template File System

When `use_template = 1`, the program-specific template file overrides namelist settings. Templates provide input not supported via namelists. The format and name vary by program.

## Common Workflows

### QM MD with external Gaussian
```
&cntrl
  imin=0, nstlim=5000, dt=0.001,
  ntpr=50, ntwx=100, ntwr=500,
  ntt=1, temp0=300.0, tautp=1.0,
  ntb=0, cut=999.0,
  ifqnt=1,
/
&qmmm
  qmmask='@*',
  qmcharge=0, spin=1,
  qm_theory='EXTERN',
  qmcut=999.0,
/
&gau
  method='B3LYP',
  basis='6-31G**',
  num_threads=8,
  mem='1GB',
/
```

### QM/MM with ADF (mechanical embedding)
```
&adf
  xc='GGA PBE',
  basis='TZP',
/
```

### QM/MM with GAMESS (16 threads)
```
&gms
  method='DFT',
  dfttyp='PBE',
  basis='6-31G**',
  num_threads=16,
/
```

### Using a template file for Gaussian
```
&gau
  use_template=1,
/
```
Template file `gau_job.tpl`:
```
#P B3LYP/6-31G* SCF=(Conver=8)
```

## Property Output Files

Properties calculated along the trajectory are written to program-specific files:

| Program | File Pattern | Property |
|---------|-------------|----------|
| ADF | `adf_job.dip` | Dipole moment |
| GAMESS | `gms_prop.dip`, `gms_prop.chg` | Dipole moment, CHELPG charges |
| Gaussian | `gau_job.dip` | Dipole moment |

Property files are deleted at the beginning of a run. Back them up if restarting a trajectory.

## MPI Parallel Execution

`sander.MPI` can be used for QM/MM with external QM programs. Each MPI process calls the external QM program. Control external program thread count via `num_threads` in the program-specific namelist. Example: 32 cores, 16 replicas -> each external QM program uses 2 threads.

MPI compatibility: MPICH and MVAPICH work well; OpenMPI does not work.

## Key Takeaways

1. **Enable external QM** via `ifqnt=1` in `&cntrl` and `qm_theory='EXTERN'` in `&qmmm`; add program-specific namelist (`&gau`, `&adf`, `&gms`, `&orc`, etc.).
2. **For QM/MM with electronic embedding**, set `qmcut` larger than the system size to include all MM charges as external field.
3. **Template files** override namelist settings when `use_template=1`; each program has specific template requirements.
4. **Properties** (dipole, charges) are written to separate files; they are deleted on restart, so back them up.
5. **MPI parallel**: each MPI process calls the external QM program independently; use `num_threads` to control per-process threading.

## Connects To

- Chapter 23: sqm (built-in semi-empirical QM)
- Chapter 10: QUICK (built-in ab initio QM)
- Chapter 24: paramfit (fits parameters to QM data)
- Chapter 25: ProPrep (generates QM input files for Gaussian)