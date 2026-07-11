# Chapter 24: paramfit -- Force Field Parameter Fitting

## Core Commands & Syntax

**Single molecule fit:**
```bash
paramfit -i Job_Control.in -p prmtop -c mdcrd -q QM_data.dat \
         -v MEDIUM --random-seed seed
```

**Multiple molecule fit:**
```bash
paramfit -i Job_Control.in -pf prmtop_list -cf mdcrd_list \
         -v MEDIUM --random-seed seed
```

**Wizard mode (recommended for first use):**
```bash
paramfit
```

### Command-line arguments

| Flag | Description |
|------|-------------|
| `-i Job_Control.in` | Job control file |
| `-p prmtop` | Molecular topology file (single molecule) |
| `-c mdcrd` | Coordinate file with multiple conformations |
| `-q QM_data.dat` | Quantum energies, one per line (single molecule) |
| `-pf prmtop_list` | File listing topologies and K values (multiple molecules) |
| `-cf mdcrd_list` | File listing coordinates, counts, and energy files |
| `-v LEVEL` | Verbosity: `LOW`, `MEDIUM`, or `HIGH` |
| `--random-seed seed` | Integer seed for reproducibility |

## Key Namelists / Input Files

The job control file uses `variable=value` format, one per line. `#` comments out lines.

### General Options

| Variable | Values | Description |
|----------|--------|-------------|
| `RUNTYPE` | `CREATE_INPUT`, `SET_PARAMS`, `FIT` | Run mode |
| `NSTRUCTURES` | integer | Number of structures in coordinate file (single molecule only) |

### CREATE_INPUT Mode (Quantum Input Files)

| Variable | Description | Default |
|----------|-------------|---------|
| `QMFILEFORMAT` | `GAUSSIAN`, `ADF`, `GAMESS` | -- |
| `QMHEADER` | File prepended to all QM input files | -- |
| `QM_SYSTEM_CHARGE` | Integral charge of system | 0 |
| `QM_SYSTEM_MULTIPLICITY` | Spin multiplicity | 1 |
| `QMFILEOUTSTART` | Output filename prefix | `Job.` |
| `QMFILEOUTEND` | Output filename suffix | `.in` |

### SET_PARAMS Mode

| Variable | Description |
|----------|-------------|
| `PARAMETER_FILE_NAME` | File to store parameter selections |

### FIT Mode

| Variable | Values | Description |
|----------|--------|-------------|
| `ALGORITHM` | `GENETIC`, `SIMPLEX`, `BOTH`, `NONE` | Minimization algorithm |
| `FUNC_TO_FIT` | `SUM_SQUARES_AMBER_STANDARD`, `AMBER_FORCES`, `DIHEDRAL_LEAST_SQUARES` | Fitting function |
| `K` | float | Difference between QM and MM energies |
| `PARAMETERS_TO_FIT` | `DEFAULT`, `K_ONLY`, `LOAD` | How parameters are selected |
| `SCEE` | float | 1-4 electrostatic scaling | 1.2 |
| `SCNB` | float | 1-4 vdW scaling | 2.0 |
| `QM_ENERGY_UNITS` | `HARTREE`, `KCALMOL`, `KJMOL` | Energy unit in QM data | `HARTREE` |
| `QM_FORCE_UNITS` | `HARTREE_BOHR`, `KCALMOL_ANGSTROM` | Force unit in QM data | `HARTREE_BOHR` |
| `WRITE_ENERGY` | filename | Write final AMBER vs QM energies |
| `WRITE_FRCMOD` | filename | Save fitted parameters as frcmod file |
| `SCATTERPLOTS` | flag | Generate bond/angle/dihedral distribution data |
| `SORT_MDCRDS` | `YES`, `NO` | Sort structures by energy | `NO` |
| `COORDINATE_FORMAT` | `TRAJECTORY`, `RESTART` | Input coordinate format | `TRAJECTORY` |
| `CHECK_BOUNDS` | `ON`, `WARN` | Bounds checking behavior | `ON` |

### Genetic Algorithm Options

| Variable | Description | Default |
|----------|-------------|---------|
| `OPTIMIZATIONS` | Population size | 50 |
| `SEARCH_SPACE` | Search range as fraction of original (±100% = 1.0) | entire range |
| `MAX_GENERATIONS` | Maximum iterations | 50000 |
| `GENERATIONS_TO_SIMPLEX` | Stagnant generations before simplex refinement | 10 |
| `GENERATIONS_WITHOUT_SIMPLEX` | Recovery generations between simplex runs | 10 |
| `GENERATIONS_TO_CONV` | Stagnant generations for convergence | 50 |
| `MUTATION_RATE` | Chance of random allele mutation | 0.05 |
| `PARENT_PERCENT` | Fraction passing to next generation | 0.25 |

### Simplex Algorithm Options

| Variable | Description | Default |
|----------|-------------|---------|
| `BONDFC_dx` | Bond force constant step size | 5.0 |
| `BONDEQ_dx` | Bond equilibrium length step size | 0.02 |
| `ANGLEFC_dx` | Angle force constant step size | 1.0 |
| `ANGLEEQ_dx` | Angle equilibrium step size | 0.05 |
| `DIHEDRALBH_dx` | Dihedral barrier height step size | 0.2 |
| `DIHEDRALN_dx` | Dihedral periodicity step size | 0.01 |
| `DIHEDRALG_dx` | Dihedral phase step size | 0.05 |
| `K_dx` | K constant step size | 10.0 |
| `CONV_LIMIT` | Convergence limit | 1.0e-15 |

### Bounds Checking

| Variable | Description | Default |
|----------|-------------|---------|
| `BOND_LIMIT` | Tolerance for bond length results (Angstrom) | 0.1 |
| `ANGLE_LIMIT` | Tolerance for angle results (radians) | 0.05*pi |
| `DIHEDRAL_SPAN` | Required sampling points per pi radian | 12 |

## Common Workflows

### Step 1: Generate QM input files
```
RUNTYPE=CREATE_INPUT
NSTRUCTURES=50
QMFILEFORMAT=GAUSSIAN
QMHEADER=Gaussian.header
```
```bash
paramfit -i Job_Control.in -p prmtop -c mdcrd
# Run Gaussian on generated Job.*.in files
# Extract energies:
$AMBERHOME/AmberTools/src/paramfit/scripts/process_gaussian.x \
    output_directory energies.dat
```

### Step 2: Select parameters to fit
```
RUNTYPE=SET_PARAMS
PARAMETER_FILE_NAME=saved_params
```
```bash
paramfit -i Job_Control.in -p prmtop
```

### Step 3: Fit K (QM/MM energy offset)
```
RUNTYPE=FIT
PARAMETERS_TO_FIT=K_ONLY
FITTING_FUNCTION=SIMPLEX
```
```bash
paramfit -i Job_Control.in -p prmtop -c mdcrd -q energies.dat
```

### Step 4: Fit parameters
```
RUNTYPE=FIT
PARAMETERS_TO_FIT=LOAD
PARAMETER_FILE_NAME=saved_params
FITTING_FUNCTION=GENETIC
ALGORITHM=GENETIC
OPTIMIZATIONS=500
GENERATIONS_TO_CONV=10
GENERATIONS_TO_SIMPLEX=2
GENERATIONS_WITHOUT_SIMPLEX=5
WRITE_FRCMOD=fitted_params.frcmod
K=50.0
```
```bash
paramfit -i Job_Control.in -p prmtop -c mdcrd -q energies.dat
```

### Force fitting (experimental)
```
FUNC_TO_FIT=AMBER_FORCES
K=0.0
QM_FORCE_UNITS=HARTREE_BOHR
QMFILEOUTSTART=output/Job.
QMFILEOUTEND=.gjf.out
```

### Multiple molecule fit
```
# prmtop_list:
molecule1.prmtop 50.0
molecule2.prmtop 100.0

# mdcrd_list:
molecule1.mdcrd 200 energy1.dat
molecule2.mdcrd 100 energy2.dat
```
```bash
paramfit -i Job_Control.in -pf prmtop_list -cf mdcrd_list
```

### Evaluate results
```bash
$AMBERHOME/AmberTools/src/paramfit/scripts/plot_energy.x energy.dat
$AMBERHOME/AmberTools/src/paramfit/scripts/scatterplots.sh
```

## Key Takeaways

1. **paramfit** fits AMBER force field parameters to QM data by minimizing `sum(wi*(E_MM - E_QM + K)^2)` over conformations.
2. Always **fit K first** (`PARAMETERS_TO_FIT=K_ONLY`) before fitting any other parameters; K is the intrinsic QM/MM energy offset.
3. Use the **genetic algorithm** for global search; the **simplex** algorithm excels at refining well-defined systems with few parameters.
4. **Input structure quality matters**: structures must sample the full range of parameters being fit (especially dihedral phases).
5. Supports **parallelization** via OpenMP (`-openmp` configure flag); set `OMP_NUM_THREADS` to control.

## Connects To

- Chapter 18: mdgx parameter fitting (alternative fitting tool)
- Chapter 23: sqm (semi-empirical QM for reference energies)
- Chapter 2: LEaP (apply frcmod output files)
- Chapter 31: External library interface (external QM programs)