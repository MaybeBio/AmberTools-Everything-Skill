# Chapter 13: Umbrella Sampling & Free Energy Methods

## Core Commands & Syntax

### NFE Umbrella Sampling (pmemd/sander)
```bash
pmemd.cuda -O -i mdin -o mdout -p prmtop -c inpcrd -r rst7
```
Set `infe=1` in `&cntrl` to enable the NFE (Nonequilibrium Free Energy) infrastructure.

### Steered Molecular Dynamics (SMD)
```bash
sander -O -i smd.in -o smd.out -p prmtop -c inpcrd
```
Activate via `&smd` namelist in mdin file.

### ABMD (Adaptively Biased MD)
```bash
pmemd.cuda -O -i abmd.in -o abmd.out -p prmtop -c inpcrd
```
Activate via `&abmd` or `&bbmd` namelist. Requires `infe=1`.

### NFE Post-Processing
```bash
nfe-umbrella-slice --help
nfe-umbrella-slice umbrella.nc
```

### WHAM Analysis (external)
```bash
wham <min> <max> <nbins> <tol> <temp> <pad> \
    metadata_file free_energy_file
```

## Key Namelists / Input Files

### &pmd namelist (Umbrella Sampling)
```
title line
&cntrl
  ..., infe = 1
/
&pmd
  output_file = 'pmd.txt'
  output_freq = 50
  cv_file = 'cv.in'
/
```

### cv.in for umbrella sampling (restrain phi/psi of dialanine)
```
&colvar ! phi
  cv_type = 'TORSION'
  cv_ni = 4, cv_i = 5, 7, 9, 15
  anchor_position = -2.05,-2.0,2.0,2.05
  anchor_strength = 500.0,500.0
/
&colvar ! psi
  cv_type = 'TORSION'
  cv_ni = 4, cv_i = 7, 9, 15, 17
  anchor_position = -1.85,-1.8,1.8,1.85
  anchor_strength = 500.0,500.0
/
```

### &smd namelist (Steered MD)
```
title line
&cntrl
  ..., infe = 1
/
&smd
  output_file = 'smd.txt'
  output_freq = 50
  cv_file = 'cv.in'
/
```

### cv.in for SMD (steering a distance from 5.0 to 3.0)
```
&colvar
  cv_type = 'DISTANCE'
  cv_ni = 2
  cv_i = 5, 9
  npath = 2, path = 5.0, 3.0, path_mode = 'LINES',
  nharm = 1, harm = 10.0
/
```

### &abmd namelist (Adaptively Biased MD)
```
title line
&cntrl
  ..., infe = 1
/
&abmd
  mode = 'FLOODING'
  monitor_file = 'abmd.txt'
  monitor_freq = 33
  cv_file = 'cv.in'
  umbrella_file = 'umbrella.nc'
  timescale = 100.0    ! ps, flooding timescale tau_F
  selection_freq = 10000
  selection_constant = 0.001
  wt_temperature = 10000.0
  wt_umbrella_file = 'wt_umbrella.nc'
/
```

### cv.in for ABMD (distance collective variable)
```
&colvar
  cv_type = 'DISTANCE'
  cv_ni = 2, cv_i = 5, 9
  cv_min = -1.0, cv_max = 10.0
  resolution = 0.5   ! required for FLOODING mode
/
```

### &bbmd namelist (Replica 0)
```
&bbmd
  exchange_freq = 100
  exchange_log_file = 'bbmd.log'
  exchange_log_freq = 25
  mt19937_seed = 123455
  mt19937_file = 'mt19937.nc'
  mode = 'FLOODING'
  ...
/
```

### &stsm namelist (Swarms-of-Trajectories String Method)
```
&stsm
  image = 1
  equilibration = 980
  release = 20
  smoothing = 0.1
  report_centers = 'ALL'
  output_file = 'stsm.001.txt'
  output_freq = 10
  cv_file = 'cv.1'
/
```

### &colvar type options
| `cv_type` | Description | `cv_ni` |
|-----------|-------------|---------|
| `'DISTANCE'` | Distance between 2 atoms | 2 |
| `'ANGLE'` | Angle between 3 atoms | 3 |
| `'TORSION'` | Dihedral angle | 4 |
| `'QUATERNION0-3'` | Orientation quaternion components | N atoms |
| `'SPINANGLE'` | Spin angle | N atoms, `cv_nr` ref coords |

### Anchor position/structure for umbrella potential
The umbrella potential U(R) is defined by `anchor_position` = (r1, r2, r3, r4) and `anchor_strength` = (k1, k2):
- U = k1*(r1-r2)*R        (R <= r1)
- U = 0.5*k1*(R-r2)^2     (r1 < R <= r2)
- U = 0                     (r2 < R <= r3)
- U = 0.5*k2*(R-r3)^2     (r3 < R <= r4)
- U = k2*(r4-r3)*R        (R > r4)

## Common Workflows

### 1. 2D Umbrella Sampling with REMD
**Goal:** Free energy surface over phi/psi dihedrals.

**Step 1:** Prepare cv.in with two TORSION collective variables, each with anchor_position defining a grid cell and anchor_strength of ~500 kcal/mol/rad^2.

**Step 2:** Create a groupfile for multi-sander/REMD:
```
-O -i window_01.mdin -p prmtop -c inpcrd
-O -i window_02.mdin -p prmtop -c inpcrd
...
```
Each mdin has `infe=1`, `&pmd` namelist, and unique cv.in with different anchor_position values.

**Step 3:** Run with H-REMD (`-rem 3`):
```bash
mpirun -np N sander.MPI -ng N -groupfile groupfile -rem 3 -remlog rem.log
```

**Step 4:** Post-process with WHAM or nfe-umbrella-slice to reconstruct free energy surface.

### 2. ABMD Flooding Followed by Umbrella Correction
**Step 1:** Run ABMD with `mode = 'FLOODING'` to build biasing potential:
```bash
pmemd.cuda -O -i abmd_flood.in -p prmtop -c inpcrd -r rst7
```
The biasing potential is saved to `umbrella.nc`.

**Step 2:** Run equilibrium umbrella sampling with `mode = 'UMBRELLA'`:
```
&abmd
  mode = 'UMBRELLA'
  umbrella_file = 'umbrella.nc'
  monitor_file = 'corrected.txt'
  monitor_freq = 50
  cv_file = 'cv.in'
/
```

**Step 3:** Correct the free energy: f(xi) = -U(xi) - k_B*T*ln(p_B(xi))

### 3. Steered MD with Jarzynski
**Step 1:** Prepare cv.in with steering path, e.g., distance from 5.0 to 3.0 Angstroms.

**Step 2:** Run SMD:
```bash
sander -O -i smd.in -o smd.out -p prmtop -c inpcrd
```

**Step 3:** Extract work from smd.txt, apply Jarzynski equality:
[exp(-W/k_B*T)] = exp(-DeltaF/k_B*T)

### 4. NFE with Replica Exchange
NFE umbrella sampling (`&pmd`) and ABMD (`&abmd`) work with all REMD types:
- `-rem 0`: Multiple walkers (same U, different trajectories)
- `-rem 3`: H-REMD (different biasing potentials per replica)
- Both compatible with SANDER and PMEMD (GPU accelerated)

## Reference Tables

### &abmd mode comparison
| Mode | Description | Key Parameters |
|------|-------------|----------------|
| `ANALYSIS` | No biasing, just dump CV values | `monitor_file`, `monitor_freq` |
| `UMBRELLA` | Static biasing potential from file | `umbrella_file` (must exist) |
| `FLOODING` | Adaptive biasing, time-dependent | `timescale`, `cv_min`, `cv_max`, `resolution` |

### ABMD extensions
| Extension | Namelist | Key Parameters |
|-----------|----------|----------------|
| Well-Tempered | `&abmd` or `&bbmd` | `wt_temperature`, `wt_umbrella_file` |
| Multiple Walker Selection | `&abmd` | `selection_freq`, `selection_constant`, `selection_epsilon` |
| Driven ABMD | `&abmd` + `&smd` | `driven_weight`, `driven_cutoff` |

### &smd steering variables
| Variable | Description | Default |
|----------|-------------|---------|
| `path` | Steering path control points | Must be specified |
| `npath` | Number of elements in path | 0 |
| `path_mode` | `'SPLINE'` or `'LINES'` | `'SPLINE'` |
| `harm` | Harmonic constant(s) | Must be specified |
| `nharm` | Number of harmonics | 0 |
| `harm_mode` | Like path_mode | `'SPLINE'` |

## Worked Example

### 2D PMF of Alanine Dipeptide phi/psi

**Step 1:** Create cv.in with two TORSION CVs:
```
&colvar ! phi: C-N-CA-C
  cv_type = 'TORSION'
  cv_ni = 4, cv_i = 5, 7, 9, 15
  anchor_position = -2.05,-2.0,2.0,2.05
  anchor_strength = 500.0,500.0
/
&colvar ! psi: N-CA-C-N
  cv_type = 'TORSION'
  cv_ni = 4, cv_i = 7, 9, 15, 17
  anchor_position = -1.85,-1.8,1.8,1.85
  anchor_strength = 500.0,500.0
/
```

**Step 2:** Create mdin for each window:
```
alanine dipeptide umbrella window
&cntrl
  imin=0, nstlim=50000, dt=0.001,
  ntc=2, ntf=2, ntt=3, gamma_ln=1.0,
  temp0=300.0, ig=-1, ioutfm=1,
  ntwx=500, ntpr=500, ntwr=5000,
  ntb=0, cut=999.0, ntp=0,
  infe = 1
/
&pmd
  output_file = 'pmd.txt'
  output_freq = 50
  cv_file = 'cv.in'
/
```

**Step 3:** For each window, adjust `anchor_position` values to cover the full (-pi,pi) x (-pi,pi) grid at ~0.3 rad spacing (36x36 = 1296 windows for dense coverage).

**Step 4:** Use WHAM to combine all windows and compute the 2D PMF.

## Key Takeaways

1. **`infe=1` is mandatory** in `&cntrl` for ALL NFE methods (umbrella, SMD, ABMD, STSM).
2. **`&pmd`** does umbrella sampling, **`&smd`** does steered MD, **`&abmd`** does adaptive biasing, **`&stsm`** does string method -- each uses its own namelist.
3. **Umbrella potential** is defined by 4 anchor points (r1,r2,r3,r4) and 2 spring constants (k1,k2), producing a flat-bottomed harmonic restraint.
4. **ABMD FLOODING** builds a time-dependent biasing potential that converges to the negative free energy; follow with UMBRELLA mode for correction.
5. **All NFE methods work with REMD** (`-rem 3` for H-REMD); SMD results can be analyzed via the Jarzynski equality.

## Connects To
- **Chapter 14 (Enhanced Sampling):** REMD, string method via path collective variables
- **Chapter 9 (Production MD):** Standard MD setup, pmemd flags
- **Chapter 15 (Constant pH & Redox):** Free energy of protonation/redox states
- WHAM analysis tools, Jarzynski equality, Metadynamics concepts