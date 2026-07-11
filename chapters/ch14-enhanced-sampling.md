# Chapter 14: Enhanced Sampling with REMD & Path Methods

## Core Commands & Syntax

### Temperature-REMD (T-REMD)
```bash
mpirun -np 8 sander.MPI -ng 4 -groupfile groupfile -rem 1 -remlog rem.log
```

### Hamiltonian-REMD (H-REMD)
```bash
mpirun -np 8 sander.MPI -ng 4 -groupfile groupfile -rem 3 -remlog rem.log
```

### pH-REMD
```bash
mpirun -np 8 sander.MPI -ng 4 -groupfile groupfile -rem 4 -remlog rem.log
```

### E-REMD (Redox Potential REMD)
```bash
mpirun -np 8 sander.MPI -ng 4 -groupfile groupfile -rem 5 -remlog rem.log
```

### Multi-Dimensional REMD
```bash
mpirun -np 4 sander.MPI -ng 4 -groupfile groupfile \
    -remd-file remd.dim -remlog rem.log
```
Use `-remd-file` instead of `-rem`. The `-rem` flag is set internally to -1.

### Generate REMD Inputs Automatically
```bash
genremdinputs.py -i ref.mdin -inputs temperatures.txt \
    -groupfile ref.groupfile -randomseed 10
```

### Targeted MD
```bash
sander -O -i tmd.in -o tmd.out -p prmtop -c inpcrd -ref ref.crd
```

### Multiply-Targeted MD (MTMD)
```bash
sander -O -i mtmd.in -o mtmd.out -p prmtop -c inpcrd -mtmd mtmd.inp
```

### AT-REMD (pmemd.cuda only)
```bash
pmemd.cuda -O -i atremd.in -p prmtop -c inpcrd
```

### cpptraj analysis of REMD trajectories
```bash
cpptraj <<EOF
parm prmtop
trajin mdcrd.000 remdtraj remdtrajvalues 300.0
trajout T300.nc netcdf
EOF
```

## Key Namelists / Input Files

### REMD mdin (common to all replicas, except exchanged variable)
```
REMD replica
&cntrl
  imin=0, nstlim=500,       ! steps BETWEEN exchange attempts
  ntx=5, irest=1, dt=0.002,
  ntc=2, ntf=2, ntt=3, gamma_ln=2.0,
  temp0=300.0,               ! varies per replica for T-REMD
  ig=-1,                      ! unique per replica
  numexchg=5000,             ! total exchange attempts
  ntwx=500, ntpr=500, ntwr=5000,
  ntb=2, ntp=1, cut=8.0,
/
```

### Groupfile example (4-replica T-REMD)
```
-O -i tremd.300.mdin -p prmtop -c rst7 -o mdout -r rst7 -x mdcrd
-O -i tremd.325.mdin -p prmtop -c rst7 -o mdout -r rst7 -x mdcrd
-O -i tremd.350.mdin -p prmtop -c rst7 -o mdout -r rst7 -x mdcrd
-O -i tremd.375.mdin -p prmtop -c rst7 -o mdout -r rst7 -x mdcrd
```
Each line = one replica. Suffixes `.000`, `.001`, etc. are auto-appended to output files.

### genremdinputs.py: temperatures.txt
```
TEMPERATURE
Temperature exchange from 300K to 400K
300.0
325.0
350.0
375.0
400.0
```
(reference mdin uses `temp0=TEMPERATURE` and `ig=RANDOMNUM`)

### Multi-dimensional REMD file (remd.dim)
```
# REPLICA EXCHANGE DIMENSION FILE
# Dimension 1
&remd
  exch_type = 'TEMPERATURE'
  title = 'Temperature REMD'
  ngroups = 2
 &end
#-----------------------------
# Dimension 2
&remd
  exch_type = 'pH'
  title = 'pH REMD'
  ngroups = 2
 &end
```
Required when running multi-dimensional REMD with `-remd-file`.

### &tgt namelist for Multiply-Targeted MD
```
&tgt
  refin = 'target1.rst7'
  mtmdrmsd = 3.0, mtmdrmsd2 = 0.0,
  mtmdforce = 1.0, mtmdforce2 = 10.0,
  mtmdstep1 = 0, mtmdstep2 = 100000,
  mtmdvari = 1, mtmdninc = 1000,
  mtmdmask = '@CA'
/
&tgt
  refin = 'target2.rst7'
  ...
/
```
Input ends when `refin = ''` is found.

### NEB input variables
```
&cntrl
  ineb = 1,       ! activate NEB
  neb_nstep = 10000,  ! NEB optimization steps
  skmin = 10.0, skmax = 10.0,  ! spring constants
  neb_nrg = 2,    ! tangent definition
  neb_spring = 1, ! variable spring constants
/
```

### AT-REMD input (pmemd.cuda)
```
&cntrl
  atremd_active = 1,
  atremd_ntgtemp = 0,
  atremd_nremd = 4,
  ...
/
```

## Common Workflows

### 1. Setting up T-REMD
**Step 1:** Prepare identical prmtop and starting coordinates for all replicas. Equilibrate each at its target temp0.

**Step 2:** Set `nstlim` = steps between exchange attempts. Total steps = `nstlim * numexchg`.

**Step 3:** `numexchg` must be identical across all mdin files, or the program hangs.

**Step 4:** Run:
```bash
mpirun -np <N_procs> sander.MPI -ng <N_replicas> \
    -groupfile groupfile -rem 1 -remlog rem.log
```
N_procs must be a multiple of N_replicas.

### 2. Analyzing REMD trajectories with cpptraj
**Step 1:** Extract the ensemble at a specific temperature:
```cpptraj
trajin mdcrd.000 remdtraj remdtrajvalues 300.0
trajout T300.nc netcdf
```

**Step 2:** For multi-dimensional REMD:
```cpptraj
trajin mdcrd.000 remdtraj remdtrajvalues 300.0,7.0
```

### 3. Targeted MD to pull structure toward a reference
**Step 1:** Create reference coordinate file (`-ref` flag).

**Step 2:** Set `itgtmd=1`, `tgtrmsd` (target RMSD), `tgtmdfrc` (force constant), `tgtfitmask`, `tgtrmsmask`:
```
&cntrl
  itgtmd=1, tgtrmsd=3.0, tgtmdfrc=1.0,
  tgtfitmask='@CA', tgtrmsmask='@CA',
  ntr=0,
/
```

**Step 3:** Use weight change (`/`) to vary `tgtrmsd` from initial RMSD to 0 during the simulation.

### 4. NEB for minimum energy path
**Step 1:** Prepare endpoint structures and a set of interpolated intermediate images.

**Step 2:** Create a groupfile for multisander with one line per image:
```
-O -i neb.in -p prmtop -c image_01.rst7
-O -i neb.in -p prmtop -c image_02.rst7
...
```

**Step 3:** Run:
```bash
mpirun -np 8 sander.MPI -ng 8 -groupfile groupfile
```

### 5. REAF (REST2-like enhanced sampling)
Enable in H-REMD with different `reaf_tau` values per replica:
```
&cntrl
  ifreaf=1, reaf_tau=0.0,    ! replica at tau=0 (unscaled)
  reaf_mask1='@1-8',
  gti_add_re=1,
/
```
```
&cntrl
  ifreaf=1, reaf_tau=0.3,    ! replica at tau=0.3 (608 K effective)
  reaf_mask1='@1-8',
  gti_add_re=1,
/
```

## Reference Tables

### REMD Types (`-rem` flag)
| -rem value | Type | Exchanged Variable |
|------------|------|--------------------|
| 1 | T-REMD | `temp0` (temperature) |
| 3 | H-REMD | Coordinates (Hamiltonian) |
| 4 | pH-REMD | `solvph` (solution pH) |
| 5 | E-REMD | `solve` (Redox Potential) |

### REMD command-line flags
| Flag | Purpose | Default |
|------|---------|---------|
| `-rem <N>` | REMD type (1,3,4,5) | None |
| `-remlog <file>` | Exchange log file | `rem.log` |
| `-remtype <file>` | Replica type info file | `rem.type` |
| `-remrandompartner <N>` | Random exchange partner | 0 (neighbors) |
| `-remd-file <file>` | Multi-dimensional REMD file | None |

### Key REMD variables in &cntrl
| Variable | Description | Notes |
|----------|-------------|-------|
| `numexchg` | Total exchange attempts | Must be identical across all replicas |
| `nstlim` | Steps between exchanges | Total steps = nstlim * numexchg |
| `gremd_acyc` | Acyclic H-REMD (no 1st-last exchange) | H-REMD only, default 0 |

### Targeted MD variables
| Variable | Description | Default |
|----------|-------------|---------|
| `itgtmd` | 0=off, 1=targeted MD, 2=MTMD | 0 |
| `tgtrmsd` | Target RMSD value | 0.0 |
| `tgtmdfrc` | Force constant (can be negative) | 0.0 |
| `tgtfitmask` | Atoms for RMS superposition | none |
| `tgtrmsmask` | Atoms for RMSD calculation | none |

### NEB key variables
| Variable | Description |
|----------|-------------|
| `ineb` | 1=activate NEB |
| `skmin`, `skmax` | Spring constant range |
| `neb_nrg` | Tangent definition (1=simple, 2=improved) |
| `neb_spring` | 0=uniform, 1=variable springs |

### REAF variables
| Variable | Description |
|----------|-------------|
| `ifreaf` | 1=enable REAF |
| `reaf_tau` | REAF parameter (0=unscaled, 1=fully scaled) |
| `reaf_temp` | Effective temperature in K |
| `reaf_mask1` | Mask for REAF region |
| `gti_add_re` | Control scaling of energy terms |

## Worked Example

### 4-Replica T-REMD of Alanine Tetrapeptide with Langevin Middle Thermostat

**Step 1:** Prepare 4 equilibration runs at 300, 325, 350, 375 K.

**Step 2:** Create mdin for each replica (identical except `temp0` and `ig`):
```
REMD with Langevin middle scheme
&cntrl
  imin=0, nstlim=500, ntx=5, irest=1, dt=0.002,
  ntc=2, ntf=2, ntt=3, gamma_ln=2.0,
  temp0=300.0, ig=1001,
  numexchg=5000, ntwx=500, ntpr=500, ntwr=5000,
  ntb=1, ntp=0, cut=8.0,
/
```

**Step 3:** Create groupfile:
```
-O -i tremd.300.mdin -p ala4.prmtop -c ala4.300.rst7
-O -i tremd.325.mdin -p ala4.prmtop -c ala4.325.rst7
-O -i tremd.350.mdin -p ala4.prmtop -c ala4.350.rst7
-O -i tremd.375.mdin -p ala4.prmtop -c ala4.375.rst7
```

**Step 4:** Run on 8 MPI processes:
```bash
mpirun -np 8 sander.MPI -ng 4 -groupfile groupfile \
    -rem 1 -remlog rem.log
```

**Step 5:** Analyze exchange success in rem.log. Extract per-temperature ensemble with cpptraj.

## Key Takeaways

1. **`nstlim` = steps between exchanges** in REMD, not total steps. Total steps = `nstlim * numexchg`.
2. **`numexchg` must be identical** across all replicas; `nstlim` should also be kept the same.
3. **Multi-dimensional REMD** uses `-remd-file` instead of `-rem`, with one `&remd` block per dimension.
4. **genremdinputs.py** automates creation of mdin files, groupfile, and remd-file for any REMD simulation.
5. **cpptraj `remdtraj` command** is essential for extracting per-ensemble trajectories from REMD output.
6. **REAF** provides REST2-like enhanced sampling via H-REMD with different `reaf_tau` values.

## Connects To
- **Chapter 13 (Umbrella Sampling & NFE):** NFE umbrella sampling and ABMD with REMD
- **Chapter 15 (Constant pH & Redox):** pH-REMD, E-REMD, cphstats/cestats analysis
- **Chapter 9 (Production MD):** Standard MD setup, thermostats
- **Chapter 10 (cpptraj Analysis):** REMD trajectory reconstruction