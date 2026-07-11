# Chapter 11: Free Energy Calculations -- Thermodynamic Integration

## Core Commands & Syntax

### TI Window Run (Sander multisander)
```bash
mpirun -np 4 sander.MPI -ng 2 -groupfile groupfile
```

Groups file format:
```
-O -i mdin -p prmtop.0 -c eq1.x -o md1.o -r md1.x -inf mdinfo
-O -i mdin -p prmtop.1 -c eq1.x -o md1b.o -r md1b.x -inf mdinfob
```

### TI Window Run (pmemd -- single prmtop)
```bash
pmemd.cuda -O -i ti.in -o ti.out -p system.prmtop \
    -c equil.rst7 -r ti.rst7 -x ti.nc
```

### TI Window Run (pmemd.cuda TI)
```bash
pmemd.cuda -O -i ti_lambda0.5.in -o ti_lambda0.5.out \
    -p system.prmtop -c equil.rst7 -r ti_lambda0.5.rst7 -x ti_lambda0.5.nc
```
Performance: approximately 70% of a regular pmemd.cuda MD simulation.

## Key Namelist / Input File

### Softcore TI input (pmemd/pmemd.cuda)
```
TI production run, lambda=0.5
 &cntrl
  imin=0, ntx=5, irest=1,
  nstlim=500000, dt=0.002,
  ntf=2, ntc=2,
  temp0=300.0, ntt=3, gamma_ln=2.0,
  ntpr=5000, ntwx=5000, ntwr=50000,
  ioutfm=1, ntxo=2,
  cut=8.0, ntb=2, ntp=1,
  ig=-1,

  icfe=1,           ! Enable TI
  ifsc=1,           ! Use softcore potentials
  clambda=0.5,      ! Lambda value for this window
  scalpha=0.5,      ! Softcore alpha parameter
  scbeta=12.0,      ! Softcore beta (electrostatics)
  klambda=1,        ! Linear mixing (k=1, default for softcore)
  logdvdl=1,        ! Print dV/dlambda summary

  timask1=':1-3',   ! Atoms unique to V0
  scmask1=':1-3',   ! Softcore atoms for V0
  timask2='',        ! Atoms unique to V1 (empty if same region)
  scmask2='',        ! Softcore atoms for V1
 /
```

### TI input flags reference
| Variable | Description | Values |
|----------|-------------|--------|
| `icfe` | Enable free energy calculation | 0=off, 1=on |
| `clambda` | Lambda value (0=V0, 1=V1) | 0.0 to 1.0 |
| `klambda` | Exponent for nonlinear mixing | 1-6 (default 1) |
| `ifsc` | Softcore potentials | 0=off, 1=on |
| `scalpha` | Softcore alpha parameter | 0.5 (default) |
| `scbeta` | Softcore beta (electrostatics) | 12.0 (default) |
| `tishake` | SHAKE handling for TI | 0=sync, 1=remove SHAKE |
| `logdvdl` | Log dV/dlambda per step | 0=off, nonzero=on |
| `dynlmb` | Dynamic lambda increment | 0.0 (off) |
| `timask1` | Unique atoms in V0 (ambmask) | `':1-3'` |
| `timask2` | Unique atoms in V1 (ambmask) | `':4-6'` |
| `scmask1` | Softcore atoms in V0 (ambmask) | `':1-3'` |
| `scmask2` | Softcore atoms in V1 (ambmask) | `':4-6'` |
| `crgmask` | Zero partial charges (ambmask) | `':225'` |
| `ntlambda` | Number of lambda windows | 12-20 |
| `clambda` | Current lambda value | 0.0 to 1.0 |

### Gaussian quadrature weights (for TI integration)
| n | lambda | 1-lambda | weight |
|---|--------|----------|--------|
| 1 | 0.5 | - | 1.0 |
| 2 | 0.21132 | 0.78867 | 0.5 |
| 3 | 0.11270 | 0.88729 | 0.27777 |
|   | 0.50000 | - | 0.44444 |
| 5 | 0.04691 | 0.95308 | 0.11846 |
|   | 0.23076 | 0.76923 | 0.23931 |
|   | 0.50000 | - | 0.28444 |
| 7 | 0.02544 | 0.97455 | 0.06474 |
|   | 0.12923 | 0.87076 | 0.13985 |
|   | 0.29707 | 0.70292 | 0.19091 |
|   | 0.50000 | - | 0.20897 |

## Common Workflows

### Workflow 1: Standard TI for ligand binding free energy
```bash
# Prepare prmtop with both endpoints using tiMerge in parmed
# Then run windows:
for lambda in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do
  cat > ti_${lambda}.in <<EOF
TI lambda=${lambda}
 &cntrl
  icfe=1, ifsc=1, clambda=${lambda},
  scalpha=0.5, scbeta=12.0,
  timask1=':LIG', scmask1=':LIG',
  timask2='', scmask2='',
  nstlim=500000, dt=0.002,
  ntt=3, gamma_ln=2.0, temp0=300.0,
  ntb=2, ntp=1, cut=8.0,
  ntpr=5000, ntwx=5000, ntwr=50000,
  logdvdl=1, ig=-1,
 /
EOF
  pmemd.cuda -O -i ti_${lambda}.in -o ti_${lambda}.out \
    -p system.prmtop -c equil.rst7 -r ti_${lambda}.rst7 -x ti_${lambda}.nc
done
```

### Workflow 2: Two-step TI (charge removal then vdW decoupling)
```bash
# Step 1: Remove charges (no softcore needed for pure charge change)
# timask1=':LIG', scmask1=':LIG', crgmask=':LIG' (zero charges)
# clambda sweeps 0.0 to 1.0

# Step 2: Remove vdW (softcore required)
# scmask1=':LIG', timask1=':LIG' (charges already zero)
# clambda sweeps 0.0 to 1.0
```

### Workflow 3: pKa TI (deprotonation free energy)
```bash
# For deprotonation: protonated -> deprotonated
# The proton being removed has zero vdW in Amber force fields,
# so only charge change matters -- no softcore needed.
# clambda can be 0.0 to 1.0 with klambda=1 (linear mixing)

cat > ti_deprot.in <<EOF
pKa deprotonation TI
 &cntrl
  icfe=1, clambda=0.5,
  timask1=':GLU@H1', scmask1=':GLU@H1',
  timask2='', scmask2='',
  nstlim=500000, dt=0.002,
  logdvdl=1, ig=-1,
 /
EOF
```

### Workflow 4: dV/dlambda analysis
```bash
# Extract dV/dlambda values from output files
# Use alchemical_analysis.py or MBAR
alchemical_analysis.py -d ti_output/ -o results/

# Or use the cpptraj ti command
cpptraj <<EOF
# Read dV/dlambda data
readdata ti_0.0.out ti_0.1.out ... ti_1.0.out
# Integrate using trapezoidal rule
EOF
```

### Workflow 5: ACES (Automated Calculation of Engineered Structures)
ACES automates the setup and execution of TI calculations for relative binding free energies. The workflow automates:
- prmtop preparation with merged endpoints
- Window scheduling
- Multi-GPU execution
- Post-processing with MBAR

## Reference Tables

### Softcore potential equations
The mixed potential for V0 (disappearing):
```
V(disappearing) = 4*epsilon*(1-lambda) * [1/(alpha*lambda + (rij/sigma)^6)^2 - 1/(alpha*lambda + (rij/sigma)^6)]
```
For V1 (appearing), replace `lambda` with `(1-lambda)`.

### TI implementation comparison
| Feature | Sander | PMEMD | PMEMD.cuda |
|---------|--------|-------|------------|
| prmtop files | 2 (same atom count) | 1 (merged) | 1 (merged) |
| Softcore | Yes | Yes | Yes |
| GB support | Yes | Yes | No |
| PME support | Yes | Yes | Yes |
| Performance | ~50% of MD | ~75% of MD | ~70% of GPU MD |
| Multi-GPU | No | No | Only REMD |

### Transformations that can skip softcore
| Transformation | Softcore needed? | Reason |
|----------------|------------------|--------|
| Charge change only | No | No vdW change |
| Proton removal (pKa) | No | H has zero vdW in Amber FF |
| vdW parameters change | Yes | Avoid singularities |
| Atom creation/annihilation | Yes | Avoid divergences at endpoints |
| Single-step (charge+vdW) | Yes | Softcore electrostatics equation |

## Worked Example

### pKa calculation of a titratable residue via TI
For a glutamate deprotonation:
```bash
cat > ti_deprot.in <<EOF
Deprotonation TI, lambda=0.5
 &cntrl
  icfe=1, clambda=0.5, klambda=1,
  timask1=':GLU@H1', scmask1=':GLU@H1',
  timask2='', scmask2='',
  nstlim=1000000, dt=0.001,
  ntc=2, ntf=2, temp0=300.0,
  ntt=3, gamma_ln=2.0,
  ntpr=1000, ntwx=5000, ntwr=50000,
  cut=8.0, ntb=2, ntp=1,
  logdvdl=1, ig=-1,
 /
EOF

# Run multiple windows
for lambda in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do
  sed "s/clambda=0.5/clambda=${lambda}/" ti_deprot.in > ti_${lambda}.in
  pmemd.cuda -O -i ti_${lambda}.in -o ti_${lambda}.out \
    -p system.prmtop -c equil.rst7 -r ti_${lambda}.rst7
done

# Analyze with MBAR
alchemical_analysis.py -d . -p ti_ -q out -o results/
```

## Key Takeaways

1. **Use softcore potentials** (`ifsc=1, scalpha=0.5, scbeta=12.0`) for any transformation involving vdW parameter changes or atom creation/annihilation.
2. **For charge-only changes** (pKa, deprotonation), you can skip softcore because the proton has zero vdW radius in Amber force fields.
3. **Use Gaussian quadrature** for lambda scheduling rather than evenly spaced windows -- it gives better convergence with fewer windows.
4. **Always use logdvdl=1** to print the full dV/dlambda time series for post-processing with MBAR.
5. **pmemd.cuda TI does NOT support GB** -- use sander or pmemd for GB-based TI calculations.

## Connects To

- **Chapter 9 (Production MD)**: TI runs use the same pmemd.cuda infrastructure with additional TI flags
- **Chapter 10 (CPPTRAJ)**: Post-process dV/dlambda data with `ti` analysis command
- **Chapter 12 (MMPBSA)**: MMPBSA is an alternative end-point free energy method; TI is more rigorous but more expensive
- **Chapter 8 (Equilibration)**: Each TI window must be equilibrated before sampling
- **Section 27 (Manual)**: Detailed theory of TI and softcore potentials