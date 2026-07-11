# Chapter 21: Advanced Tools & Free Energy Methods

## Core Commands & Syntax

### FEW (Free Energy Workflow Tool)

```bash
# Basic FEW invocation
perl FEW.pl <procedure> <command-file>

# Procedure: MMPBSA, MMGBSA, LIE, or TI
perl FEW.pl MMPBSA command_file.txt
```

Command file first line must match the procedure (`@WAMM`, `@LIEW`, `@TIW`). Key flags in command files:

```
# MD Setup (common to all procedures)
lig_struct_path /path/to/ligands
rec_structure /path/to/receptor.pdb
output_path /path/to/output
am1_lig_charges 1      # or resp_lig_charges 1
prepare_leap_input 1
setup_MDsimulations 1
traj_setup_method 3     # 1=1-trajectory, 3=3-trajectory
total_MDequil_time 1000  # equilibration time in ps
total_MDprod_time 5000   # production time in ps
```

### BAR/PBSA Post-Processing

```bash
# Stage 1: Strip solvent and ions from TI trajectories
python bar_pbsa.py strip strip_input.yaml

# Stage 2: Prepare sander PBSA input files
python bar_pbsa.py prep prep_input.yaml

# Stage 3: Run sander in parallel (ligand and complex separately)
python bar_pbsa.py run lig_input.yaml -n 8
python bar_pbsa.py run com_input.yaml -n 8

# Stage 4: Calculate decharging energies via BAR
python bar_pbsa.py calc lig_input.yaml
python bar_pbsa.py calc com_input.yaml
```

YAML input example:
```yaml
dest_path: '1C5X'
ligand_res: 'DRG'
istrng: 150
epsin: 1.0
radiscale: 1.0
protscale: 1.0
```

### GIST (Grid Inhomogeneous Solvation Theory)

```bash
cpptraj -p prmtop << EOF
trajin trajectory.nc
gist doorder doeij gridcntr 30.0 25.0 20.0 griddim 40 30 50 \
     gridspacn 0.5 refdens 0.0334 temp 300.0 \
     prefix gist_out
go
EOF
```

Key GIST grid properties output:
| Dataset | Description |
|---|---|
| `gO`, `gH` | O/H normalized density (gz/rho0) |
| `Esw` | Solute-water interaction energy density |
| `Eww` | Water-water interaction energy density |
| `dTStrans` | Translational entropy density |
| `dTSorient` | Orientational entropy density |
| `dTSsix` | Six-dimensional entropy density |
| `dipole` | Mean dipole moment magnitude |
| `neighbor` | Mean neighbor count |
| `order` | Tetrahedral order parameter |

### EMIL (Absolute Free Energy via Einstein Molecule)

```bash
# mdin setup
&cntrl
  ntt=3, gamma_ln=1.0, ntc=1, ntf=1, dt=0.001,
  emil_do_calc=1, ntp=0,
/
&emil_cntrl
  emil_paramfile = "emilParameters.in",
  emil_logfile = "emil.log",
/
```

emilParameters.in:
```
seed 2325
lambda 0.0
solidRes DC,DG,DA,DT
liquidRes WAT
liquidRes NA
liquidRes CL
epsilonWell 1.0
rWell 0.5
epsilonTrap 0.5
rTrap 5.0
swapTriesPerChain 0.1
printEvery 1000
```

### 3D-RISM

```bash
# 1D-RISM: compute bulk solvent susceptibility
rism1d cSPCE_NaCl > cSPCE_NaCl.out

# 3D-RISM with sander
sander -O -i mdin.rism -o output.r3d \
       -p solute.parm7 -c solute.rst7 \
       -xvv solvent.xvv -guv output
```

mdin.rism:
```
&cntrl
  ntx=1, nstlim=0, irism=1,
/
&rism
  periodic='pme', closure='kh', tolerance=1e-6,
  grdspc=0.35,0.35,0.35, centering=0,
  mdiis_del=0.4, mdiis_nvec=20, maxstep=5000,
  solvcut=9.0, verbose=2, npropagate=0,
  volfmt='dx', ntwrism=1,
/
```

### MoFT (Volumetric Analysis)

```bash
# Laplacian analysis of density
metatwist --dx density.dx --species O \
          --convolve 4 --sigma 1.0 --odx laplacian.dx > output.lp

# Place waters at Laplacian centers
metatwist --dx density.dx --ldx laplacian.dx \
          --map blobsper --species O WAT \
          --bulk 55.55 --threshold 0.5 > output.blobs

# Add hydrogens
gwh -p system.parm7 -w wats.pdb < protein.pdb >> system.wat.pdb
```

### edgembar (MBAR/BAR Analysis)

```bash
edgembar [options] input.xml
edgembar-amber2dats.py --odir ./dats $(ls *.mdout)
```

## Key Namelists / Input Files

### GIST cpptraj action

```
gist [doorder] [doeij] [skipE] [skipS]
     [gridcntr <x> <y> <z>] [griddim <nx> <ny> <nz>]
     [gridspacn <spacing>] [refdens <value>] [temp <T>]
     [prefix <prefix>] [rmsfit <fitmask>]
     [pme cut <cutoff>] [neighborcut <ncut>]
     [solute <mask>] [solventmols <MOLS>]
```

### EMIL &emil_cntrl namelist

```
emil_paramfile = "emilParameters.in"
emil_logfile = "emil.log"
emil_model_infile = "wellsIn.dat"
emil_model_outfile = "wellsOut.dat"
```

### 3D-RISM &rism namelist

```
periodic='pme' | 'nonperiodic'
closure='kh' | 'pse2' | 'pse3' | 'pse4'
grdspc=<dx>,<dy>,<dz>
solvcut=<cutoff>
tolerance=<tol>
mdiis_del=<step>
mdiis_nvec=<nvectors>
maxstep=<maxiters>
ntwrism=<freq>
volfmt='dx' | 'xplor'
```

## Common Workflows

### FEW: Automated MM-PBSA for Factor Xa

```bash
# Step 1A: Parameter file preparation
perl FEW.pl MMPBSA cfiles/leap_am1

# Step 1B: MD simulation setup
perl FEW.pl MMPBSA cfiles/md_setup

# Run MD simulations on the cluster
# (submit generated job scripts)

# Step 2: MM-PBSA setup and analysis
perl FEW.pl MMPBSA cfiles/mm_pbsa_setup
perl FEW.pl MMPBSA cfiles/mm_pbsa_analysis
```

### GIST: Water Thermodynamics in Factor Xa Binding Site

```bash
# 1. Run restrained MD production (protein restrained, waters free)
# 2. Strip counterions and cosolutes
cpptraj -p system.parm7 << EOF
trajin prod.nc
strip :Na+,Cl-
trajout stripped.nc
go
EOF

# 3. Run GIST analysis
cpptraj -p stripped.parm7 << EOF
trajin stripped.nc
gist doorder doeij gridcntr 15.0 20.0 10.0 \
     griddim 40 40 40 gridspacn 0.5 \
     refdens 0.0334 temp 300.0 \
     prefix gist_xa
go
EOF

# 4. Analyze with GISTPP
gistpp --input gist_xa-output.dat
```

### 3D-RISM + MoFT: Place Waters in Crystal

```bash
# 1. Run 3D-RISM
sander -O -i mdin.rism -o out.r3d -p 1aho.parm7 \
       -c 1aho.rst7 -xvv cSPCE_kh.xvv -guv 1aho.kh

# 2. Laplacian analysis
metatwist --dx 1aho.kh.O.0.dx --species O \
          --convolve 4 --sigma 1.0 --odx 1aho.kh.O.dx > 1aho.lp

# 3. Place water molecules
metatwist --dx 1aho.kh.O.0.dx --ldx 1aho.kh.O.dx \
          --map blobsper --species O WAT \
          --bulk 55.55 --threshold 0.5 > 1aho.blobs

# 4. Add hydrogens
grep -v TER *blobs-centroid.pdb > wats.pdb
sed 's/END/TER/' 1aho.amber.pdb > 1aho.wat.pdb
gwh -p 1aho.parm7 -w wats.pdb < 1aho.amber.pdb >> 1aho.wat.pdb
```

### EMIL: Absolute Free Energy of Dialanine

```bash
# 1. Equilibrate at constant pressure
# 2. Find average box size, scale to constant volume
# 3. Prepare emilParameters.in files at lambda=0.0, 0.1, ..., 1.0
# 4. Run pmemd at each lambda
pmemd -O -i mdin -o emil.out -p system.parm7 \
      -c equil.rst7 -r emil.rst7

# 5. Extract dHdL from emil.log files
# 6. Integrate dHdL vs lambda, subtract from EMIL ref free energy
```

## Reference Tables

### FEW Procedures

| Procedure | Module | Key Phrase | Requires |
|-----------|--------|------------|----------|
| MM-PBSA | WAMM | `@WAMM` | MD trajectories |
| MM-GBSA | WAMM | `@WAMM` | MD trajectories |
| LIE | LIEW | `@LIEW` | MD trajectories |
| TI | TIW | `@TIW` | TI window input |

### GIST Grid Parameters

| Parameter | Default | Recommendation |
|-----------|---------|---------------|
| gridspacn | 0.5 A | 0.5 A for detail, 0.75 A for convergence |
| griddim | 40,40,40 | Cover entire binding site + 2-3 A |
| refdens | 0.0334 | Use model-specific bulk density |
| neighborcut | 3.5 A | Water O-O first shell |
| Minimum trajectory | -- | 10-20 ns for room-temp binding site |

### EMIL Key Parameters

| Parameter | Units | Typical Value |
|-----------|-------|--------------|
| epsilonWell | kT | 1.0 |
| rWell | A | 0.5 |
| epsilonTrap | kT | 0.5 |
| rTrap | A | 5.0 |
| swapTriesPerChain | -- | 0.1 |
| scalpha (pmemd) | -- | 0.3 |
| scbeta (pmemd) | -- | 16.0 |

### BAR/PBSA Radii Scaling

| Keyword | Applies To | Purpose |
|---------|-----------|---------|
| radiscale | Ligand atoms | Match PBSA/explicit electrostatic FE |
| protscale | Protein atoms | Match PBSA/explicit electrostatic FE |
| epsin | Solute interior | Screens Coulombic interactions |

## Worked Example

### GIST Analysis of Factor Xa Active Site

From the tutorial, a GIST calculation on Factor Xa's binding pocket reveals water thermodynamics that guide ligand design. Key commands:

```bash
# Strip ions, autoimage, RMS-fit to protein backbone
cpptraj -p fxa.parm7 << EOF
trajin prod.nc 1 last 10
strip :Na+,Cl-
autoimage
rms first :1-282@CA
trajout stripped.nc
go
EOF

# GIST on the binding pocket
cpptraj -p stripped.parm7 << EOF
trajin stripped.nc
gist doorder doeij \
     gridcntr 10.0 15.0 5.0 \
     griddim 40 40 40 \
     gridspacn 0.5 \
     refdens 0.0334 temp 300.0 \
     solute :1-282 \
     prefix gist_fxa
go
EOF
```

This produces `gist_fxa-output.dat` with per-voxel thermodynamics and `.dx` grid files for visualization. Voxels with favorable `Esw` (strong solute-water interactions) and high `gO` (stable occupancy) indicate water positions that contribute favorably to binding. Voxels with high `dTSsix` (entropic penalty) suggest displaceable waters -- potential targets for ligand design.

## Key Takeaways

1. **FEW automates the entire free energy pipeline** -- from MD setup through MM-PBSA/GBSA, LIE, or TI analysis -- for multiple ligands against one receptor. Always use the minimalistic command files as starting templates.
2. **GIST maps water thermodynamics onto a 3D grid** from explicit solvent MD. It identifies favorable vs. displaceable water sites. Use at least 10-20 ns of restrained-solute simulation for convergence.
3. **BAR/PBSA post-processes alchemical TI trajectories** to incorporate electronic polarization via continuum dielectric screening. Calibrate radiscale and protscale before varying epsin.
4. **EMIL computes absolute free energies** by integrating from an analytically tractable reference. Use Langevin thermostat (gamma_ln ~1.0), disable SHAKE, and use dt=0.001.
5. **3D-RISM + MoFT** provides an efficient implicit-solvent alternative for placing water and ions in crystal structures or binding sites without running explicit solvent MD.

## Connects To

- [ch11](ch11-free-energy-ti.md) -- Thermodynamic Integration with softcore potentials
- [ch12](ch12-mmpbsa.md) -- MM-PBSA and MMPBSA.py binding free energy calculations
- [ch13](ch13-umbrella-sampling-nfe.md) -- Umbrella sampling and NFE toolkit
- [ch14](ch14-enhanced-sampling.md) -- Enhanced sampling (REMD, steered MD, WESTPA)
- [ch20](ch20-implicit-solvent.md) -- Implicit solvent models (GB, PBSA, RISM)