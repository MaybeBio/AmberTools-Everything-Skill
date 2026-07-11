# Chapter 20: Implicit Solvent Models

## Core Commands & Syntax

```bash
# GB Simulation in sander
sander -O -i gb.in -o gb.out -p prmtop -c inpcrd -r gb.rst -x gb.nc

# GB Simulation in pmemd (GPU)
pmemd.cuda -O -i gb.in -o gb.out -p prmtop -c inpcrd -r gb.rst -x gb.nc

# GBNSR6 standalone (single-point solvation energy)
gbnsr6 -i gbnsr6.in -o gbnsr6.out -p prmtop -c inpcrd

# PBSA standalone (PB solvation energy)
pbsa -O -i pbsa.in -o pbsa.out -p prmtop -c inpcrd
pbsa -O -i pbsa.in -o pbsa.out -pqr system.pqr

# PBSA GPU-accelerated
pbsa.cuda -O -i pbsa.in -o pbsa.out -p prmtop -c inpcrd

# 3D-RISM in sander
sander -O -i rism.in -o rism.out -p prmtop -c inpcrd \
  -xvv xvvfile -guv guv -huv huv -cuv cuv

# rism1d (1D-RISM for bulk solvent)
rism1d -i rism1d.in -o rism1d.out

# rism3d.snglpnt (single-point 3D-RISM)
rism3d.snglpnt -i rism3d.in -o rism3d.out -p prmtop -c inpcrd

# Parallel 3D-RISM
rism3d.MPI -i rism3d.in -o rism3d.out -p prmtop -c inpcrd

# rism3d_ng (network-graph 3D-RISM)
rism3d_ng -i rism3d.in -o rism3d.out -p prmtop -c inpcrd
rism3d_ng.MPI -i rism3d.in -o rism3d.out -p prmtop -c inpcrd

# PQR file generation
ambpdb -p prmtop -pqr -c inpcrd > system.pqr

# Electrostatic size estimation
elsize system.pqr
```

## Key Namelists & Input Files

### GB/SA Model: &cntrl Variables

```fortran
&cntrl
  ! GB model selection
  igb=1,              ! GB model:
                      !   0=no GB, 1=HCT (Hawkins/Cramer/Truhlar)
                      !   2=OBC model I (alpha=0.8, beta=0, gamma=2.909125)
                      !   5=OBC model II (alpha=1.0, beta=0.8, gamma=4.85)
                      !   6=vacuum, 7=GBn, 8=GBneck2 (recommended)
                      !   10=Numerical PB (PBSA)

  ! Dielectric constants
  intdiel=1.0,        ! Interior dielectric constant (default: 1.0)
  extdiel=78.5,       ! Exterior solvent dielectric constant (default: 78.5)

  ! Salt
  saltcon=0.0,        ! Monovalent salt concentration in M (default: 0.0)

  ! Surface area
  gbsa=1,             ! 0=no SA, 1=LCPO (CPU), 2=icosahedra (MM-GBSA),
                      ! 3=fast pairwise (GPU), 4=exact geometric dSASA
  surften=0.005,      ! Surface tension kcal/mol/A^2 (default: 0.005)

  ! Born radius control
  rgbmax=25.0,        ! Max distance for Born radius calc (default: 25.0)
  offset=0.09,        ! Dielectric radius offset in A (default: 0.09)

  ! ALPB correction
  alpb=0,             ! 0=off, 1=ALPB correction (recommended with igb>=1)

  ! Statistics
  rbornstat=0,        ! 0=off, 1=print effective Born radii statistics

  ! Non-periodic
  ntb=0,              ! GB requires non-periodic (ntb=0)
  cut=999.0,          ! Cutoff should be > system size for GB
/
```

### GB Model Radii Sets (LEaP)

| igb | Recommended Radii | LEaP Command |
|-----|-------------------|--------------|
| 1 | mbondi | `set default PBradii mbondi` |
| 2, 5 | mbondi2 | `set default PBradii mbondi2` |
| 7 | bondi | `set default PBradii bondi` |
| 8 | mbondi3 | `set default PBradii mbondi3` |

### GB Model Comparison

| igb | Model | Description | Best For |
|-----|-------|-------------|----------|
| 1 | GBHCT | Hawkins-Cramer-Truhlar pairwise | Most tested, nucleic acids |
| 2 | GBOBC I | Onufriev-Bashford-Case model I | Proteins |
| 5 | GBOBC II | OBC model II (recommended) | Proteins, better PB agreement |
| 7 | GBn | Mongan et al. neck correction | Proteins, best PB agreement |
| 8 | GBneck2 | Element-specific GBn modification | Proteins + nucleic acids (ff14SBonlysc) |

### GBneck2 (igb=8) Default Parameters

```fortran
! Scaling parameters
Sh=1.425952, Sc=1.058554, Sn=0.733599,
So=1.061039, Ss=-0.703469, Sp=0.5,
offset=0.195141, gbneckscale=0.826836,

! Protein element-specific alpha/beta/gamma
gbalphaH=0.788440, gbbetaH=0.798699, gbgammaH=0.437334,
gbalphaC=0.733756, gbbetaC=0.506378, gbgammaC=0.205844,
gbalphaN=0.503364, gbbetaN=0.316828, gbgammaN=0.192915,
gbalphaOS=0.867814, gbbetaOS=0.876635, gbgammaOS=0.387882,

! Nucleic acid parameters (end with _nu)
screen_hnu=1.69654, screen_cnu=1.26890,
screen_nnu=1.425974, screen_onu=0.18401,
gb_alpha_hnu=0.53705, gb_alpha_cnu=0.33167,
gb_alpha_nnu=0.68631, gb_alpha_onu=0.60634,
```

### GBION: Implicit Solvent with Explicit Ions

```fortran
&cntrl
  igb=8,              ! Must use GBneck2 with GBION
  gbion=3,            ! 0=off, 1=KGB only, 2=both no ion-type, 3=both (recommended)
  gbsa=3,             ! Non-polar via SASA

  ! GBION parameters for NaCl:
  gi_coef_1_n=0.05, gi_coef_2_pn=0.05,
  intdiel_ion_1_p=54, intdiel_ion_1_n=10,
  intdiel_ion_2_pp=54, intdiel_ion_2_pn=10,
  intdiel_ion_2_nn=10,

  ! GBION parameters for KCl:
  ! gi_coef_1_n=0.05, gi_coef_2_pn=0.05,
  ! intdiel_ion_1_p=36, intdiel_ion_1_n=10,
  ! intdiel_ion_2_pp=36, intdiel_ion_2_pn=10,
  ! intdiel_ion_2_nn=10,
/
```

### PBSA: &cntrl and &pb Namelists

```fortran
&cntrl
  ipb=2,              ! 0=no PB, 1=geometric, 2=level-set (default),
                      ! 4=IIM, 6=IIM+analytical, 7=harmonic-avg, 8=2nd-order
  inp=2,              ! 0=no non-polar, 1=PARSE (SASA only),
                      ! 2=cavity+dispersion (default)
  ntx=1,              ! 1=formatted, 2=unformatted
/
&pb
  ! Physical constants
  epsin=1.0,          ! Solute dielectric constant (default: 1.0)
  epsout=80.0,        ! Solvent dielectric constant (default: 80.0)
  epsmem=1.0,         ! Membrane dielectric (only if membraneopt>0)
  istrng=0,           ! Ionic strength in mM (default: 0, 145=physiological)
  pbtemp=300.0,       ! Temperature for PB equation (K)

  ! Dielectric smoothing
  smoothopt=1,        ! 0=harmonic avg, 1=weighted harmonic (default),
                      ! 2=midpoint rule

  ! Radii and surface
  radiopt=1,          ! 0=prmtop radii, 1=Tan-Luo optimized (default)
  dprob=1.4,          ! Solvent probe radius (default: 1.4 A)
  iprob=2.0,          ! Ion probe radius for Stern layer (default: 2.0 A)

  ! Surface type
  sasopt=0,           ! 0=SES, 1=SAS/VDW, 2=revised density function,
                      ! 3=MLSES (machine-learned SES)

  ! Numerical solvers
  npbopt=0,           ! 0=linear PB (default), 1=nonlinear PB
  solvopt=1,          ! 1=Modified ICCG (default), 2=multigrid, 3=CG,
                      ! 4=SOR, 5=adaptive SOR, 6=damped SOR
  accept=0.001,       ! Convergence criterion (default: 0.001)
  maxitn=100,         ! Max iterations (default: 100)

  ! Grid control
  fillratio=2.0,      ! Ratio grid/solute size (default: 2.0, use 4.0 for ligands)
  space=0.5,          ! Grid spacing in A (default: 0.5)
  nbuffer=0,          ! Grid buffer distance (default: auto)
  nfocus=2,           ! Focusing levels (default: 2, max)
  fscale=8,           ! Coarse/fine grid ratio (default: 8)
  npbgrid=1,          ! Grid regeneration frequency (default: 1)

  ! Surface area
  saopt=0,            ! 0=off, 1=field-view method
  triopt=1,           ! 0=off, 1=trimer arc dots (default)
  arcres=0.25,        ! Arc resolution in A (default: 0.25)

  ! Torch PBSA / machine learning
  mlsestype=1,        ! 1=GENIUSES (CPU/GPU), 2=Con2SES-2D, 3=Con2SES-3D

  ! Implicit membrane
  membraneopt=0,      ! 0=off, 1=uniform, 2=heterogeneous (PCHIP), 3=heterogeneous (Spline)
  mprob=2.70,         ! Membrane probe radius (default: 2.70 A)
  mthick=40.0,        ! Membrane thickness (default: 40.0 A)
  mctrdz=0.0,         ! Membrane center z position
  poretype=0,         ! 0=off, 1=auto pore detection
/
```

### GBNSR6: &gb Namelist

```fortran
&cntrl
  inp=0,              ! 0=no non-polar, 1=compute non-polar
/
&gb
  epsin=1.0,          ! Solute dielectric
  epsout=78.5,        ! Solvent dielectric
  istrng=0,           ! Ionic strength (mM)
  dprob=1.4,          ! Solvent probe radius
  space=0.5,          ! Grid spacing (A)
  arcres=0.2,         ! Arc resolution (A)
  alpb=1,             ! ALPB correction (recommended)
  chagb=0,            ! 0=off, 1=CHAGB charge hydration asymmetry
  cavity_surften=0.005, ! Surface tension (only if inp=1)
  rbornstat=0,        ! Print Born radii
  dgij=0,             ! 0=off, 1=print pairwise electrostatic energies
  radiopt=0,          ! 0=hardcoded CHAGB radii, 1=topology radii
  ROH=0.586,          ! CHAGB RzOH parameter (TIP3P/SPC/E)
  tau=1.47,           ! CHAGB tau parameter
/
```

### 3D-RISM: &cntrl and &rism Namelists

```fortran
&cntrl
  irism=1,            ! 0=off, 1=enable 3D-RISM
  ntwrism=100,        ! Output frequency for RISM files
/
&rism
  ! Closure
  closure='KH',       ! KH/HNC/PSEn (e.g., PSE3)

  ! Convergence
  tolerance=1e-12,    ! Target residual tolerance
  maxstep=10000,      ! Max iterations

  ! Grid
  grdspc=0.5,0.5,0.5, ! Linear grid spacing (A)
  buffer=14.0,        ! Solvent box buffer distance (A)

  ! Periodic boundaries
  periodic='pme',     ! none/pme/ewald

  ! Free energy corrections
  gfCorrection=0,     ! Gaussian fluctuation (0/1)
  pcpluscorrection=0, ! PC+/3D-RISM (0/1)
  uccoeff=0,0,0,0,   ! UC correction coefficients

  ! Long-range asymptotics (open boundary)
  asympcorr=.true.,   ! Use asymptotic corrections
  treeDCF=.true.,     ! Treecode for DCF asymptotics
  treeTCF=.true.,     ! Treecode for TCF asymptotics
  treeCoulomb=.false.,! Treecode for Coulomb

  ! Solvent (from rism1d output)
  xvvfile='solvent.xvv',  ! Required: bulk solvent properties
/
```

### rism1d: &parameters Namelist

```fortran
&PARAMETERS
  THEORY='DRISM',     ! DRISM/XRISM/ARISM
  CLOSURE='KH',       ! KH/HNC/PSE-n
  NR=16384,           ! Grid points
  DR=0.025,           ! Grid spacing (A)
  OUTLIST='x',        ! Output options
  ROUT=384,           ! Max real-space output (A)
  KOUT=0,             ! Max reciprocal-space output (A^-1)
  MDIIS_NVEC=20,      ! MDIIS vectors
  MDIIS_DEL=0.3,      ! MDIIS step size
  TOLERANCE=1e-12,    ! Target tolerance
  TEMPERATURE=298.15, ! Temperature (K)
  DIEPS=78.4,         ! Dielectric constant
  NSP=1,              ! Number of species
/
&SPECIES
  DENSITY=55.5,       ! Density in M
  UNITS='M',          ! M/mM/1/A^3/g/cm^3/kg/m^3
  MODEL='water.mdl',  ! Solvent model file
/
```

## Common Workflows

### GBneck2 MD Simulation

```bash
# 1. Prepare system with GB radii in LEaP
cat > build.in << 'EOF'
source leaprc.protein.ff14SBonlysc
prot = loadPdb protein.pdb
set default PBradii mbondi3
saveAmberParm prot system.prmtop system.inpcrd
quit
EOF
tleap -f build.in

# 2. Minimize in GB
cat > gb_min.in << 'EOF'
&cntrl
  imin=1, maxcyc=5000, ncyc=2500,
  ntb=0, igb=8, cut=999.0,
  ntc=1, ntf=1,
/
EOF
sander -O -i gb_min.in -o gb_min.out -p system.prmtop -c system.inpcrd \
  -r gb_min.rst

# 3. GB MD
cat > gb_md.in << 'EOF'
&cntrl
  imin=0, nstlim=1000000, dt=0.002,
  ntt=3, gamma_ln=2.0, tempi=0.0, temp0=300.0,
  ntb=0, igb=8, cut=999.0,
  ntc=2, ntf=2, ntwx=5000, ntpr=5000,
  saltcon=0.0,
/
EOF
sander -O -i gb_md.in -o gb_md.out -p system.prmtop -c gb_min.rst \
  -r gb_md.rst -x gb_md.nc
```

### PBSA Single-Point Solvation Energy

```bash
cat > pb.in << 'EOF'
&cntrl
  ipb=2, inp=2,
/
&pb
  epsin=1.0, epsout=80.0,
  istrng=0, radiopt=1,
  fillratio=4.0, space=0.5,
/
EOF

pbsa -O -i pb.in -o pb.out -p system.prmtop -c system.inpcrd
```

### 3D-RISM Solvation Free Energy

```bash
# 1. Compute bulk solvent properties with rism1d
cat > rism1d.in << 'EOF'
&PARAMETERS
  THEORY='DRISM', CLOSURE='KH',
  NR=16384, DR=0.025,
  OUTLIST='x', ROUT=384, KOUT=0,
  MDIIS_NVEC=20, MDIIS_DEL=0.3, TOLERANCE=1.e-12,
  TEMPERATURE=298.15, DIEPS=78.4,
  NSP=1,
/
&SPECIES
  DENSITY=55.5, UNITS='M', MODEL='spce.mdl',
/
EOF
rism1d -i rism1d.in -o rism1d.out
# Produces: xvvfile (solvent properties)

# 2. Run 3D-RISM single-point
cat > rism3d.in << 'EOF'
  closure='KH', tolerance=1e-12,
  buffer=14.0, grdspc=0.5,0.5,0.5,
  asympcorr=.true.,
/
EOF
rism3d.snglpnt -i rism3d.in -o rism3d.out -p system.prmtop \
  -c system.inpcrd -xvv xvvfile

# 3. Run 3D-RISM MD in sander
cat > rism_md.in << 'EOF'
&cntrl
  imin=0, nstlim=10000, dt=0.002,
  ntt=3, gamma_ln=2.0, tempi=0.0, temp0=300.0,
  ntb=0, igb=6, cut=999.0,
  irism=1, ntwrism=100,
  ntc=1, ntf=1,
/
&rism
  closure='KH', tolerance=1e-12,
  buffer=14.0, grdspc=0.5,0.5,0.5,
  asympcorr=.true.,
/
EOF
sander -O -i rism_md.in -o rism_md.out -p system.prmtop -c system.inpcrd \
  -r rism_md.rst -x rism_md.nc -xvv xvvfile \
  -guv guv -huv huv -cuv cuv
```

### GBION: DNA with Implicit Solvent and Explicit Ions

```bash
# Requires GBneck2 (igb=8) + mbondi3 radii
cat > gbion_md.in << 'EOF'
&cntrl
  imin=0, nstlim=5000000, dt=0.002,
  ntt=3, gamma_ln=2.0, tempi=0.0, temp0=300.0,
  ntb=0, igb=8, cut=999.0,
  gbion=3, gbsa=3, saltcon=0.0,
  gi_coef_1_n=0.05, gi_coef_2_pn=0.05,
  intdiel_ion_1_p=54, intdiel_ion_1_n=10,
  intdiel_ion_2_pp=54, intdiel_ion_2_pn=10,
  intdiel_ion_2_nn=10,
  ntc=2, ntf=2, ntwx=5000, ntpr=5000,
/
EOF
sander -O -i gbion_md.in -o gbion_md.out -p dna.prmtop -c dna.inpcrd \
  -r gbion_md.rst -x gbion_md.nc
```

### Torch PBSA with ML-SES

```bash
# Requires LibTorch compiled: -DLIBTORCH=ON
cat > torch_pb.in << 'EOF'
&cntrl
  ipb=2, inp=2,
/
&pb
  epsin=1.0, epsout=80.0,
  sasopt=3,                ! Use MLSES (machine-learned surface)
  mlsestype=1,             ! GENIUSES model (GPU ~28x speedup)
/
EOF
pbsa -O -i torch_pb.in -o torch_pb.out -p system.prmtop -c system.inpcrd
```

## Reference Tables

### 3D-RISM Closure Approximations

| closure | Full Name | Description |
|---------|-----------|-------------|
| KH | Kovalenko-Hirata | Default, robust for biomolecules |
| HNC | Hyper-netted chain | More accurate for polar systems |
| PSE-n | Partial series expansion | PSE-1, PSE-2, etc. |

### 3D-RISM Output Files

| Flag | Output | Description |
|------|--------|-------------|
| `-guv` | GUV(r) | Pair distribution function |
| `-huv` | HUV(r) | Total correlation function |
| `-cuv` | CUV(r) | Direct correlation function |
| `-quv` | Charge density | Solute-solvent 3D charge density |
| `-exchem` | Excess chemical potential | 3D map |
| `-solvene` | Solvation energy | 3D map |
| `-entropy` | Solvation entropy | 3D map |

### PBSA Solver Options (solvopt)

| solvopt | Solver | Notes |
|---------|--------|-------|
| 1 | Modified ICCG | Default, fast convergence |
| 2 | Geometric multigrid | 4-level v-cycle, best with nfocus=1 |
| 3 | Conjugate gradient | Needs large maxitn |
| 4 | SOR | Needs large maxitn |
| 5 | Adaptive SOR | Only with npbopt=1 (nonlinear) |
| 6 | Damped SOR | Only with npbopt=1 (nonlinear) |

### GBNSR6 GB Equations

| Equation | Description |
|----------|-------------|
| Canonical GB | Original Still et al. (Eqs. 4.2, 4.3) |
| ALPB | Analytical Linearized PB correction |
| CHAGB | Charge hydration asymmetry (water-model specific) |

## Worked Example: GBneck2 Protein Relaxation

```bash
# 1. Clean PDB and prepare with LEaP
pdb4amber -i protein.pdb -o protein_clean.pdb --reduce

cat > build.in << 'EOF'
source leaprc.protein.ff14SBonlysc
prot = loadPdb protein_clean.pdb
set default PBradii mbondi3
saveAmberParm prot system.prmtop system.inpcrd
quit
EOF
tleap -f build.in

# 2. Multi-stage GB relaxation
# Stage 1: Strong restraint minimization
cat > min1.in << 'EOF'
&cntrl
  imin=1, maxcyc=1000, ncyc=500,
  ntb=0, igb=8, cut=999.0,
  ntr=1, restraint_wt=10.0, restraintmask='@CA',
/
EOF
sander -O -i min1.in -o min1.out -p system.prmtop -c system.inpcrd \
  -r min1.rst -ref system.inpcrd

# Stage 2: Weaker restraint minimization
cat > min2.in << 'EOF'
&cntrl
  imin=1, maxcyc=5000, ncyc=2500,
  ntb=0, igb=8, cut=999.0,
  ntr=1, restraint_wt=2.0, restraintmask='@CA,N,C',
/
EOF
sander -O -i min2.in -o min2.out -p system.prmtop -c min1.rst \
  -r min2.rst -ref min1.rst

# Stage 3: Heating
cat > heat.in << 'EOF'
&cntrl
  imin=0, nstlim=25000, dt=0.002,
  ntt=3, gamma_ln=2.0, tempi=0.0, temp0=300.0,
  ntb=0, igb=8, cut=999.0,
  ntc=2, ntf=2,
  ntr=1, restraint_wt=1.0, restraintmask='@CA,N,C',
  ntwx=500, ntpr=500,
/
EOF
sander -O -i heat.in -o heat.out -p system.prmtop -c min2.rst \
  -r heat.rst -x heat.nc -ref min2.rst

# Stage 4: Production
cat > prod.in << 'EOF'
&cntrl
  imin=0, nstlim=5000000, dt=0.002,
  ntt=3, gamma_ln=2.0, tempi=0.0, temp0=300.0,
  ntb=0, igb=8, cut=999.0,
  ntc=2, ntf=2,
  ntwx=5000, ntpr=5000, ntwr=50000,
/
EOF
sander -O -i prod.in -o prod.out -p system.prmtop -c heat.rst \
  -r prod.rst -x prod.nc
```

## Key Takeaways

1. **GBneck2 (igb=8) with mbondi3 radii and ff14SBonlysc is the current recommended GB model** for proteins and nucleic acids. Always use `ntb=0` (non-periodic) and `cut >= system size`.
2. **PBSA provides numerical PB solvation energies** with `ipb=2` (level-set) and `inp=2` (cavity+dispersion). Use `fillratio=4.0` for small molecules (ligands) and `radiopt=1` for Tan-Luo optimized radii.
3. **3D-RISM** computes solvation structure and thermodynamics from first principles. Requires rism1d first for bulk solvent properties (xvvfile). Use `closure='KH'` for robust convergence.
4. **GBION (igb=8 + gbion=3)** extends GBneck2 with explicit ions, capturing ion-specific effects (Na+, K+, Cl-, CoHex3+) around DNA and nucleosomes. Requires non-default KGB and Kepsilon parameters.
5. **Torch PBSA** provides GPU-accelerated ML-based SES construction (GENIUSES: ~28x speedup) via LibTorch. Enable with `-DLIBTORCH=ON` during CMake configuration. Use `sasopt=3, mlsestype=1`.

## Connects To

- **Chapter 3: Force Fields** - ff14SBonlysc recommended for GB, ff19SB for explicit solvent
- **Chapter 5: System Setup with LEaP** - Setting PBradii (mbondi, mbondi2, mbondi3, bondi)
- **Chapter 12: Running MD** - sander/pmemd integration for GB simulations
- **Chapter 18: Membrane Systems** - PBSA implicit membrane (`membraneopt`)
- **Chapter 19: QM/MM** - GB for QM/MM (`qmgb`), PBSA (`igb=10`)
- AMBER tutorials: `07_Creating_Stable_Systems_and_Running_MD/03-2_Relaxation of Implicit Solvent System GB.md`, `09_Case_Studies/05-5_Simulation of DNA with Implicit Solvent and Explicit Ions (GBION).md`