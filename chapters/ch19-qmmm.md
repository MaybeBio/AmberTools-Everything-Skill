# Chapter 19: QM/MM Calculations

## Core Commands & Syntax

```bash
# QM/MM with built-in semiempirical methods (sander)
sander -O -i qmmm.in -o qmmm.out -p prmtop -c inpcrd -r qmmm.rst -x qmmm.nc

# QM/MM with external QM program (Gaussian)
sander -O -i qmmm_ext.in -o qmmm_ext.out -p prmtop -c inpcrd -r ext.rst

# QM/MM with QUICK (GPU-accelerated, built-in API)
sander -O -i qmmm_quick.in -o qmmm_quick.out -p prmtop -c inpcrd

# QUICK standalone QM calculation
quick -i quick.in -o quick.out -p prmtop -c inpcrd

# QM/MM with TeraChem (TCPB client-server, recommended)
sander -O -i qmmm_tc.in -o qmmm_tc.out -p prmtop -c inpcrd

# QM/MM with xTB
sander -O -i qmmm_xtb.in -o qmmm_xtb.out -p prmtop -c inpcrd

# QM/MM with DFTB+
sander -O -i qmmm_dftbp.in -o qmmm_dftbp.out -p prmtop -c inpcrd

# Build reference system with antechamber
antechamber -i ligand.pdb -fi pdb -o ligand.mol2 -fo mol2 \
  -c bcc -s 2
parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod

# Parallel QM/MM
mpirun -np 8 sander.MPI -O -i qmmm.in -o qmmm.out -p prmtop -c inpcrd
```

## Key Namelist Variables

### &cntrl Namelist

```
ifqnt=1     ! Enable QM/MM coupled potential (0=off, 1=on)
```

### &qmmm Namelist - General QM/MM Control

```fortran
&qmmm
  ! QM region definition (choose one)
  qmmask=':753',              ! Mask string for QM atoms
  iqmatoms='1001,1002,1003',  ! Atom numbers (prmtop ordering) for QM atoms

  ! QM method
  qm_theory='PM3',            ! PM3/AM1/MNDO/PDDG-PM3/PM3CARB1/PM6/DFTB/DFTB2/DFTB3/EXTERN
  qmcharge=-2,                ! Total charge of QM region (default: 0)
  spin=1,                     ! Spin multiplicity (default: 1)

  ! Electrostatic interaction
  qmcut=8.0,                  ! QM/MM electrostatic cutoff (Angstroms)
  qmmm_int=1,                 ! QM-MM interaction scheme:
                              !   0=no electrostatic, 1=default, 2=Gaussian extra,
                              !   3=PM3/MM* reformulated, 5=mechanical embedding
  qm_ewald=1,                 ! 0=cutoff, 1=PME/Ewald (default), 2=frozen image charges
  qm_pme=1,                   ! 0=full Ewald, 1=PME approach (default)

  ! SHAKE
  qmshake=1,                  ! 0=no SHAKE on QM, 1=SHAKE QM H atoms (default)

  ! SCF convergence
  scfconv=1.0d-10,            ! SCF convergence criterion
  tight_p_conv=1,             ! 0=off, 1=use tight convergence
  diag_routine=1,             ! Diagonalization routine (0-3)
  verbosity=0,                ! Output verbosity level
  print_charges=0,            ! Print Mulliken charges
  itrmax=1000,                ! Maximum SCF iterations

  ! QM-QM analysis
  qmqmdx=1,                   ! QM-QM nonbonded exclusion list
  qmqm_anal=0,                ! QM-QM energy decomposition
  qm_mm_anal=0,               ! QM-MM pairwise energy decomposition

  ! Ewald control
  kmaxqx=8, kmaxqy=8, kmaxqz=8,  ! k-space vectors
  ksqmaxq=100,                    ! K^2 max for reciprocal space

  ! Switching function (NDDO methods)
  qmmm_switch=0,              ! 0=off, 1=use switching function
  r_switch_hi=8.0,            ! Upper boundary (defaults to qmcut)
  r_switch_lo=6.0,            ! Lower boundary (defaults to r_switch_hi-2)

  ! GB for QM/MM
  qmgb=2,                     ! 2=GB polarization in Fock matrix, 3=debugging

  ! Diagnostics
  writepdb=0,                 ! 0=off, 1=write qmmm_region.pdb
  printdipole=0,              ! 0=off, 1=QM dipole, 2=total dipole

  ! Adaptive solvent
  vsolv=0,                    ! 0=off, 1=simple, 2=fixed N, 3=fixed size
/
```

### Link Atom Variables

```fortran
&qmmm
  lnk_dis=1.09,               ! Distance from QM atom to link atom (default: 1.09 A)
  lnk_atomic_no=1,            ! Atomic number of link atom (default: 1 = Hydrogen)
  lnk_method=1,               ! 1=include MM valence terms, 2=exclude MM link pair
  adjust_q=2,                 ! 0=no adjustment, 1=nearest nlink, 2=all MM (default)
/
```

### External QM Program Interface

```fortran
&cntrl
  ifqnt=1,
/
&qmmm
  qmmask='@*',                ! All atoms QM (for pure QM)
  qmcharge=0,
  spin=1,
  qm_theory='EXTERN',         ! Enable external QM interface
  qmcut=999.0,                ! Include all MM charges (non-periodic)
/
&gau                          ! Gaussian-specific
  method='B3LYP',
  basis='6-31G*',
  scf_conv=1.0d-08,
  use_template=0,
/
&orc                          ! ORCA-specific
  method='B3LYP',
  basis='def2-SVP',
  scf_conv=1.0d-08,
/
&gms                          ! GAMESS-US-specific
  method='B3LYP',
  basis='6-31G*',
/
&tc                           ! TeraChem-specific
  method='B3LYP',
  basis='6-31G*',
/
```

### DFTB-specific Variables

```fortran
&qmmm
  qm_theory='DFTB3',          ! DFTB, DFTB2, or DFTB3
  dftb_disper=1,              ! Dispersion correction
  dftb_3rd_order='PA',        ! DFTB3 third-order extension
  dftb_chg=0,                 ! SCC charge
  dftb_telec=300.0,           ! Electronic temperature
  dftb_maxiter=100,           ! Max SCC iterations
/
```

## QM/MM Theory: Embedding Schemes

### Electrostatic Embedding (default)
The QM region is polarized by MM point charges. The MM charges enter the QM Hamiltonian:

```
H_QM/MM = sum_{m} Q_m [h_electron(R) - Z_q h_core(R)] + A/r^12 - B/r^6
```

### Mechanical Embedding (qmmm_int=5)
Classical point-charge interaction between QM and MM regions. No QM polarization.

### PM3/MM* Reformulated Interface (qmmm_int=3, qm_theory=PM3)
Reformulated QM core-MM charge potential with scaling, available for H, C, N, O QM atoms.

## Common Workflows

### Standard Semiempirical QM/MM MD

```bash
cat > qmmm_md.in << 'EOF'
&cntrl
  imin=0, nstlim=10000, dt=0.002,
  ntt=1, tempi=0.0, temp0=300.0,
  ntb=1, ntf=2, ntc=2,
  cut=8.0,
  ifqnt=1
/
&qmmm
  qmmask=':1-2',
  qmcharge=0,
  qm_theory='PM6',
  qmcut=8.0,
  qmshake=1,
/
EOF

sander -O -i qmmm_md.in -o qmmm_md.out -p system.prmtop -c system.inpcrd \
  -r qmmm_md.rst -x qmmm_md.nc
```

### QM/MM with External Gaussian

```bash
# Ensure 'g16' is in PATH
cat > qmmm_gau.in << 'EOF'
&cntrl
  imin=0, nstlim=5000, dt=0.001,
  ntt=1, tempi=0.0, temp0=300.0,
  ntb=0, ntf=1, ntc=1,
  cut=999.0,
  ifqnt=1
/
&qmmm
  qmmask=':1',
  qmcharge=0, spin=1,
  qm_theory='EXTERN',
  qmcut=999.0,
/
&gau
  method='B3LYP',
  basis='6-31G*',
  scf_conv=1.0d-08,
  scf_maxcyc=128,
/
EOF

sander -O -i qmmm_gau.in -o qmmm_gau.out -p system.prmtop -c system.inpcrd
```

### QM/MM with QUICK (GPU-Accelerated)

```bash
cat > qmmm_quick.in << 'EOF'
&cntrl
  imin=0, nstlim=10000, dt=0.001,
  ntb=1, ntf=1, ntc=1,
  cut=8.0,
  ifqnt=1
/
&qmmm
  qmmask=':1-2',
  qmcharge=0, spin=1,
  qm_theory='QUICK',
  qmcut=8.0,
/
&quick
  method='B3LYP',
  basis='6-31G*',
/
EOF

sander -O -i qmmm_quick.in -o qmmm_quick.out -p system.prmtop -c system.inpcrd
```

### Link Atom Setup for QM/MM Boundary

QM/MM boundary cut rules:
1. Cut non-polar bonds (C-C single bonds) -- avoid unsaturated or polar bonds
2. Link atoms are NOT placed between bonds to hydrogen
3. Only one link atom allowed per MM link pair atom
4. By default, hydrogen link atoms placed at 1.09 A along the QM-MM bond vector

```fortran
&qmmm
  ! Manual link atom control
  lnk_dis=1.09,               ! Default for C-H bond
  lnk_atomic_no=1,            ! Hydrogen link atom
  lnk_method=1,                ! Include MM valence terms crossing boundary
  adjust_q=2,                  ! Distribute charge excess over all MM atoms
/
```

## Reference Tables

### Built-in QM Methods

| qm_theory | Method | Notes |
|-----------|--------|-------|
| AM1 | Austin Model 1 | Basic semiempirical |
| PM3 | Parameterized Model 3 | Default |
| PM6 | Parameterized Model 6 | Improved PM3 |
| PM3CARB1 | PM3 with carbohydrate corrections | |
| MNDO | Modified Neglect of Diatomic Overlap | |
| PDDG-PM3 | PDDG reparameterized PM3 | |
| DFTB | SCC-DFTB | Requires SK files |
| DFTB2 | SCC-DFTB2 | mio-1-1 parameters |
| DFTB3 | DFTB3 with third-order | Requires dftb_3rd_order |
| EXTERN | External QM program | Gaussian, ORCA, etc. |
| QUICK | GPU-accelerated DFT | Built-in via API |

### Supported External QM Programs

| Program | Embedding | License | Interface |
|---------|-----------|---------|-----------|
| Gaussian | Elec+Mech | Commercial | File-based |
| ORCA | Elec+Mech | Free (academic) | File-based |
| GAMESS-US | Mechanical | Free (academic) | File-based |
| Q-Chem | Elec+Mech | Commercial | File-based |
| TeraChem | Elec+Mech | Commercial/Demo | File/TCPB client-server |
| QUICK | Elec+Mech | Free (AmberTools) | API library |
| NWChem | Mechanical | Free | File-based |
| ADF | Mechanical | Commercial | File-based |
| MRCC | Elec+Mech | Free (academic) | File-based |
| Fireball | Elec+Mech | Free | Linked library |
| xTB | Both | Free | Command-line |
| DFTB+ | Both | Free | Command-line |

### QM/MM Interaction Scheme (qmmm_int)

| Value | Description |
|-------|-------------|
| 0 | No electrostatic interaction (only VDW) |
| 1 | Default: QMcore-MMcharge + QMelectron-MMcharge |
| 2 | Same as 1 + Gaussian core-core terms (CHARMM style) |
| 3 | PM3/MM* reformulated core-charge potential |
| 4 | PM3/MMX2 with QM electron-fin-tuning |
| 5 | Mechanical embedding (classical point charges) |

### qm_ewald Options

| Value | Description |
|-------|-------------|
| 0 | Real-space cutoff (non-periodic default) |
| 1 | PME or Ewald sum (periodic default) |
| 2 | Frozen image charges (faster, approximate) |

## Worked Example: Enzyme Active Site QM/MM

```bash
# 1. Build reference force field (MM system)
antechamber -i substrate.pdb -fi pdb -o substrate.mol2 -fo mol2 -c bcc -s 2
parmchk2 -i substrate.mol2 -f mol2 -o substrate.frcmod

# 2. Prepare in LEaP
cat > build.in << 'EOF'
source leaprc.protein.ff19SB
source leaprc.gaff
loadAmberParams substrate.frcmod
SUB = loadMol2 substrate.mol2
prot = loadPdb protein.pdb
sys = combine {prot SUB}
solvateOct sys TIP3PBOX 12.0
addIons sys Na+ 0
saveAmberParm sys system.prmtop system.inpcrd
quit
EOF
tleap -f build.in

# 3. Run QM/MM MD (substrate + catalytic residues in QM region)
cat > qmmm_prod.in << 'EOF'
&cntrl
  imin=0, nstlim=100000, dt=0.001,
  ntt=3, gamma_ln=2.0, tempi=0.0, temp0=300.0,
  ntb=1, ntf=1, ntc=1,
  cut=9.0,
  ifqnt=1, ntwx=500, ntpr=100,
/
&qmmm
  qmmask=':1-2,:156',
  qmcharge=-1, spin=1,
  qm_theory='DFTB3',
  qmcut=9.0,
  qm_ewald=1, qm_pme=1,
  qmshake=0,
  dftb_3rd_order='PA',
  dftb_disper=1,
  writepdb=1, printdipole=0,
/
EOF

sander -O -i qmmm_prod.in -o qmmm_prod.out -p system.prmtop -c system.inpcrd \
  -r qmmm_prod.rst -x qmmm_prod.nc
```

## Key Takeaways

1. **Enable QM/MM with `ifqnt=1`** in `&cntrl` and define QM region via `qmmask` or `iqmatoms` in `&qmmm` namelist.
2. **Electrostatic embedding is default** -- the QM electron density is polarized by MM point charges. Use `qmmm_int=5` for mechanical embedding.
3. **Link atoms** are automatically added at QM/MM boundary cuts. Cut non-polar C-C bonds; avoid cutting polar/unsaturated bonds or C-H bonds.
4. **For external QM programs**, set `qm_theory='EXTERN'` and provide the program-specific namelist (`&gau`, `&orc`, `&tc`, etc.). The QUICK library interface is recommended for GPU-accelerated DFT.
5. **PME is not supported with external QM programs** -- use non-periodic simulations with a large cutoff (`qmcut` > system size) for external ab initio/DFT QM/MM.

## Connects To

- **Chapter 3: Force Fields** - MM force field parameters for QM/MM regions
- **Chapter 5: System Setup with LEaP** - Building prmtop files with antechamber
- **Chapter 6: Antechamber** - Creating reference MM parameters for QM regions
- **Chapter 12: Running MD** - Integration with sander MD engine
- **Chapter 17: Metal Ion Modeling** - MCPB.py as alternative to QM/MM for metal sites
- **Chapter 20: Implicit Solvent** - GB for QM/MM (qmgb), PB implicit solvent
- AMBER tutorials: Built-in QM/MM, external QM/MM interfaces