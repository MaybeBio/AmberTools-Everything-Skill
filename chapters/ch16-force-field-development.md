# Chapter 16: Force Field Development with mdgx, py_resp, pyMSMT

## Core Commands & Syntax

### mdgx: Charge Fitting (RESP/IPolQ)
```bash
mdgx -O -i mdin -o mdout -p prmtop -c inpcrd -r rst7
```

### mdgx: Bonded Term Fitting
```bash
mdgx -O -i mdin -o mdout -p prmtop -c inpcrd -d fitted_params.frcmod
```

### mdgx: Configuration Sampling
```bash
mdgx -O -i configs.in -o configs.out -p prmtop -c inpcrd
```

### py_resp: RESP Charge Derivation
```bash
# Single conformation
py_resp.py -i input.gesp -o output.chg -t charge -e esp.dat -s esp.out

# With mol2 input
py_resp.py -i mol.mol2 -o mol.resp.chg -t charge -e esp.dat

# With polarizabilities
py_resp.py -i input.gesp -o output.chg -ip polariz.in -t charge -e esp.dat
```

### py_resp: Two-Stage RESP
```bash
# Stage 1
py_resp.py -i resp1.in -o resp1.chg -t resp1.chg -e esp.dat -q qin

# Stage 2 (with equivalencing)
py_resp.py -i resp2.in -o resp2.chg -t resp2.chg -e esp.dat -q resp1.chg
```

### pyresp_gen.py: Generate py_resp Input
```bash
pyresp_gen.py -i mol.mol2 -o mol.gesp
```

### antechamber RESP workflow
```bash
antechamber -i mol.pdb -fi pdb -o mol.mol2 -fo mol2 -c bcc -s 2
antechamber -i mol.mol2 -fi mol2 -o mol.ac -fo ac
prepgen -i mol.ac -o mol.prepin -m mol.prepc -rn MOL
parmchk2 -i mol.prepin -f prepi -o mol.frcmod
```

### pyMSMT: Metal Site Modeling
```bash
MCPB.py -i MCPB.in -s 1   # Step 1: generate initial files
MCPB.py -i MCPB.in -s 2   # Step 2: QM calculation
MCPB.py -i MCPB.in -s 3   # Step 3: force field parameterization
MCPB.py -i MCPB.in -s 4   # Step 4: generate tleap input
```

## Key Namelists / Input Files

### mdgx &fitq namelist (Charge Fitting)
```
&fitq
  RespPhi   Conf12/pcm12.cube, prmtop, 1.0,
  RespPhi   Conf13/pcm13.cube, prmtop, 1.0,
  pnrg    2.0,
  nfpt    15000,
  minqwt 175.0,
  EqualizeQ '@H1,H2'
  EqualizeQ '@Cl2,Cl3'
  MinimizeQ = '@E*'
  EPRules     frag.xpt
  ConfFile    f6xp.pdb
&end
```

### &fitq key variables
| Variable | Description | Default |
|----------|-------------|---------|
| `RespPhi` | Gaussian cubegen file, topology, weight | None |
| `IPolQPhi` | Vacuum cubegen, solvent cubegen, topology, weight | None |
| `TotalQ (qtot)` | Total charge constraint | 0.0 |
| `MinimizeQ (minq)` | Mask of atoms restrained to zero | None |
| `EqualizeQ (equalq)` | Mask of atoms with equal charges | None |
| `MinQWeight (minqwt)` | Weight for charge restraints | None |
| `FitPoints (nfpt)` | Number of fitting points | 1000 |
| `ProbeSig (psig)` | LJ sigma of solvent probe | 3.16435 (TIP4P O) |
| `ProbeEps (peps)` | LJ epsilon of solvent probe | 0.16275 (TIP4P O) |
| `ProbeArm (parm)` | Probe arm reach | 0.9572 (TIP O-H) |
| `StericLimit (pnrg)` | Max LJ energy for point inclusion | 3.0 kcal/mol |
| `Proximity (flim)` | Min distance between fitting points | 0.4 A |
| `EPRules (eprules)` | Output Virtual Sites rule file | None |

### mdgx &ipolq namelist (IPolQ Polarized Charge Development)
```
&ipolq
  SoluteMol = ':1-5',
  FrameRate = 1000,
  FrameCount = 100,
  EqStepCount = 50000,
  Blocks = 4,
  Verbose = 1,
/
```

### &ipolq key variables
| Variable | Description | Default |
|----------|-------------|---------|
| `SoluteMol (solute)` | Mask for solute molecule | Required |
| `FrameRate (ntqs)` | Steps between SRFP snapshots | 1000 |
| `FrameCount (nqframe)` | Number of frames for SRFP | 10 |
| `EqStepCount (nsteqlim)` | Equilibration steps before collection | 10000 |
| `Blocks (nblock)` | Blocks for convergence estimation | 4 |
| `EConverge (econv)` | Convergence tolerance | not implemented |

### mdgx &param namelist (Bonded Term Fitting)
```
&param
  System    prmtop, coords.rst7, -150.234567,
  System    prmtop, coords2.rst7, -150.236789,
  FitBonds = 1,
  FitAngles = 1,
  FitTorsions = 1,
  FitCmaps = 1,
  BondRest = 100.0,
  AngleRest = 50.0,
  DihedralRest = 100.0,
  CmapRest = 500.0,
  ReportAll = 1,
  ShowProgress = 1,
  EnergyUnits = 'Hartree',
  ParmTitle = 'My Custom Force Field',
/
```
Note: `System` keywords are followed by `prmtop_file coord_file energy_value`.

### &param key variables
| Variable | Description |
|----------|-------------|
| `System (sys)` | prmtop, coords, target energy (Hartrees) |
| `FitBonds (bonds)` | 1=fit bond stiffnesses |
| `FitAngles (angles)` | 1=fit angle stiffnesses |
| `FitTorsions (torsions)` | 1=fit torsion amplitudes |
| `FitCmaps (cmaps)` | 1=fit CMAP surfaces |
| `FitBondEq (bondeq)` | 1=fit bond equilibrium lengths |
| `FitAnglEq (angleq)` | 1=fit angle equilibrium values |
| `BondRest (brst)` | Harmonic restraint on bond stiffness |
| `AngleRest (arst)` | Harmonic restraint on angle stiffness |
| `DihedralRest (hrst)` | Harmonic restraint on torsion amplitudes |
| `CmapRest (mrst)` | Restraint on CMAP surface smoothness |
| `EnergyUnits (eunits)` | 'Hartree', 'kJ', or 'j' |
| `ReportAll (repall)` | 1=write all params (parm.dat), 0=frcmod style |
| `AccReport (accrep)` | Write accuracy report (MatLab format) |
| `ElimOutliers (elimsig)` | 1=remove outlier conformations |
| `ConfTol (ctol)` | Sigma threshold for outlier removal |

### Targeted restraint sub-commands
| Command | Purpose |
|---------|---------|
| `RestrainB (sbrst)` | Specific restraint on a bond: `type1 type2 Keq <val> Leq <val>` |
| `RestrainA (sarst)` | Specific restraint on an angle: `type1 type2 type3 Keq <val> Leq <val>` |
| `RestrainH (shrst)` | Specific restraint on torsion: `type1-4 period <N> weight <val> target <val>` |

### py_resp.py input file format (-i)
```
TITLE (line 1)
&cntrl
  nmol = 1,           ! number of structures
  iqopt = 1,          ! 1=zero init, 2=read from -q
  ihfree = 1,         ! 1=H not restrained (default)
  irstrnt = 1,        ! 0=harmonic, 1=hyperbolic (default)
  qwt = 0.0005,       ! charge restraint weight
  ioutopt = 1,        ! 1=write new ESP to -s
  ipol = 0,           ! 0=additive RESP, 5=pGM damping
  ipol = 5,           ! for polarizable models
  igdm = 1,           ! 1=use distributed pGM charges
  exc12 = 0,          ! include 1-2 interactions
  exc13 = 0,          ! include 1-3 interactions
  ipermdip = 0,       ! 0=RESP-ind, 1=RESP-perm
  pwt = 0.0005,       ! permanent dipole restraint weight
  virtual = 0,        ! 1=enable virtual dipoles
&end
wtmol (line 3, e.g., 1.0)
subtitle (line 4)
charge natom [natom_p] (line 5, e.g., 0 3)
atom1: element ivary [ivary_p] (line 6+)
...
(blank line)
intra-molecular charge constraints
(blank line)
inter-molecular charge constraints
(blank line)
multi-structure equivalencing
(blank line)
```

### py_resp.py key &cntrl variables
| Variable | Description | Default |
|----------|-------------|---------|
| `nmol` | Number of structures | 1 |
| `iqopt` | 1=zero init, 2=read from -q | 1 |
| `ihfree` | 0=restrain all, 1=H not restrained | 1 |
| `irstrnt` | 0=harmonic, 1=hyperbolic, 2=analysis only | 1 |
| `qwt` | Charge restraint weight | 0.0005 |
| `ipol` | 0=additive, 1-5=polarizable models | 0 |
| `ipermdip` | 0=RESP-ind, 1=RESP-perm | 0 |
| `pwt` | Permanent dipole restraint weight | 0.0005 |
| `virtual` | 1=enable virtual dipoles | 0 |
| `iquad` | 0=Buckingham, 1=Gaussian quadrupole | 1 |

## Common Workflows

### 1. RESP Charge Fitting with py_resp.py
**Step 1:** Generate ESP data from QM (Gaussian, GAMESS, etc.):
```bash
# Gaussian example
g16 < mol.gjf > mol.log
espgen -i mol.log -o mol.esp
```

**Step 2:** Generate py_resp input:
```bash
pyresp_gen.py -i mol.mol2 -o mol.gesp
```

**Step 3:** Run two-stage RESP:
```bash
# Stage 1: fit charges with weak restraint (qwt=0.0005)
py_resp.py -i mol.gesp -o stage1.chg -t stage1.chg -e mol.esp

# Stage 2: refit with stronger restraint (qwt=0.001), methyl H equivalenced
# Modify the .gesp file to equivalence methyl hydrogens (ivary=n)
py_resp.py -i mol_resp2.gesp -o final.chg -t final.chg -e mol.esp -q stage1.chg
```

### 2. Full Force Field Development with mdgx
**Step 1:** Generate conformations via &configs:
```
&configs
  count = 1000,
  maxcyc = 500, ncyc = 100,
  ...
  # Dihedral scanning restraints
&end
```
```bash
mdgx -O -i sample.in -p prmtop -c inpcrd -o sampled.out
```

**Step 2:** Run QM single-points on each conformation.

**Step 3:** Extract target energies and fit parameters:
```
&param
  System  prmtop, conf_0001.rst7, -150.234567,
  System  prmtop, conf_0002.rst7, -150.236789,
  ...
  FitBonds = 1, FitAngles = 1, FitTorsions = 1,
  BondRest = 100.0, AngleRest = 50.0, DihedralRest = 100.0,
  EnergyUnits = 'Hartree',
  ParmTitle = 'Custom Dihedral Parameters',
/
```
```bash
mdgx -O -i fit.in -p prmtop -c inpcrd -d fitted.frcmod
```

**Step 4:** Validate: load fitted.frcmod in LEaP and test against QM reference data.

### 3. Full RESP Charge Pipeline (antechamber + py_resp)
```bash
# 1. Generate mol2 and ESP
antechamber -i mol.pdb -fi pdb -o mol.mol2 -fo mol2 -c bcc
g16 < mol.gjf > mol.log
espgen -i mol.log -o mol.esp

# 2. Generate py_resp input
pyresp_gen.py -i mol.mol2 -o mol.gesp

# 3. Two-stage RESP
py_resp.py -i mol.gesp -o mol_resp1.chg -t mol_resp1.chg -e mol.esp
# Edit mol.gesp -> mol_resp2.gesp (equivalence methyl H, set qwt=0.001)
py_resp.py -i mol_resp2.gesp -o mol_resp2.chg -t mol_resp2.chg \
    -e mol.esp -q mol_resp1.chg

# 4. Generate prepi and frcmod
antechamber -i mol.mol2 -fi mol2 -o mol_charges.ac -fo ac \
    -c rc -cf mol_resp2.chg
prepgen -i mol_charges.ac -o mol.prepin -m mol.prepc -rn MOL
parmchk2 -i mol.prepin -f prepi -o mol.frcmod
```

### 4. pyMSMT Metal Site Parameterization
```bash
# Step 1: Generate initial files
MCPB.py -i MCPB.in -s 1

# Step 2: Run QM optimization/frequency (Gaussian)
MCPB.py -i MCPB.in -s 2

# Step 3: Force field parameterization
MCPB.py -i MCPB.in -s 3

# Step 4: Generate tleap input
MCPB.py -i MCPB.in -s 4
```

## Reference Tables

### mdgx run modes
| Namelist | Mode | Purpose |
|----------|------|---------|
| `&fitq` | Charge fitting | RESP/IPolQ charge derivation from QM ESP |
| `&ipolq` | IPolQ protocol | Automate solvent reaction field potential calculation |
| `&param` | Bonded fitting | Fit bond/angle/torsion/CMAP to QM energies |
| `&configs` | Conformation sampling | Generate thousands of restrained conformations |

### mdgx output files
| Flag | Purpose |
|------|---------|
| `-d` | Parameter output (frcmod or parm.dat format) |
| `-o` | mdout: extensive fitting diagnostics |
| `EPRules` | Virtual Sites rule file with fitted charges |
| `ConfFile` | PDB of first conformation with virtual sites |
| `AccReport` | MatLab-format accuracy report |

### RESP stages
| Stage | qwt | Description |
|-------|-----|-------------|
| Stage 1 | 0.0005 | Weak restraint, fit all atoms |
| Stage 2 | 0.001 | Stronger restraint, equivalence methyl H/C |

### py_resp polarizable models (ipol)
| ipol | Model |
|------|-------|
| 0 | Additive RESP (default) |
| 1 | Applequist without damping |
| 2 | Tinker exponential damping |
| 3 | Exponential damping |
| 4 | Linear damping |
| 5 | pGM damping (default for polarizable) |

### py_resp permanent dipole models (ipermdip)
| ipermdip | Model |
|----------|-------|
| 0 | RESP-ind (no permanent dipoles) |
| 1 | RESP-perm (calculate permanent dipoles) |

## Worked Example

### Water RESP Charge Fitting

**RESP (additive model):**
```
resp for water
&cntrl
 nmol = 1, iqopt = 1, ihfree = 1,
 qwt = 0.0005, ioutopt = 1, ipol = 0, ipermdip = 0
&end
1.0
water
    0    3
    8    0     ! oxygen, independent fit
    1    0     ! hydrogen 1, independent fit
    1    2     ! hydrogen 2, equivalenced to atom 2
```

**RESP-perm (polarizable with permanent dipoles):**
```
resp-perm for water
&cntrl
 nmol = 1, iqopt = 1, ihfree = 1,
 qwt = 0.0005, ioutopt = 1, ipol = 5, igdm = 1,
 exc12 = 0, exc13 = 0, ipermdip = 1, pwt = 0.0005, virtual = 0
&end
1.0
water
     0   3 4           ! 4 permanent dipoles total
     8   0 0 1         ! O: charge ivary=0, dip ivary=0,1
     1   0 0           ! H1: charge ivary=0, dip ivary=0
     1   2 3           ! H2: charge ivary=2 (equiv to H1), dip ivary=3 (equiv)
```

**Run:**
```bash
py_resp.py -i water.in -o water.chg -t water.chg -e water.esp
```

## Key Takeaways

1. **mdgx** is the primary AmberTools program for force field development: charge fitting via `&fitq`, IPolQ via `&ipolq`, bonded parameter fitting via `&param`, and configuration sampling via `&configs`.
2. **py_resp.py** performs two-stage RESP fitting: Stage 1 (qwt=0.0005) fits all atoms weakly; Stage 2 (qwt=0.001) equivalences methyl groups and refits with stronger restraint.
3. **mdgx &param** fits bond/angle/torsion/CMAP parameters to QM target energies by solving a linear least-squares problem with optional harmonic restraints.
4. **pyMSMT (MCPB.py)** provides a 4-step pipeline for metal site force field parameterization: generate input, QM calculation, FF parameterization, and tleap input generation.
5. The **antechamber + py_resp + prepgen + parmchk2** pipeline produces a complete residue definition (prepi + frcmod) with RESP charges for any molecule.

## Connects To
- **Chapter 2 (Force Fields):** Force field parameter formats, frcmod, prepi, lib files
- **Chapter 3 (LEaP):** Loading custom residue libraries and frcmod files
- **Chapter 17 (Metal Ion Modeling):** MCPB.py, pyMSMT, metal site parameters
- Quantum chemistry software: Gaussian, GAMESS, ORCA, QUICK