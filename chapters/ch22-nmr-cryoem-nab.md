# Chapter 22: NMR, CryoEM, SAXS & NAB

## Core Commands & Syntax

### NMR Restraint Setup (DISANG file)

```bash
# sander with NMR restraints
sander -O -i mdin -o mdout -p prmtop -c inpcrd -r rst7 \
       -x mdcrd -ref refcrd
```

mdin entry:
```
&cntrl
  nmropt=1,      ! Enable NMR restraint reading
  ntx=1, ntb=0,
  nstlim=10000, dt=0.001,
  ntt=3, gamma_ln=5.0,
  cut=12.0,
/
```

### NMR Restraint Format (DISANG file)

```
&rst
  iat(1)=15, iat(2)=48,    ! Atoms defining the restraint
  r1=0.0, r2=1.5, r3=2.5, r4=4.0,  ! Distance well (A)
  rk2=20.0, rk3=20.0,      ! Force constants (kcal/mol-A2)
  nstep1=0, nstep2=10000,  ! Steps to apply
  irstyp=0,                 ! 0=absolute, 1=relative
  ialtd=0,                  ! 0=standard, 1=flattened
/
&rst iat(1)=0, /           ! End of restraints
```

Restraint types determined by atom count:
- `iat(3)=0` -- Distance restraint
- `iat(4)=0` -- Angle restraint
- `iat(5)=0` -- Torsional or J-coupling restraint
- `iat(6)=0` -- Plane-point angle restraint
- `iat(7)=0` -- Generalized distance (6 atoms)
- All 8 atoms nonzero -- Plane-plane angle or generalized distance (8 atoms)

### J-Coupling Restraints

```
&rst
  iat(1)=7, iat(2)=9, iat(3)=15, iat(4)=17,
  rjcoef(1)=6.5, rjcoef(2)=-1.5, rjcoef(3)=1.8,  ! Karplus: A,B,C
  r2=5.0, r3=5.0,
  rk2=1.0, rk3=1.0,
/
```

### Chemical Shift Restraints

```
&cntrl
  nmropt=1,
  icshft=1,     ! Chemical shift restraints
  icshft_radius=0.5,  ! Radius for shift calculation
/
```

### Natural-Language Restraint Syntax (sander only)

```
restraint = "distance( (45) (49) )"
restraint = "angle(:21@C5' :21@C4' 108)"
restraint = "torsion[-1,-1,-1, com(67, 68, 69)]"
restraint = "angle( -1, plane(81, 85, 87, 90) )"
restraint = "coordinate(distance(:5@C3',:6@O5'),-1.0,distance(134,-1),1.0)"
```

### NOESY Volume Restraints

```
&cntrl
  nmropt=2,     ! Enable NOESY volume restraints
/
```

```
&noeexp
  npeak(1)=150,  ! # peaks at mixing time 1
  emix(1)=0.05,  ! Mixing time in seconds
  omega=600.0,   ! Spectrometer frequency (MHz)
  taurot=5.0,    ! Rotational tumbling time (ns)
  taumet=0.02,   ! Methyl correlation time (ns)
  id2o=0,        ! 0=include exchangeable, 1=exclude
  oscale=1.0e-5, ! Scaling factor
/
```

### CryoEM: EMAP Restraints

```bash
# EMAP-constrained SGLD simulation
sander -O -i emap.in -o emap.out -p prmtop -c inpcrd
```

emap.in:
```
&cntrl
  ntx=1, ntb=0, nstlim=100000, imin=0,
  ntc=2, ntf=2, cut=9.0, dt=0.001,
  ntt=3, gamma_ln=10.0,
  igb=0, isgld=1,              ! SGLD
  iemap=1,                      ! Turn on EMAP
/
&emap
  mapfile='data/map.ccp4',      ! EM density map
  atmask=':1-20',               ! Restrained atoms
  fcons=0.1,                    ! Force constant (kcal/g)
  move=1,                       ! Map can move with domain
  ifit=1,                       ! Rigid fit first
  resolution=2.0,               ! Map resolution (A)
  mapfit='final_map.ccp4',
  molfit='final_struct.pdb',
/
```

### FRET Restraints

```bash
# Step 1: Add pseudo atoms to topology
python3 placeAVmp.py -p protein_noH.pdb -o protein_PA.pdb \
        -j fret_config.json --chi2 'C3 chi2'

# Step 2: Generate DISANG restraints
python3 FRETrest.py -t system.prmtop -r equil.rst7 \
        -j fret_config.json --chi2 'C3 chi2' \
        --fout prod_0001.f --restout prod_0000.restrt \
        --force 50 --resoffset -1
```

### KMMD: Kernel Machine MD

```bash
# mdin setup
iextpot = 1,
/
&extpot
  extprog='kmmd',
  json='../kmmd_test.json',
/
```

KMMD JSON configuration defines dihedral angles to track:
```json
{
"DB_fileList": "cache_ffenes.txt",
"ref_pdb": "system.pdb",
"sigma2": 0.1,
"dh_atnames": [
  {"alpha": "O3' P O5' C5'"},
  {"nu0": "C4' O4' C1' C2'"}
],
"dh_byres": [
  {"nu0 1": ""},
  {"nu0 2": ""}
]
}
```

### NAB: Nucleic Acid Builder

```bash
# Compile NAB code to C
nabc script.nab

# Compile and link
gcc -o script script.c -I$AMBERHOME/include -L$AMBERHOME/lib -lsff -lm

# Run NAB program
./script
```

Basic NAB program structure:
```c
molecule m;
m = bdna("gcgc");
putpdb("dna.pdb", m);
```

## Key Namelists / Input Files

### &rst Namelist (per restraint)

| Variable | Description | Default |
|----------|-------------|---------|
| iat(1-8) | Atom numbers (or negative for groups) | 0 |
| r1,r2,r3,r4 | Restraint well boundaries (A or deg) | 0 |
| rk2,rk3 | Force constants (kcal/mol-A2 or kcal/mol-rad2) | 0 |
| nstep1,nstep2 | Step range for restraint | 0,0 |
| irstyp | 0=absolute, 1=relative | 0 |
| ialtd | 0=standard, 1=flattened (large viol.) | 0 |
| ifvari | 0=constant, >0=varying | 0 |
| rjcoef(1-3) | Karplus A,B,C for J-coupling | 0,0,0 |
| r0,k0 | Simplified parabolic well | -- |
| restraint | Natural-language restraint string | '' |

### &emap Namelist (per EM restraint)

| Variable | Description | Default |
|----------|-------------|---------|
| mapfile | Map file (.map, .ccp4, .mrc) or '' | '' |
| atmask | Atom mask for restrained atoms | ':*' |
| fcons | Restraint force constant (kcal/g) | 0.05 |
| move | 0=fixed, >0=translate+rotate map | 0 |
| resolution | Map resolution (A), <0=boundary | 2.0 |
| ifit | 0=none, 1=rigid fit map, 2=fit coords | 0 |
| grids | Grid search params for fitting | 1,1,1,1,1,1 |

### NOESY &noeexp Namelist

| Variable | Description |
|----------|-------------|
| npeak(imix) | Number of peaks per mixing time |
| emix(imix) | Mixing time (seconds) |
| ihp, jhp | Atom numbers for cross-peak |
| aexp | Experimental integrated intensity |
| omega | Spectrometer frequency (MHz) |
| taurot | Rotational tumbling time (ns) |
| taumet | Methyl correlation time (ns) |
| id2o | Exclude exchangeable protons |
| oscale | Scaling factor |

## Common Workflows

### Distance Restraints for NMR Refinement

```bash
# 1. Prepare DISANG file
cat > RST.dist << 'EOF'
&rst
  iat(1)=15, iat(2)=48,
  r1=0.0, r2=1.5, r3=2.5, r4=4.0,
  rk2=20.0, rk3=20.0,
  nstep1=0, nstep2=10000,
/
&rst iat(1)=0, /
EOF

# 2. Run sander with NMR
sander -O -i mdin -o mdout \
  -p prmtop -c inpcrd -r rst7 \
  -ref refcrd -x mdcrd
```

### EM Fitting with CryoEM Map

```
# 1. Rigid fitting then flexible refinement
&cntrl
  imin=0, nstlim=50000, dt=0.001,
  ntt=3, gamma_ln=5.0, ntb=0,
  cut=9.0, igb=1,
  iemap=1,
/
&emap
  mapfile='experimental.ccp4',
  atmask=':1-150',
  fcons=0.1,
  ifit=1,         ! Rigid fit first
  resolution=8.0, ! Low-res map
  mapfit='fitted.ccp4',
  molfit='fitted.pdb',
/
```

### NAB: Build DNA Duplex

```c
// dna_duplex.nab
molecule m;
int i;

// Build B-DNA
m = bdna("gcgcaattcgcg");

// Save to PDB
putpdb("dna_duplex.pdb", m);

// Minimize with sff
m = minimize(m, "sff_min.in");
```

## Reference Tables

### NMR Restraint Types

| Type | Atoms | Variables | Functional Form |
|------|-------|-----------|----------------|
| Distance | 2 | r1-r4, rk2, rk3 | Flat-bottom parabolic well |
| Angle | 3 | r1-r4, rk2, rk3 | Flat-bottom parabolic well |
| Torsion | 4 | r1-r4, rk2, rk3 | Flat-bottom parabolic well |
| J-coupling | 4 | rjcoef, r2-r3, rk2-rk3 | Karplus relation |
| Plane-point | 5 | r1-r4, rk2, rk3 | Angle between plane normal and vector |
| Plane-plane | 8 | r1-r4, rk2, rk3 | Angle between two plane normals |
| Generalized distance | 4,6,8 | rstwt + r1-r4 | Weighted linear combination |

### NAB Key Functions

| Function | Purpose |
|----------|---------|
| `bdna(seq)` | Build B-DNA duplex |
| `arna(seq)` | Build A-RNA duplex |
| `putpdb(file, mol)` | Write PDB file |
| `getpdb(file)` | Read PDB file |
| `minimize(mol, input)` | Energy minimization |
| `dynamics(mol, input)` | Molecular dynamics |
| `putbnd()`, `putang()` | Set bond/angle parameters |
| `putvdw()`, `putels()` | Set VDW/electrostatic params |
| `addions(mol, "Na+", 10)` | Add counterions |
| `setbox(mol, "rect", 10.0)` | Set periodic box |
| `setsugar(mol, "dna")` | Set sugar pucker |

### AI/ML Methods in Amber

| Method | Program | Use Case |
|--------|---------|----------|
| KMMD | sander | Dihedral energy correction via kernel machine |
| DPRc | sander | QM/MM+ML neural network correction |
| Torch PBSA | pbsa | ML-accelerated PB surface construction |
| External library | sander | Custom ML potentials via iextpot |

## Worked Example

### EMAP Restrained Flexible Fitting

From the manual, fitting a protein domain into a cryo-EM density map with separate domain restraints:

```bash
# Domain 1 (rigid map): residues 1-20
&emap
  mapfile='data/1gb1.ccp4',
  atmask=':1-20',
  fcons=0.1,
  move=1,                    ! Map follows domain
  ifit=1,                    ! Rigid fit first
  mapfit='scratch/gb1n_1.ccp4',
  molfit='scratch/gb1n_1.pdb',
/
# Domain 2 (fixed map): residues 22-37
&emap
  mapfile='data/1gb1.pdb',
  atmask=':22-37',
  fcons=0.1, move=0,         ! Map stays fixed
  ifit=1,
  mapfit='scratch/gb1h_1.ccp4',
  molfit='scratch/gb1h_1.pdb',
/
# Domain 3 (self-restraint): residues 41-56
&emap
  mapfile='',                 ! Auto-generated from coords
  atmask=':41-56',
  fcons=0.1, move=1,
  ifit=1,
/
# Boundary map for finite system
&emap
  mapfile='', atmask=':*',
  fcons=-1.0,                 ! Negative = boundary mode
  resolution=-4,              ! Negative = create boundary
  grids=20,20,20,1,1,1,      ! Rectangular shape
/
```

## Key Takeaways

1. **NMR restraints use the DISANG file format** with `&rst` namelists. Each restraint defines a flat-bottom parabolic well controlled by r1-r4 and rk2-rk3. Use `makeDIST_RST` to generate from simpler input.
2. **EMAP fits structures into cryo-EM density maps** via energy restraints. Use `move=1` for domain motion, `ifit=1` for initial rigid fitting, and `isgld=1` for large conformational changes.
3. **FRET restraints** are added via pseudo-atom placement (placeAVmp.py) and DISANG restraint generation (FRETrest.py). Requires the Olga software for FRET configuration.
4. **KMMD adds ML corrections to force fields** by learning energy differences from ab initio data. Define dihedral angles of interest in a JSON file and use `iextpot=1`.
5. **NAB is a legacy scripting language** for nucleic acid building and molecular mechanics. Moved to a separate GitHub repository (https://github.com/dacase/nabc). Use `nabc` to compile, `libsff` for force field routines.

## Connects To

- [ch07](ch07-md-engines-input.md) -- MD engines and mdin reference
- [ch08](ch08-minimization-relaxation.md) -- Minimization and restraint setup
- [ch14](ch14-enhanced-sampling.md) -- Enhanced sampling (SGLD used with EMAP)
- [ch19](ch19-qmmm.md) -- QM/MM (DPRc ML corrections)
- [ch20](ch20-implicit-solvent.md) -- Implicit solvent (used in NAB sff)