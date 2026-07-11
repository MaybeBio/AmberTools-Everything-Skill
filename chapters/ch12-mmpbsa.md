# Chapter 12: MM-PBSA and MM-GBSA Binding Free Energy

## Core Commands & Syntax

### Serial
```bash
MMPBSA.py -O -i mmpbsa.in -o FINAL_RESULTS_MMPBSA.dat \
    -sp com_solvated.prmtop \
    -cp com.prmtop -rp rec.prmtop -lp lig.prmtop \
    -y prod.nc
```

### Parallel (MPI)
```bash
mpirun -np 8 MMPBSA.py.MPI -O -i mmpbsa.in \
    -o FINAL_RESULTS_MMPBSA.dat \
    -sp com_solvated.prmtop \
    -cp com.prmtop -rp rec.prmtop -lp lig.prmtop \
    -y prod.nc
```

### Preparatory: ante-MMPBSA.py
```bash
ante-MMPBSA.py -p com_solvated.prmtop \
    -c com.prmtop -r rec.prmtop -l lig.prmtop \
    -s ":WAT,:Na+,:Cl-" -m ":LIG"
```

### Command-line options
| Flag | Purpose |
|------|---------|
| `-i` | Input file |
| `-o` | Output file (default: FINAL_RESULTS_MMPBSA.dat) |
| `-sp` | Solvated complex topology |
| `-cp` | Dry complex topology |
| `-rp` | Receptor topology |
| `-lp` | Ligand topology |
| `-y` | Input trajectory (comma-separated list) |
| `-yr` | Receptor trajectory (multiple trajectory protocol) |
| `-yl` | Ligand trajectory (multiple trajectory protocol) |
| `-do` | Decomposition output (default: FINAL_DECOMP_MMPBSA.dat) |
| `-eo` | CSV energy output for every frame |
| `-deo` | CSV decomposition energy output for every frame |
| `-make-mdins` | Generate mdin files and quit |
| `-use-mdins` | Use existing mdin files |
| `-rewrite-output` | Reparse existing output files |
| `-xvvfile` | XVV file for 3D-RISM |
| `-prefix` | Prefix for temp files |
| `--clean` | Clean temp files from previous run |

## Key Namelist / Input File

### MMPBSA input file format
```
&general
    startframe=5, endframe=100, interval=5,
    verbose=2, keep_files=0,
    strip_mask=":WAT:CL:CIO:CS:IB:K:LI:MG:NA:RB",
    netcdf=0,
    entropy=0,
/
&gb
    igb=5, saltcon=0.150,
    surften=0.0072, surfoff=0.0,
    molsurf=0,
/
&pb
    istrng=0.15, fillratio=4.0,
    indi=1.0, exdi=80.0,
    radiopt=1, prbrad=1.4,
    linit=1000, scale=2.0,
    inp=2,
    cavity_surften=0.0378, cavity_offset=-0.5692,
/
&decomp
    idecomp=2, dec_verbose=3,
    print_res="20, 40-80, 200",
    csv_format=1,
/
&alanine_scanning
    mutant_only=0,
/
&nmode
    nmstartframe=1, nmendframe=100, nminterval=1,
    maxcyc=10000, drms=0.001,
    nmode_igb=1, nmode_istrng=0.0,
    dielc=1.0,
/
```

### &general namelist variables
| Variable | Description | Default |
|----------|-------------|---------|
| `startframe` | First frame to extract | 1 |
| `endframe` | Last frame to extract | 9999999 |
| `interval` | Frame stride | 1 |
| `verbose` | 0=diff only, 1=all terms, 2=+bonded | 1 |
| `keep_files` | 0=none, 1=traj+mdout, 2=all temp | 1 |
| `netcdf` | Use NetCDF temp trajectories | 0 |
| `entropy` | Quasi-harmonic entropy via ptraj | 0 |
| `strip_mask` | Atoms to strip from solvated traj | `:WAT:CL:CIO:CS:IB:K:LI:MG:NA:RB` |
| `receptor_mask` | Receptor residues in complex | Auto-guessed |
| `ligand_mask` | Ligand residues in complex | Auto-guessed |
| `use_sander` | Force sander over mmpbsa_py_energy | 0 |
| `full_traj` | Write full combined trajectories (parallel) | 0 |
| `debug_printlevel` | 0=silent, 1=full traceback | 0 |

### &gb namelist variables (MM-GBSA)
| Variable | Description | Default |
|----------|-------------|---------|
| `igb` | GB model: 1, 2, 5, 7, 8, 66 | 66 |
| `saltcon` | Salt concentration (M) | 0.0 |
| `surften` | Surface tension (kcal/mol/A^2) | 0.0072 |
| `surfoff` | Surface area offset | 0.0 |
| `molsurf` | Use molsurf (1) vs LCPO (0) | 0 |
| `probe` | Probe radius for molsurf | 1.4 |
| `msoffset` | Atomic radii offset for molsurf | 0 |
| `ifqnt` | QM/MM: 0=off, 1=on | 0 |
| `qm_residues` | QM residues (comma-delimited) | None |
| `qm_theory` | Semi-empirical Hamiltonian | None |
| `qmcharge_com` | QM charge for complex | 0 |
| `qmcharge_lig` | QM charge for ligand | 0 |
| `qmcharge_rec` | QM charge for receptor | 0 |
| `qmcut` | QM/MM charge interaction cutoff | 9999.0 |

### &pb namelist variables (MM-PBSA)
| Variable | Description | Default |
|----------|-------------|---------|
| `inp` | Nonpolar optimization: 1, 2 | 2 |
| `istrng` | Ionic strength (M) | 0.0 |
| `fillratio` | Grid/solute dimension ratio | 4.0 |
| `scale` | Grid spacing reciprocal | 2.0 |
| `indi` | Internal dielectric constant | 1.0 |
| `exdi` | External dielectric constant | 80.0 |
| `radiopt` | Radii setup: 0=prmtop, 1=pre-computed | 1 |
| `prbrad` | Solvent probe radius (A) | 1.4 |
| `linit` | Max linear PB iterations | 1000 |
| `cavity_surften` | Cavity surface tension | 0.0378 |
| `cavity_offset` | Cavity offset | -0.5692 |
| `sander_apbs` | Use APBS instead of PBSA | 0 |
| `memopt` | Membrane protein support | 0 |
| `emem` | Membrane dielectric | 1.0 |
| `mthick` | Membrane thickness (A) | 40.0 |
| `mctrdz` | Membrane center Z (A) | 0.0 |
| `poretype` | Auto pore finding | 1 |

### &decomp namelist variables
| Variable | Description | Default |
|----------|-------------|---------|
| `idecomp` | 1=per-residue (1-4 to internal), 2=per-residue (1-4 split), 3=pairwise (1-4 internal), 4=pairwise (1-4 split) | None (required) |
| `dec_verbose` | 0=delta total, 1=+sidechain, 2=all total, 3=all+sidechain | 0 |
| `print_res` | Residues to print (e.g., "1,3-10,15") | All |
| `csv_format` | 0=ASCII, 1=CSV with SEM | 1 |

### &alanine_scanning namelist
| Variable | Description | Default |
|----------|-------------|---------|
| `mutant_only` | 0=do mutant+original, 1=mutant only | 0 |

### &nmode namelist (Normal Mode Entropy)
| Variable | Description | Default |
|----------|-------------|---------|
| `nmstartframe` | Start frame for nmode | 1 |
| `nmendframe` | End frame for nmode | 1000000 |
| `nminterval` | Frame stride for nmode | 1 |
| `maxcyc` | Max minimization cycles | 10000 |
| `drms` | Energy gradient convergence | 0.001 |
| `nmode_igb` | GB model: 0=vacuum, 1=HCT | 1 |
| `nmode_istrng` | Ionic strength for nmode (M) | 0.0 |
| `dielc` | Distance-dependent dielectric | 1.0 |

## Common Workflows

### Workflow 1: Standard MM-GBSA binding free energy
```bash
# 1. Build topology files in tleap
cat > mmpbsa_leap.in <<EOF
source leaprc.protein.ff14SB
source leaprc.water.tip3p
loadAmberParams LIG.frcmod
LIG = loadMol2 LIG.mol2
receptor = loadPDB receptor.pdb
complex = combine {receptor LIG}
set default PBRadii mbondi2
saveAmberParm LIG lig.prmtop lig.crd
saveAmberParm receptor rec.prmtop rec.crd
saveAmberParm complex com.prmtop com.crd
solvateOct complex TIP3PBOX 15.0
saveAmberParm complex com_solvated.prmtop com_solvated.crd
quit
EOF
tleap -f mmpbsa_leap.in

# 2. Run MD, then MMPBSA.py
cat > mmpbsa.in <<EOF
&general
    startframe=5, endframe=1000, interval=5,
    verbose=2, keep_files=0,
/
&gb
    igb=5, saltcon=0.150,
/
EOF

MMPBSA.py -O -i mmpbsa.in -o FINAL_RESULTS.dat \
    -sp com_solvated.prmtop \
    -cp com.prmtop -rp rec.prmtop -lp lig.prmtop \
    -y prod.nc
```

### Workflow 2: MM-PBSA with decomposition analysis
```bash
cat > mmpbsa_pb.in <<EOF
&general
    startframe=1, endframe=500, interval=5,
    verbose=2,
/
&pb
    istrng=0.15, fillratio=4.0,
    radiopt=1, indi=1.0, exdi=80.0,
    inp=2,
/
&decomp
    idecomp=2, dec_verbose=3,
    print_res="1-300",
    csv_format=1,
/
EOF

MMPBSA.py -O -i mmpbsa_pb.in \
    -cp com.prmtop -rp rec.prmtop -lp lig.prmtop \
    -y prod.nc
```

### Workflow 3: MM-GBSA with normal mode entropy
```bash
cat > mmpbsa_nmode.in <<EOF
&general
    startframe=5, endframe=100, interval=5,
    verbose=2, keep_files=2,
/
&gb
    igb=5, saltcon=0.150,
/
&nmode
    nmstartframe=2, nmendframe=20, nminterval=2,
    maxcyc=50000, drms=0.0001,
/
EOF

MMPBSA.py -O -i mmpbsa_nmode.in \
    -cp com.prmtop -rp rec.prmtop -lp lig.prmtop \
    -y prod.nc
```

### Workflow 4: Alanine scanning
```bash
cat > mmpbsa_ala.in <<EOF
&general
    verbose=2,
/
&gb
    igb=2, saltcon=0.10,
/
&alanine_scanning
/
EOF

MMPBSA.py -O -i mmpbsa_ala.in \
    -cp com.prmtop -rp rec.prmtop -lp lig.prmtop \
    -mc mutant_complex.prmtop \
    -mr mutant_receptor.prmtop \
    -ml mutant_ligand.prmtop \
    -y prod.nc
```

### Workflow 5: Stability calculation (single trajectory)
```bash
# Only specify complex topology; omit receptor and ligand
MMPBSA.py -O -i mmpbsa.in -o stability.dat \
    -cp com.prmtop \
    -y prod.nc
```

### Workflow 6: Membrane protein MM-PBSA
```bash
cat > mmpbsa_mem.in <<EOF
&general
    use_sander=1,
    startframe=1, endframe=100, interval=1,
    keep_files=0, debug_printlevel=2,
/
&pb
    radiopt=0, indi=20.0, istrng=0.150,
    fillratio=1.25, ipb=1, nfocus=1,
    bcopt=10, eneopt=1, cutfd=7.0, cutnb=99.0,
    npbverb=1, solvopt=2, inp=1,
    memopt=1, emem=7.0, mctrdz=-10.383, mthick=36.086, poretype=1,
    maxarcdot=15000,
/
EOF

MMPBSA.py -O -i mmpbsa_mem.in \
    -cp com.prmtop -rp rec.prmtop -lp lig.prmtop \
    -y prod.nc
```

## Reference Tables

### IGb model comparison
| igb | Model | Description |
|-----|-------|-------------|
| 1 | HCT | Hawkins, Cramer, Truhlar pairwise GB |
| 2 | OBC-I | Onufriev, Bashford, Case; alpha=1.0, beta=0.8, gamma=4.85 |
| 5 | OBC-II | Same as 2 with alpha=1.0, beta=0.8, gamma=4.85 (optimized) |
| 7 | GBn | GB-neck, better for nucleic acids |
| 8 | GBn2 | GB-neck with modified parameters |
| 66 | GBNSR6 | New generalized Born model (default in MMPBSA.py) |

### Decomposition schemes (idecomp)
| idecomp | Type | 1-4 Treatment |
|---------|------|---------------|
| 1 | Per-residue | 1-4 NB added to internal (bond, angle, dihedral) |
| 2 | Per-residue | 1-4 EEL added to EEL, 1-4 VDW added to VDW |
| 3 | Pairwise | 1-4 NB added to internal |
| 4 | Pairwise | 1-4 EEL added to EEL, 1-4 VDW added to VDW |

### Calculation types
| Type | Required topologies | Trajectories |
|------|-------------------|--------------|
| Binding free energy | `-cp -rp -lp` | Complex (receptor/ligand from complex if not given) |
| Stability | `-cp` only | Complex only |
| Alanine scanning | `-cp -rp -lp -mc -mr -ml` | Complex |
| Entropy (nmode) | Same as binding | Same as binding |
| Entropy (quasi-harmonic) | Same as binding | Same as binding |
| QM/MM-GBSA | `-cp -rp -lp` | Complex |
| 3D-RISM | `-cp -rp -lp` + `-xvvfile` | Complex |
| Membrane protein | `-cp -rp -lp` | Complex |

### Output file header contents
```
Ggas = EEL + EVDW + (EBOND + EANGLE + EDIHED)  (gas phase)
Gsolv = EGB/EPB + ESURF/ECAVITY               (solvation)
Delta G = Ggas + Gsolv - T*Delta S             (total)
```

## Worked Example

### MM-GBSA binding free energy of protein-ligand complex
```bash
# 1. Build topologies
tleap -f build.in <<EOF
source leaprc.protein.ff14SB
source leaprc.gaff
source leaprc.water.tip3p
LIG = loadMol2 lig.mol2
loadAmberParams lig.frcmod
receptor = loadPDB protein.pdb
complex = combine {receptor LIG}
set default PBRadii mbondi2
saveAmberParm LIG lig.prmtop lig.crd
saveAmberParm receptor rec.prmtop rec.crd
saveAmberParm complex com.prmtop com.crd
solvateOct complex TIP3PBOX 12.0
addIons complex Na+ 0
saveAmberParm complex com_solvated.prmtop com_solvated.crd
quit
EOF

# 2. Run production MD (see Chapter 9 for details)
pmemd.cuda -O -i prod.in -o prod.out -p com_solvated.prmtop \
    -c equil.rst7 -r prod.rst7 -x prod.nc

# 3. Run MMPBSA.py
cat > mmpbsa.in <<EOF
MM-GBSA binding free energy
&general
    startframe=10, endframe=1000, interval=10,
    verbose=2, keep_files=0,
/
&gb
    igb=5, saltcon=0.150,
/
EOF

MMPBSA.py -O -i mmpbsa.in -o FINAL_RESULTS.dat \
    -sp com_solvated.prmtop \
    -cp com.prmtop -rp rec.prmtop -lp lig.prmtop \
    -y prod.nc
```

## Key Takeaways

1. **Build all prmtop files in the same tleap session** to ensure compatibility (same charges, force field, PBRadii set).
2. **Use `set default PBRadii mbondi2`** before saving prmtops -- this is the recommended radius set for GB models.
3. **The default GB model is igb=66 (GBNSR6)** in MMPBSA.py, not igb=5. Explicitly set `igb=5` if you want the OBC model.
4. **For pairwise decomposition (idecomp=3,4)**, the output scales as O(N^2) -- limit `print_res` to avoid massive files.
5. **Use `-rewrite-output`** to change analysis parameters (e.g., entropy) without rerunning the expensive energy calculations.

## Connects To

- **Chapter 9 (Production MD)**: MMPBSA.py post-processes production MD trajectories
- **Chapter 10 (CPPTRAJ)**: molsurf in cpptraj is used internally for surface area calculations
- **Chapter 11 (Free Energy TI)**: TI is a more rigorous but more expensive alternative to MM-PBSA
- **Chapter 4 (GB/SA)**: The igb models used in MMPBSA.py are the same generalized Born models from sander
- **Chapter 6 (PBSA)**: The PB solver used in MMPBSA.py is the same pbsa program