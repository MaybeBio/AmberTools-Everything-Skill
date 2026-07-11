# Chapter 18: Membrane Systems

## Core Commands & Syntax

```bash
# Basic membrane-protein packing
packmol-memgen --pdb protein.pdb --lipids POPC --ratio 1

# Binary lipid mixture with ratios
packmol-memgen --pdb protein.pdb --lipids DOPE:DOPG --ratio 3:1

# Multi-component membrane
packmol-memgen --pdb protein.pdb --lipids POPC:POPE:CHL1 --ratio 4:3:3

# Membrane-only system (no protein)
packmol-memgen --lipids POPC --ratio 1 --distxy_fix 100

# Asymmetric bilayer (different lipid composition per leaflet)
packmol-memgen --pdb protein.pdb --lipids DOPC:DOPE//DOPE:DOPS --ratio 3:1//3:1

# Multiple bilayers
packmol-memgen --pdb protein.pdb --lipids DOPC:DOPE --ratio 3:1 \
  --lipids DOPC:DOPS --ratio 3:1

# Add salt and solvent
packmol-memgen --pdb protein.pdb --lipids POPC --ratio 1 \
  --salt --saltconc 0.15 --dist 10 --dist_wat 10

# Auto-parameterize and minimize
packmol-memgen --pdb protein.pdb --lipids POPC --ratio 1 \
  --parametrize --minimize

# Keep intermediate files
packmol-memgen --pdb protein.pdb --lipids POPC --ratio 1 --keep

# Prevent reprotonation
packmol-memgen --pdb protein.pdb --lipids POPC --ratio 1 \
  --notprotonate --nottrim

# OPM/PPM orientation
packmol-memgen --pdb protein.pdb --lipids POPC --ratio 1 --ppm

# Periodic boundary conditions for packing
packmol-memgen --pdb protein.pdb --lipids POPC --ratio 1 --pbc

# Solvate only (no membrane)
packmol-memgen --pdb protein.pdb --solvate --cubic \
  --solute LIG.pdb --solute_con 0.15M

# Mixed solvents
packmol-memgen --pdb protein.pdb --solvate --cubic \
  --solvents WAT:CL3 --solvent_ratio 4:1

# Custom lipid parameters
packmol-memgen --lipids DOPE:TOPC --ratio 3:1 --distxy_fix 100 \
  --keep --parametrize --memgen_parm my_parm.dat

# Curved/buckled membrane
packmol-memgen --xygauss 50 5000 40 --dims 300 50 121 \
  --tight_box --parametrize

# GUI interface
packmol-memgen-gui
packmol-memgen-web

# List available lipids
packmol-memgen --available_lipids
packmol-memgen --available_lipids_all
packmol-memgen --available_solvents

# Help
packmol-memgen --help
packmol-memgen -h
```

## Key Flag Reference

| Flag | Description |
|------|-------------|
| `--pdb` | Input protein/solute PDB file |
| `--lipids` | Colon-separated lipid list (e.g., POPC:POPE) |
| `--ratio` | Colon-separated ratio list (matching `--lipids` order) |
| `--salt` | Add ions to neutralize + achieve salt concentration |
| `--saltconc` | Salt concentration in M (default: 0.15) |
| `--dist` | Water padding distance above/below membrane (Angstrom) |
| `--dist_wat` | Water padding in xy directions (Angstrom) |
| `--distxy_fix` | Fixed xy box dimensions for membrane-only systems |
| `--dims` | Explicit x, y, z box dimensions |
| `--parametrize` | Auto-generate LEaP input and run tleap |
| `--minimize` | Auto-minimize the packed system |
| `--keep` | Keep intermediate files |
| `--notprotonate` | Do not reprotonate protein |
| `--nottrim` | Do not trim hydrogens |
| `--ppm` | Use OPM/PPM for membrane orientation |
| `--pbc` | Enable periodic boundary conditions during packing |
| `--solvate` | Solvate only mode (no membrane) |
| `--cubic` | Use cubic box for solvation |
| `--solvents` | Colon-separated solvent list |
| `--solvent_ratio` | Ratio for solvents in v/v |
| `--solute` | Add additional solute PDB |
| `--solute_con` | Solute concentration (number, M, or %) |
| `--ligand_param` | Lib and frcmod files for custom ligands |
| `--memgen_parm` | Custom lipid parameter file |
| `--premol` | Prebuilt molecule for placement |
| `--xygauss` | Gaussian membrane shape parameters (c, d, h) |
| `--tight_box` | Use tight box dimensions |

## Lipid Naming Convention

```
<sn-1 tail><sn-2 tail><headgroup>
```

### Acyl Chain Abbreviations

| Abbreviation | Full Name | Formula |
|-------------|-----------|---------|
| L | Lauric acid | 12:0 |
| M | Myristoic acid | 14:0 |
| P | Palmitic acid | 16:0 |
| S | Stearic acid | 18:0 |
| O | Oleic acid | 18:1(9) |
| A | Arachidonic acid | 20:4 |
| D/H | Docosahexaenoic acid | 22:6 |

### Head Group Abbreviations

| Abbreviation | Full Name |
|-------------|-----------|
| PC | Phosphatidylcholine |
| PE | Phosphatidylethanolamine |
| PG | Phosphatidylglycerol |
| PA | Phosphatidic acid |
| PS | Phosphatidylserine |
| CL | Cardiolipin |

### Common Lipid Examples

| Lipid | Description |
|-------|-------------|
| POPC | 1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine |
| DOPC | 1,2-dioleoyl-sn-glycero-3-phosphocholine |
| POPE | 1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine |
| DOPE | 1,2-dioleoyl-sn-glycero-3-phosphoethanolamine |
| POPS | 1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoserine |
| POPG | 1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoglycerol |
| CHL1 | Cholesterol |
| TOCL | Tetraoleoyl cardiolipin |

## Common Workflows

### Standard Membrane Protein Setup

```bash
# 1. Pack membrane-protein system
packmol-memgen --pdb protein.pdb --lipids POPC:POPE --ratio 3:1 \
  --salt --saltconc 0.15 --dist 15 --dist_wat 15 --parametrize

# 2. Inspect leap output and generated system
# Outputs: protein_membrane.pdb, leap.in, leap.log
# Protein charmm-gui compatible: protein_g.pdb

# 3. Minimize the system
sander -O -i min.in -o min.out -p system.prmtop -c system.inpcrd \
  -r min.rst -ref system.inpcrd

# 4. Equilibrate with restraints on lipids and protein
# 5. Production MD
```

### Lipid21 Force Field LEaP Setup

```bash
# tleap script for Lipid21 membrane system
cat > leap.in << 'EOF'
source leaprc.protein.ff19SB
source leaprc.lipid21
source leaprc.water.tip3p
loadAmberParams frcmod.ionsjc_tip3p

# Load protein
prot = loadPdb protein_g.pdb

# Load membrane (from packmol-memgen output)
memb = loadPdb bilayer.pdb

# Combine
sys = combine { prot memb }
solvateOct sys TIP3PBOX 15.0
addIonsRand sys Na+ 0
addIonsRand sys Cl- 0

saveAmberParm sys system.prmtop system.inpcrd
quit
EOF

tleap -f leap.in
```

### Membrane Analysis with cpptraj

```bash
# RMSD of protein
cpptraj << EOF
parm system.prmtop
trajin md.nc
rmsd @CA out rmsd_ca.dat
run
EOF

# Lipid order parameters (deuterium order SCD)
cpptraj << EOF
parm system.prmtop
trajin md.nc
lipidorder :POPC out order.dat
run
EOF

# Membrane thickness
cpptraj << EOF
parm system.prmtop
trajin md.nc
membrane thickness out thickness.dat
run
EOF

# Area per lipid
cpptraj << EOF
parm system.prmtop
trajin md.nc
membrane area out area.dat
run
EOF

# Autoimage for periodic boundary handling
cpptraj << EOF
parm system.prmtop
trajin md.nc
autoimage
trajout autoimaged.nc
run
EOF
```

### Custom Lipid Setup

```bash
# 1. Generate custom lipid parameters
packmol-memgen --lipids DOPE:TOPC --ratio 3:1 --distxy_fix 100 \
  --keep --parametrize

# 2. Edit leap.in manually if needed, then run
tleap -f leap.in

# 3. Custom memgen.parm file format:
# RES HEAD_A TAIL_A APL_FF APL VOLUME CHARGE HEAD_PLANE TAIL_PLANE CH3:CH2:CH HEAD_VOLUME CHARMM NAME
```

### SIRAH Coarse-Grained Membrane

```bash
packmol-memgen --pdb protein.pdb --lipids siPOPC --ratio 1 --sirah
```

## Reference Tables

### Lipid Force Field Compatibility

| Force Field | Source | Notes |
|-------------|--------|-------|
| Lipid17 | Amber17 | Previous lipid FF |
| Lipid21 | Amber21 | Current recommended lipid FF |
| Lipid_ext | packmol-memgen | Extended: lysophospholipids, PI, cardiolipins, sterols |
| GLYCAM_06j | Amber | Inositol/phosphate parameters |
| SIRAH | Coarse-grained | si- prefix lipids |

### Extended Lipid Types (Lipid_ext)

| Type | Residue Names | Examples |
|------|--------------|----------|
| Lysophospholipids | PE1, PE2, PG1, PG2 | 2LPC (lyso-PC) |
| Cardiolipins | CLI | TOCL, PODOCL |
| Sterols | ERG, STI, SIT, CAM | Ergosterol, stigmasterol |
| Phosphatidylinositols | PI3, PI4, PI5, P2A-P2F, P3A-P3H | Multiple phosphorylation states |

### memgen.parm File Fields

| Field | Description |
|-------|-------------|
| RES | Residue name and PDB prefix |
| HEAD_A | Atom indices constrained to membrane surface |
| TAIL_A | Atom indices constrained to membrane center |
| APL_FF | Area per lipid in DOPC 4:1 mixture |
| APL | Experimental APL in Angstrom^2 |
| VOLUME | Experimental volume in Angstrom^3 |
| CHARGE | Net lipid charge |
| HEAD_PLANE | Distance to membrane center (~18 A) |
| TAIL_PLANE | Distance to membrane center (~4 A) |
| CH3:CH2:CH | Tail carbon counts |
| HEAD_VOLUME | Headgroup volume in Angstrom^3 |
| CHARMM | Y/N for CHARMM compatibility |
| NAME | Full lipid name |

## Worked Example: Bacterial Membrane Protein

```bash
# 1. Prepare protein PDB
pdb4amber -i 1abc.pdb -o 1abc_clean.pdb --reduce

# 2. Pack into bacterial-like membrane (DOPE:DOPG 3:1)
packmol-memgen --pdb 1abc_clean.pdb --lipids DOPE:DOPG --ratio 3:1 \
  --salt --saltconc 0.15 --dist 15 --dist_wat 15 \
  --parametrize --keep

# 3. Minimize
cat > min.in << 'EOF'
&cntrl
  imin=1, maxcyc=5000, ncyc=2500,
  ntb=1, ntc=1, ntf=1,
  cut=10.0,
  ntr=1, restraint_wt=10.0, restraintmask='!@H=',
/
EOF

sander -O -i min.in -o min.out -p system.prmtop -c system.inpcrd \
  -r min.rst -ref system.inpcrd

# 4. Equilibrate (NVT then NPT with lipid restraints)
# 5. Production NPT
# 6. Analyze
cpptraj << EOF
parm system.prmtop
trajin prod.nc
autoimage
rmsd @CA out rmsd.dat mass
lipidorder :DOPE,DOPG out scd.dat
membrane thickness out thickness.dat
run
EOF
```

## Key Takeaways

1. **packmol-memgen is the primary tool** for building Amber-ready membrane systems, using Memembed for orientation, pdbremix for volume estimation, and Packmol as the packing engine.
2. **Lipid naming**: `<sn-1 tail><sn-2 tail><headgroup>`. Asymmetric bilayers use `//` to separate leaflets. Multiple bilayers use repeated `--lipids` flags.
3. **Lipid21 is the current recommended lipid force field**, loaded via `source leaprc.lipid21` in tleap. Lipid_ext extends this with lysophospholipids, cardiolipins, and additional sterols.
4. **Analysis with cpptraj**: Use `lipidorder` for SCD order parameters, `membrane thickness` for bilayer thickness, and `membrane area` for area per lipid.
5. **Always minimize after packing** due to potential initial clashes. Use `--minimize` flag for auto-minimization, or run sander/pmemd manually with moderate restraints.

## Connects To

- **Chapter 5: System Setup with LEaP** - tleap integration for membrane system building
- **Chapter 2: Force Fields** - Lipid21, force field selection
- **Chapter 12: Running MD** - Equilibration protocols for membrane systems
- **Chapter 10: Analysis with cpptraj** - Membrane-specific analysis
- **Chapter 20: Implicit Solvent** - PBSA implicit membrane models
- AMBER tutorials: `05_building_systems/01-9_Building Membrane Systems-Overview.md`