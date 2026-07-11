# Chapter 17: Metal Ion Modeling with pyMSMT

## Core Commands & Syntax

```bash
# MCPB.py - main metal center parameter builder
MCPB.py -i input.in -s 1       # Step 1: generate PDB/fingerprint/Gaussian input files
MCPB.py -i input.in -s 2       # Step 2: generate frcmod file (bond/angle parameters)
MCPB.py -i input.in -s 3       # Step 3: charge fitting, generate mol2 files
MCPB.py -i input.in -s 4       # Step 4: generate LEaP input for bonded model
MCPB.py -i input.in -s 4n1     # Step 4n1: nonbonded model with refitted charges
MCPB.py -i input.in -s 4n2     # Step 4n2: nonbonded model, no charge refitting
MCPB.py -i input.in -s 2s      # Step 2s: Seminario method (default)
MCPB.py -i input.in -s 2z      # Step 2z: Z-matrix method
MCPB.py -i input.in -s 2e      # Step 2e: Empirical method (Zn2+ only)
MCPB.py -i input.in -s 2b      # Step 2b: blank (zero force constants)

# Log file and fchk file overrides
MCPB.py -i input.in -s 2 --logf my_sidechain_fc.log --fchk my_sidechain_opt.fchk

# IPMach.py - ion parameterization machine
IPMach.py -i inputfile

# Other pyMSMT tools
metalpdb2mol2.py input.pdb     # Convert PDB to mol2 for metal ions
CartHess2FC.py                 # Force constants from Cartesian Hessian
OptC4.py                       # Optimize C4 terms for 12-6-4 potential
```

## Key Input File: MCPB.py (`input.in`)

MCPB.py model names: `sidechain_model` (ligating residues only), `standard_model` (sidechain + backbone), `large_model` (large QM region).

### Required Variables

```
original_pdb 4zaf_H.pdb           # PDB file name (single chain, with H atoms)
group_name MOL                     # Prefix for modeling files (default: MOL)
ion_ids 2401                       # PDB atom ID(s) of central metal ion(s)
cut_off 2.8                        # Metal-ligand bond cutoff (Angstroms, default: 2.8)
ion_mol2files ZN.mol2              # Pre-built mol2 files for metal ions
ion_info ZN ZN Zn 2                # Only for step 4n2: resName atomName element charge
```

### Optional Variables

```
force_field ff19SB                 # ff94/ff99/ff99SB/ff03/ff10/ff14SB/ff14SBonlysc/ff19SB/ff15ipq
gaff 1                             # 0=no, 1=GAFF, 2=GAFF2
software_version gau               # g03/g09/g16/gau/gms
large_opt 1                        # 0=no opt, 1=H-only opt, 2=full opt
water_model OPC                    # TIP3P/SPCE/TIP4PEW/OPC3/OPC/FB3/FB4
naa_mol2files LIG.mol2            # Non-amino-acid mol2 files (ligands, hydroxyl)
frcmod_files LIG.frcmod           # Parameter modification files for nonstandard residues
scale_factor 0.81                  # Force constant scale factor (default: 1.0)
ion_paraset 12_6_4                 # HFE/CM/IOD/12_6/12_6_4 (default: 12_6)
add_bonded_pairs 1001-1320         # Additional Metal-C bonds (not auto-detected)
additional_resids 100 101           # Extra residues for models
chgfix_resids 100                  # Residues with fixed charges during fitting
smmodel_chg 2                      # Small model charge (auto-detected, override if wrong)
smmodel_spin 1                     # Small model spin multiplicity (default: 1 or 2)
lgmodel_chg 0                      # Large model charge
lgmodel_spin 1                     # Large model spin
smmodel_frz_ids 1320 1380          # Atom IDs to freeze in small model opt
sqm_opt 2                          # 0=none, 1=sidechain, 2=large, 3=both
xstru 1                            # 0=QM structure, 1=original PDB structure
anglefc_avg 0                      # Average angle force constants (0 or 1)
bondfc_avg 0                       # Average bond force constants (0 or 1)
add_redcrd 2                       # Redundant coordinates for Z-matrix (0/1/2)
```

## Common Workflows

### Bonded Model Parameterization (4 Steps)

```bash
# Step 1: Generate models (1a = auto-rename atom types, default)
MCPB.py -i input.in -s 1a

# --- Run Gaussian calculations manually ---
# - Sidechain model: geometry optimization + frequency (force constants)
# - Large model: RESP charge calculation
# (Use Opt=CalcAll keyword, default since 2026)

# --- Override method/basis in Gaussian input files as needed ---
# method: B3LYP, basis_set: 6-31G*, chg: -2, spin: 1, multiplicity: 1

# Step 2: Generate frcmod (Seminario method default)
MCPB.py -i input.in -s 2s --fchk MOL_sidechain_opt.fchk

# Step 3: Charge fitting (3b = restrain backbone heavy atoms, default)
MCPB.py -i input.in -s 3b --logf MOL_large_mk.log

# Step 4: Generate LEaP input for bonded model
MCPB.py -i input.in -s 4b
```

### Nonbonded Model Parameterization

```bash
# Option A: With refitted charges (3 steps)
MCPB.py -i input.in -s 1m       # Rename metal ion atom types only
MCPB.py -i input.in -s 3b       # RESP charge fitting
MCPB.py -i input.in -s 4n1      # Nonbonded model with refitted charges

# Option B: Without charge fitting (1 step)
MCPB.py -i input.in -s 4n2      # Nonbonded model, default charges
```

### 12-6-4 LJ Nonbonded Model Setup

The 12-6-4 model adds a C4/r^4 term to the standard 12-6 LJ potential for improved metal ion coordination. Use `ion_paraset 12_6_4` in the MCPB input file. For zinc parameters specifically, the 12-6-4 model captures ion-induced dipole interactions critical for Zn2+ sites.

```bash
# Optimize C4 terms to reproduce experimental structure
OptC4.py -i optimized_system.prmtop -c optimized_system.inpcrd
```

### Parameter File Outputs

After MCPB.py, the key output files are:
- `MOL_mcpbpy.frcmod` - bonded parameters (bond, angle, dihedral)
- `MOL_mcpbpy.mol2` - residue mol2 files with fitted charges
- `MOL_mcpbpy.in` - LEaP input file for system building
- `MOL_sidechain_fc.log` - QM force constant calculation log
- `MOL_large_mk.log` - QM RESP charge calculation log

### LEaP Integration

```bash
tleap -f MOL_mcpbpy.in
# Loads force fields, frcmod files, mol2 files, builds system
# Outputs: prmtop and inpcrd for MD simulation
```

## Reference Tables

### MCPB Step Options

| Step | Option | Description |
|------|--------|-------------|
| 1 | 1a (default) | Auto-rename metal ion + ligating atom types |
| 1 | 1m | Rename metal ion atom types only |
| 1 | 1n | No renaming, manual editing |
| 2 | 2s (default) | Seminario method (requires fchk file) |
| 2 | 2z | Z-matrix method (requires log file) |
| 2 | 2e | Empirical method (Zn2+ only) |
| 2 | 2b | Blank method (zero force constants) |
| 2 | 2ms | Modified Seminario method |
| 3 | 3b (default) | RESP: restrain backbone heavy atoms |
| 3 | 3a | RESP: no restrictions |
| 3 | 3c | RESP: restrain backbone (all atoms) |
| 3 | 3d | RESP: restrain backbone + C-beta |
| 3 | 3f | Fluctuating charge model (3d metals only) |
| 4 | 4b (default) | Bonded model |
| 4 | 4n1 | Nonbonded model with charge refitting |
| 4 | 4n2 | Nonbonded model without charge refitting |

### Supported QM Software

| Option | Software | Format |
|--------|----------|--------|
| gau/g03/g09/g16 | Gaussian | Log for ESP, fchk for Hessian |
| gms | GAMESS-US | Log for both ESP and Hessian |

### Nonbonded Model Parameter Sets

| Set | Description | Recommended For |
|-----|-------------|-----------------|
| HFE | Hydration free energy optimized | +1 and -1 ions |
| CM | Cation parameters (Merz group) | +2 ions |
| IOD | Ion-oxygen distance optimized | +3 and +4 ions |
| 12_6 | Equivalent to CM set | Default |
| 12_6_4 | 12-6-4 LJ-type (C4 term) | Polarizable metal ions |

### Zinc LJ Parameters (12-6-4 model)

| Parameter | Description |
|-----------|-------------|
| radius | VDW radius (Angstrom) |
| epsilon | Well depth (kcal/mol) |
| rmin | Rmin = rmin/2 * 2 |
| sigma | sigma = rmin/2 / 2^(1/6) |
| C4 | C4 term for ion-induced dipole |
| pol | Polarizability |
| alpha | Polarizability (alternate) |

## Worked Example: Zn2+ Metalloprotein

```bash
# 1. Clean PDB with pdb4amber, add hydrogens with H++ server
pdb4amber -i 4zaf.pdb -o 4zaf_H.pdb --reduce

# 2. Prepare metal ion mol2 file
# Use antechamber to create ZN.mol2, manually set atom type and charge

# 3. Create MCPB input file (mcppb.in):
# original_pdb 4zaf_H.pdb
# group_name ZN1
# ion_ids 2401
# cut_off 2.8
# ion_mol2files ZN.mol2
# force_field ff19SB
# software_version g16
# large_opt 1
# water_model OPC

# 4. Run MCPB steps
MCPB.py -i mcppb.in -s 1a

# 5. Submit Gaussian calculations (edit *_sidechain_opt.gjf and *_large_mk.gjf)
# Sidechain: #p B3LYP/6-31G* Opt Freq
# Large: #p B3LYP/6-31G* Pop=MK IOp(6/33=2)

# 6. After Gaussian completes:
MCPB.py -i mcppb.in -s 2s --fchk ZN1_sidechain_opt.fchk
MCPB.py -i mcppb.in -s 3b --logf ZN1_large_mk.log
MCPB.py -i mcppb.in -s 4b

# 7. Build system in LEaP
tleap -f ZN1_mcpbpy.in

# 8. Run MD simulation
sander -O -i md.in -o md.out -p ZN1_mcpbpy.prmtop -c ZN1_mcpbpy.inpcrd \
  -r md.rst -x md.nc
```

## Key Takeaways

1. **Four-step workflow**: 1a (generate models) -> 2s (Seminario) -> 3b (RESP charge fitting) -> 4b (LEaP input) for bonded model parameterization.
2. **Seminario is default**: Uses Cartesian Hessian sub-matrices from QM frequency calculations; requires Gaussian fchk file. Z-matrix is alternative.
3. **Between steps 1 and 2**: Run QM calculations manually (sidechain opt+freq, large RESP). You can customize method, basis set, charge, spin, and multiplicity in the Gaussian input files.
4. **Nonbonded model**: One-step with `-s 4n2`; requires `ion_info` variable. For 12-6-4 LJ, set `ion_paraset 12_6_4`.
5. **Supports 80+ metal ions**: Charges/oxidation states +1 to +8, across multiple AMBER force fields (ff94 through ff19SB, GAFF, GAFF2).

## Connects To

- **Chapter 5: System Setup with LEaP** - MCPB.py output integrates with tleap
- **Chapter 6: Antechamber** - Generate mol2 files for ligands and nonstandard residues
- **Chapter 3: Force Fields** - ff19SB, ff14SB, ff14SBonlysc, GAFF/GAFF2
- **Chapter 12: Running MD** - sander/pmemd for production simulations
- **Chapter 19: QM/MM** - QM/MM as alternative to MM metal ion models
- AMBER tutorials: `06_Developing_Nonstandard_Parameters/02-4_Metal Ion Modeling Tutorial.md`