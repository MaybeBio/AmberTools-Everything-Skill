# Chapter 4: PDB Preparation

## Core Commands & Syntax

### pdb4amber

```
pdb4amber [-h] [-i FILE] [-o FILE] [-y] [-d] [-s STRIP_ATOM_MASK]
          [-m MUTATION_STRING] [-p] [--constantph] [--most-populous]
          [--keep-altlocs] [--reduce] [--no-reduce-db] [--pdbid]
          [--add-missing-atoms] [--model MODEL] [-l FILE] [-v]
          [--leap-template] [--no-conect] [--noter]
          [input]
```

Basic usage:

```
pdb4amber -i input.pdb -o output.pdb
pdb4amber -i input.pdb -o output.pdb --reduce --add-missing-atoms
pdb4amber -i input.pdb -o output.pdb -d -y            # dry + no hydrogens
pdb4amber -i input.pdb -o output.pdb -p               # keep only Amber-compatible residues
```

### reduce

The `reduce` program adds hydrogens to PDB files. It can be invoked standalone or via pdb4amber:

```
pdb4amber -i input.pdb -o output.pdb --reduce
```

The `--no-reduce-db` flag skips running reduce on heteroatoms.

### pdb2pqr

```
pdb2pqr --ff=amber input.pdb output.pqr
```

### gwh (guess water hydrogens)

Used to add hydrogens to crystal water molecules. Invoked via LEaP or standalone.

## Key Namelists / Input Files

### pdb4amber options

| Flag | Description |
|------|-------------|
| `-i FILE, --in FILE` | Input PDB file (default: stdin) |
| `-o FILE, --out FILE` | Output PDB file (default: stdout) |
| `-y, --nohyd` | Remove all hydrogen atoms |
| `-d, --dry` | Remove all water molecules (WAT, HOH) |
| `-s MASK, --strip MASK` | Strip atoms matching Amber mask |
| `-m STR, --mutate STR` | Mutate a residue (e.g., `ALA-42-GLY`) |
| `-p, --prot` | Keep only Amber-compatible residues |
| `--constantph` | Rename GLU, ASP, HIS for constant pH |
| `--most-populous` | Keep most populous alternative conformation |
| `--keep-altlocs` | Keep all alternative conformations |
| `--reduce` | Run reduce to add hydrogens first |
| `--no-reduce-db` | Skip reduce for heteroatoms |
| `--pdbid` | Fetch structure by PDB ID |
| `--add-missing-atoms` | Use tleap to add missing atoms |
| `--model MODEL` | Model number from multi-model PDB (negative = keep all) |
| `-l FILE, --logfile FILE` | Log file |
| `--leap-template` | Write LEaP template for easy adaptation |
| `--no-conect` | Do not write S-S CONECT records |
| `--noter` | Do not write TER records |
| `-v, --version` | Show version |
| `-h, --help` | Show help |

### Automatic output files from pdb4amber

Running pdb4amber creates:
- The cleaned PDB file (specified by `-o`)
- `outputname_renum.txt` -- residue renumbering table (maps original to 1-based sequential numbering)
- `outputname_nonprot.pdb` -- if `-p` is used, contains non-protein atoms
- `outputname_sslink` -- detected disulfide bonds

### PDB format requirements for Amber

1. **Heavy atoms only** for amino acids (LEaP adds hydrogens). Delete all explicit hydrogens -- they may have names LEaP does not recognize.
2. **ATOM records only** -- remove HETATM records for standard residues.
3. **TER records** between chains, across gaps, and between individual water molecules.
4. **Unique residue names** reflecting protonation states (see naming conventions below).
5. **No connectivity records** (CONECT) -- disulfide bonds are made explicitly with the `bond` command in LEaP.
6. **Delete non-essential components** -- ligands, ions, co-factors, or save them separately.
7. **Single conformation** -- choose one alternate location. The default is "A" conformation.

## Common Workflows

### Standard PDB cleanup pipeline

```bash
# 1. Basic cleanup: remove hydrogens, water, keep only protein
pdb4amber -i raw.pdb -o clean.pdb -y -d -p

# 2. Add hydrogens with reduce, then add missing atoms
pdb4amber -i clean.pdb -o final.pdb --reduce --add-missing-atoms

# 3. Inspect the renumbering table
cat final_renum.txt
```

### Handling ligands

```bash
# 1. Save ligand separately (in a viewer or text editor)
# 2. Remove ligand from PDB file
# 3. Clean the protein PDB
pdb4amber -i protein_only.pdb -o protein_clean.pdb --reduce

# 4. Parameterize ligand separately with antechamber
antechamber -i ligand.sdf -fi sdf -o ligand.mol2 -fo mol2 -c bcc -s 2
parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod
```

### Handling disulfide bonds

```bash
# 1. In your PDB file, rename cysteines in S-S bridges to CYX
#    (CYS -> CYX for each cysteine in a disulfide)

# 2. Clean with pdb4amber
pdb4amber -i protein.pdb -o protein_clean.pdb --reduce

# 3. In LEaP, create the bond explicitly:
#    bond mol.12.SG mol.35.SG
```

### Multi-model PDB files

```bash
pdb4amber -i nmr_ensemble.pdb -o model1.pdb --model 1
pdb4amber -i nmr_ensemble.pdb -o all_models.pdb --model -1
```

### Constant pH preparation

```bash
pdb4amber -i input.pdb -o output.pdb --constantph
# Renames: GLU->GL4, ASP->AS4, HIS->HIP (adds alternative protonation states)
```

### Fetching a structure by PDB ID

```bash
pdb4amber --pdbid -i 1CRN -o 1crn_clean.pdb --reduce
```

### Handling crystal waters

```bash
# 1. Decide which crystal waters are important (e.g., ligand-binding)
# 2. Keep them in the PDB as WAT or HOH residues
# 3. Delete the rest
# 4. Separate each water with a TER record
# 5. In LEaP, crystal waters are recognized by residue name WAT or HOH
```

## Reference Tables

### Residue naming conventions for protonation states

| Standard name | Amber name | State |
|---------------|------------|-------|
| HIS | HIE | epsilon-protonated (default) |
| HIS | HID | delta-protonated |
| HIS | HIP | doubly protonated (positive) |
| CYS | CYS | Free cysteine (-SH) |
| CYS | CYX | Disulfide-bonded cysteine (no hydrogen on SG) |
| ASP | ASP | Deprotonated (negative, default) |
| ASP | ASH | Protonated (neutral) |
| GLU | GLU | Deprotonated (negative, default) |
| GLU | GLH | Protonated (neutral) |
| LYS | LYS | Protonated (positive, default) |
| LYS | LYN | Deprotonated (neutral) |

### Terminal cap residues

| Cap | PDB residue name | Position | Required atoms |
|-----|-----------------|----------|----------------|
| Acetyl | ACE | N-terminal | CH3, C, O |
| NH2 | NHE | C-terminal | N |
| N-methylamide | NME | C-terminal | N, C |

**ACE PDB format:**
```
ATOM    1  CH3 ACE     1       x       y       z
ATOM    2  C   ACE     1       x       y       z
ATOM    3  O   ACE     1       x       y       z
```

**NHE PDB format:**
```
ATOM    1  N   NHE     1       x       y       z
```

**NME PDB format:**
```
ATOM    1  N   NME     1       x       y       z
ATOM    2  C   NME     1       x       y       z
```

Hydrogens are added automatically. Do not include them.

### PDB cleanup checklist

| Step | Action | Tool |
|------|--------|------|
| 1 | Visual inspection | PyMOL, VMD, ChimeraX |
| 2 | Delete non-protein/non-water components | Text editor |
| 3 | Choose single unit cell image | Text editor |
| 4 | Extract ligand to separate file | Viewer or text editor |
| 5 | Delete irrelevant water molecules | Text editor |
| 6 | Remove explicit hydrogens | `pdb4amber -y` or manually |
| 7 | Set protonation state names | Text editor (HIE/HID/HIP, CYX, ASH, etc.) |
| 8 | Add TER records between chains/gaps | Text editor |
| 9 | Remove CONECT records | Text editor |
| 10 | Run pdb4amber | `pdb4amber -i in.pdb -o out.pdb --reduce` |
| 11 | Verify renumbering | Check `_renum.txt` file |

### PDB numbering in LEaP

LEaP renumbers residues sequentially from 1, ignoring original PDB numbers and gaps. Example:

```
Original: 21 22 23 24 25 ... 31 32 33 ... 80
            (gap: residues 26-30 missing)
LEaP:      1  2  3  4  5 ...  6  7  8 ... 55
```

**Always use the LEaP numbering** for masks, restraints, and selection in subsequent AMBER simulations. Use the `_renum.txt` file from pdb4amber to map between original and LEaP numbering.

## Worked Example

Preparing crambin (1CRN) for MD simulation:

```bash
# 1. Fetch and clean the PDB
pdb4amber --pdbid -i 1CRN -o 1crn_clean.pdb --reduce -p

# 2. Inspect: check for disulfide bonds
# 1CRN has disulfides: 3-40, 4-32, 16-26
# Rename those cysteines CYS->CYX in the PDB file

# 3. In LEaP:
# source leaprc.protein.ff19SB
# source leaprc.water.opc
# mol = loadPdb 1crn_clean.pdb
# bond mol.3.SG mol.40.SG
# bond mol.4.SG mol.32.SG
# bond mol.16.SG mol.26.SG
# check mol
# addIons mol Na+ 0
# solvateOct mol OPCBOX 12.0
# saveAmberParm mol prmtop rst7
# quit
```

## Key Takeaways

1. **Always clean your PDB before LEaP.** This is the most common source of errors. Remove hydrogens, non-standard residues, and connectivity records. Set proper residue names for protonation states.
2. **TER records are critical.** They prevent LEaP from creating spurious bonds across chain breaks, gaps, and between water molecules. Each chain, gap, and individual water should be separated by TER.
3. **Residue numbering changes.** LEaP renumbers from 1 sequentially. Use the `_renum.txt` file from pdb4amber and always reference the final LEaP numbering in masks and restraints.
4. **Disulfide bonds are explicit.** Rename CYS to CYX in the PDB and use the `bond` command in LEaP. Do not rely on CONECT records.
5. **Ligands must be prepared separately.** Remove them from the protein PDB, keep their heavy-atom coordinates, and parameterize with antechamber/GAFF.

## Connects To

- Chapter 3: LEaP System Building -- loading cleaned PDB files into LEaP
- Chapter 2: Force Fields -- protonation states and residue naming
- Section 16: Antechamber and GAFF -- ligand parameterization
- Section 37.1: ambpdb -- converting Amber output back to PDB format