# Chapter 2: Molecular Mechanics Force Fields

## Core Commands & Syntax

### Loading force fields in LEaP

The recommended preamble for loading force fields (order matters):

```
source leaprc.protein.ff19SB
source leaprc.DNA.OL24
source leaprc.lipid21
source leaprc.water.opc
source leaprc.gaff2
```

**Critical ordering rule:** The water model must be loaded after the protein force field. Explicit solvent simulations now require an explicit `leaprc.water.xxxx` -- TIP3P is no longer loaded by default.

### Loading user-defined parameter modifications

```
loadamberparams user-defined-file.frcmod
```

### Force field file locations

Standard leaprc files: `$AMBERHOME/dat/leap/cmd/`
Older (deprecated) force fields: `$AMBERHOME/dat/leap/cmd/oldff/`

## Key Namelists / Input Files

### Force field equation (Amber functional form)

```
E_total = sum_bonds Kb*(b - b0)^2
        + sum_angles Ktheta*(theta - theta0)^2
        + sum_dihedrals (Vn/2)*[1 + cos(n*phi - gamma)]
        + sum_nonbond(i<j) [A_ij/R_ij^12 - B_ij/R_ij^6 + q_i*q_j/(epsilon*R_ij)]
```

- **Bond term:** Harmonic spring with force constant `Kb`, equilibrium length `b0`
- **Angle term:** Harmonic spring with force constant `Ktheta`, equilibrium angle `theta0`
- **Dihedral term:** Periodic with barrier height `Vn`, periodicity `n`, phase `gamma`
- **VDW term:** Lennard-Jones 12-6 potential with `A_ij` and `B_ij` parameters
- **Electrostatic term:** Coulomb potential with partial charges `qi`, `qj`

Amber does not have a separate improper dihedral term -- out-of-plane motions use the same formulation as torsion terms. No Urey-Bradley term (unlike CHARMM).

### frcmod file format

A `frcmod` file contains modifications or additions to force field parameters:

```
MASS
atom_type   atomic_mass

BOND
atom_type1-atom_type2   Kb   b0

ANGLE
atom_type1-atom_type2-atom_type3   Ktheta   theta0

DIHEDRAL
atom_type1-atom_type2-atom_type3-atom_type4   IDIVF   Vn   gamma   n

IMPROPER
atom_type1-atom_type2-atom_type3-atom_type4   Vn   gamma   n

NONBOND
atom_type   Rmin   epsilon
```

### Water model leaprc files

| Water model | LEaP command | Residue name | Points | Charges |
|-------------|-------------|--------------|--------|---------|
| TIP3P | `source leaprc.water.tip3p` | TP3 | 3 | 3 |
| TIP4P | `source leaprc.water.tip4p` | TP4 | 4 | 3 |
| TIP4P-Ew | `source leaprc.water.tip4pew` | T4E | 4 | 3 |
| TIP5P | `source leaprc.water.tip5p` | TP5 | 5 | 3 |
| OPC | `source leaprc.water.opc` | OPC | 4 | 3 |
| OPC3 | `source leaprc.water.opc3` | OP3 | 3 | 3 |
| OPC3-pol | `source leaprc.water.opc3pol` | O3P | 3 | 3 |
| POL3 | `source leaprc.water.pol3` | PL3 | 3 | 3 |
| SPC/E | `source leaprc.water.spce` | SPC | 3 | 3 |
| SPC/Eb | `source leaprc.water.spceb` | SPC | 3 | 3 |
| TIP3PFB | `source leaprc.water.tip3pfb` | FB3 | 3 | 3 |
| TIP4PFB | `source leaprc.water.tip4pfb` | FB4 | 4 | 3 |

To override the default residue name for water in a PDB file:

```
WAT = OPC                          # residues named WAT become OPC
source leaprc.water.opc
```

## Common Workflows

### Standard protein + water + ions setup

```
source leaprc.protein.ff19SB
source leaprc.water.opc
prot = loadPdb protein.pdb
check prot
addIons prot Na+ 0                # neutralize with Na+
solvateOct prot OPCBOX 12.0
saveAmberParm prot prmtop rst7
quit
```

### Protein + DNA + lipid combined system

```
source leaprc.protein.ff19SB
source leaprc.DNA.OL24
source leaprc.lipid21
source leaprc.water.opc
source leaprc.gaff2
```

### Using a flexible water model

```
WAT = SPG
loadAmberParams frcmod.qspcfw
set default FlexibleWater on
```

### Building a simple peptide in methanol

```
source leaprc.protein.ff14SB
loadAmberParams frcmod.meoh
peptide = sequence { ACE VAL NME }
solvateBox peptide MEOHBOX 12.0 0.8
saveAmberParm peptide prmtop prmcrd
quit
```

## Reference Tables

### Protein force fields

| Force field | leaprc | Highlights | Recommended water |
|-------------|--------|------------|-------------------|
| ff19SB | `leaprc.protein.ff19SB` | Amino-acid specific CMAP, new XC atom type | OPC (strongly) |
| ff14SB | `leaprc.protein.ff14SB` | Standard SB backbone, CX atom type | TIP3P or OPC |
| ff14SBonly | `leaprc.protein.ff14SBonly` | ff14SB + frcmod.ff99SB14 | -- |
| ff15ipq | `leaprc.protein.ff15ipq` | IPolQ charges, ~1200 parameters | SPC/E-b |
| fb15 | `leaprc.protein.fb15` | Force balance model | TIP3P-FB |
| ff03.r1 | `leaprc.protein.ff03.r1` | ff03 variant | -- |

### Nucleic acid force fields

| Force field | leaprc | Description |
|-------------|--------|-------------|
| OL24 | `leaprc.DNA.OL24` | Current DNA recommended |
| OL21 | `leaprc.DNA.OL21` | Previous DNA (kept for backup) |
| OL3 | `leaprc.RNA.OL3` | RNA (Chi-OL3) |
| OL15 | `leaprc.DNA.OL15` | Older DNA |

### Lipid force fields

| Force field | leaprc | Key files | Residues |
|-------------|--------|-----------|----------|
| lipid21 | `leaprc.lipid21` | `lipid21.lib`, `lipid21.dat` | 8 tail groups, 6 head groups, cholesterol |

LIPID21 tail residues: LAL (12:0), MY (14:0), PA (16:0), SA (16:1), OL (18:1 n-9), ST (18:0), AR (20:4), DHA (22:6)
LIPID21 head residues: PC, PE, PS, PGR, PGS, PH-, SPM
Other: CHL (cholesterol)

LIPID11, LIPID14, and LIPID17 are deprecated to `oldff/`.

### Water model comparison (bulk properties)

| Model | Dielectric constant | Notes |
|-------|--------------------|-------|
| OPC | 78.4 (exp: 78.4) | 0.76% avg error across 11 properties |
| TIP3P | 94 | Widely used, seriously limited with ff19SB |
| TIP4P-Ew | 63.9 | Common 4-point alternative |
| OPC3 | -- | 3-point OPC variant |

### Ion parameter sets (with OPC)

| Set | Purpose |
|-----|---------|
| 12-6 normal usage | General simulations (HFE monovalent + CM polyvalent) |
| 12-6 HFE | Reproduce experimental hydration free energy |
| 12-6 IOD | Reproduce experimental ion-oxygen distance |
| 12-6-4 | LJ-type with C4 term; best transferability |

## Worked Example

Building a protein-DNA-lipid system with ff19SB, OPC water, and neutralizing ions:

```
# Load force fields
source leaprc.protein.ff19SB
source leaprc.DNA.OL24
source leaprc.lipid21
source leaprc.water.opc

# Load structures
prot = loadPdb protein_cleaned.pdb
dna  = loadPdb dna_cleaned.pdb

# Combine
sys = combine { prot dna }
check sys

# Neutralize and solvate
addIons sys Na+ 0
solvateOct sys OPCBOX 12.0

# Verify
check sys

# Save
saveAmberParm sys system.prmtop system.rst7
savePdb sys system.pdb
quit
```

## Key Takeaways

1. **Order matters when loading leaprc files.** Load water last. Load protein before nucleic acids.
2. **ff19SB + OPC is the current recommended combination.** ff19SB uses amino-acid-specific CMAP corrections with a new XC atom type for C-alpha. Do NOT use ff19SB with TIP3P.
3. **Explicit water model loading is required.** The `leaprc.water.xxxx` must be explicitly sourced -- TIP3P is no longer the default.
4. **Mix-and-match is risky.** The recommended force field combinations have been extensively tested together. Ad-hoc combinations require deeper knowledge.
5. **The 12-6-4 LJ-type nonbonded model is preferred for multivalent ions.** It adds a C4 term to the standard 12-6 potential and shows excellent transferability.

## Connects To

- Chapter 3: LEaP System Building -- how to use these force fields in practice
- Chapter 4: PDB Preparation -- preparing PDB files before loading into LEaP
- Section 15: Reading and modifying Amber parameter files
- Section 16: Antechamber and GAFF -- small molecule force field parameterization