# Chapter 3: LEaP System Building

## Core Commands & Syntax

### Invocation

```
tleap [-I<dir>] [-f <file>|-]
tleap -f leap_input_file
```

`-I` adds search directories. `-f` reads input from a file (or stdin with `-f -`). Without `-f`, input is interactive.

### Command grammar

```
command argument1 argument2 argument3 ...
variable = command argument1 argument2 ...
```

Commands are case-insensitive (`loadPdb` = `loadpdb`). Variable names are case-sensitive.

### Standard LEaP session structure

```
source leaprc.protein.ff19SB    # load force field
source leaprc.water.opc         # load water model
mol = loadPdb my_protein.pdb    # load structure
check mol                       # verify integrity
addIons mol Na+ 0               # neutralize
solvateOct mol OPCBOX 12.0      # solvate
check mol                       # verify again
saveAmberParm mol prmtop rst7   # save topology + coordinates
savePdb mol system.pdb          # save PDB
quit
```

## Key Namelists / Input Files

### Full LEaP command reference

| Command | Syntax | Purpose |
|---------|--------|---------|
| `source` | `source filename` | Execute LEaP commands from a file |
| `loadPdb` | `var = loadPdb filename` | Load PDB structure into a UNIT |
| `loadPdbUsingSeq` | `loadPdbUsingSeq filename unitlist` | Load PDB with user-defined sequence |
| `sequence` | `var = sequence { UNIT1 UNIT2 ... }` | Link UNITs into a single chain |
| `combine` | `var = combine { UNIT1 UNIT2 ... }` | Combine UNITs without linking |
| `copy` | `newvar = copy variable` | Duplicate an object (deep copy) |
| `createAtom` | `var = createAtom name type charge` | Create a new ATOM |
| `createResidue` | `var = createResidue name` | Create an empty RESIDUE |
| `createUnit` | `var = createUnit name` | Create an empty UNIT |
| `add` | `add container item` | Add ATOM to RESIDUE, or RESIDUE to UNIT |
| `remove` | `remove container item` | Remove ATOM or RESIDUE from container |
| `bond` | `bond atom1 atom2 [order]` | Create bond between atoms (order: `-`, `=`, `#`, `:`) |
| `deleteBond` | `deleteBond atom1 atom2` | Delete bond between atoms |
| `bondByDistance` | `bondByDistance unit [maxBond]` | Auto-bond atoms within distance |
| `set` | `set container property value` | Set properties on objects |
| `desc` | `desc variable` | Print description of an object |
| `list` | `list` | List all defined variables |
| `check` | `check unit [parms]` | Check UNIT for errors (charge, bonds, atom types, close contacts) |
| `edit` | `edit UNIT` | Edit a UNIT interactively |
| `saveAmberParm` | `saveAmberParm unit topfile crdfile` | Write prmtop and inpcrd files |
| `savePdb` | `savePdb unit filename` | Write PDB file |
| `saveOff` | `saveOff object filename` | Save in OFF (Object File Format) library |
| `loadOff` | `loadOff filename` | Load OFF library |
| `loadAmberParams` | `var = loadAmberParams filename` | Load parameter set (frcmod/dat) |
| `loadAmberPrep` | `loadAmberPrep filename [prefix]` | Load AMBER PREP file |
| `loadMol2` | `var = loadMol2 filename` | Load Tripos MOL2 file |
| `saveMol2` | `saveMol2 unit filename type-flag` | Write MOL2 (type-flag: 0=Sybyl, 1=Amber) |
| `addIons` | `addIons unit ion1 numIon1 [ion2 numIon2]` | Add ions using Coulombic grid |
| `addIons2` | `addIons2 unit ion1 numIon1 [ion2 numIon2]` | addIons with solvent+solute treated equally |
| `addIonsRand` | `addIonsRand unit ion1 num1 [ion2 num2] [separation]` | Add ions by replacing random solvent |
| `solvateBox` | `solvateBox solute solvent dist ["iso"] [closeness]` | Add cuboid solvent box |
| `solvateOct` | `solvateOct solute solvent dist ["iso"] [closeness]` | Add truncated octahedron solvent box |
| `solvateCap` | `solvateCap solute solvent pos radius [closeness]` | Add spherical solvent cap |
| `solvateShell` | `solvateShell solute solvent thickness [closeness]` | Add solvent shell |
| `setBox` | `setBox solute enclosure [distance]` | Create periodic box without solvent |
| `alignAxes` | `alignAxes unit` | Align UNIT to principal axes |
| `charge` | `charge unit` | Compute total charge of UNIT |
| `addAtomTypes` | `addAtomTypes { {type element hybrid} ... }` | Define element/hybridization for atom types |
| `addPdbAtomMap` | `addPdbAtomMap { {pdbName leapName} ... }` | Map PDB atom names to LEaP names |
| `addPdbResMap` | `addPdbResMap { {0/1 pdbName leapName} ... }` | Map PDB residue names to LEaP variables |
| `addC4Pairwise` | `addC4Pairwise unit atom1 atom2 C4Value` | Add atom-specific 12-6-4 C4 term |
| `addC4Type` | `addC4Type unit type1 type2 C4Value` | Add type-based 12-6-4 C4 term |
| `impose` | `impose unit seqlist internals` | Set internal coordinates |
| `measureGeom` | `measureGeom a1 a2 [a3 [a4]]` | Measure distance/angle/torsion |
| `transform` | `transform atoms matrix` | Apply symmetry transformation |
| `translate` | `translate atoms direction` | Translate atoms by vector |
| `zMatrix` | `zMatrix object zmatrix` | Define coordinates from internal coords |
| `groupSelectedAtoms` | `groupSelectedAtoms unit name` | Create named atom group |
| `alias` | `alias [string1 [string2]]` | Add/remove/list command aliases |
| `logFile` | `logFile filename` | Open log file |
| `verbosity` | `verbosity level` | Set output verbosity (0=default, 1, 2) |
| `help` | `help [string]` | Show command help |
| `quit` | `quit` | Exit LEaP |

### Object types and their properties

**ATOM properties** (set via `set ATOM property value`):
- `name` -- STRING, atom identifier
- `type` -- STRING, AMBER force field atom type
- `charge` -- NUMBER, electrostatic point charge
- `position` -- LIST of 3 NUMBERs: (X, Y, Z)
- `element` -- STRING, atomic element

**RESIDUE properties**:
- `connect0` -- N-terminal connection atom (head ATOM)
- `connect1` -- C-terminal connection atom (tail ATOM)
- `connect2` -- disulfide bridge connection atom
- `restype` -- "undefined", "solvent", "protein", "nucleic", "saccharide"
- `name` -- STRING, residue name

**UNIT properties**:
- `head` -- connection ATOM for joining to previous UNIT
- `tail` -- connection ATOM for joining to next UNIT
- `box` -- null, NUMBER (cube), or LIST of 3 NUMBERs (orthorhombic)
- `cap` -- LIST of 4 NUMBERs: (X, Y, Z, radius)

**Global defaults** (set via `set default parameter value`):
- `PBRadii` -- "bondi", "mbondi" (default), "mbondi2", "mbondi3", "amber6"
- `Dielectric` -- "distance" (default), "constant"
- `PdbWriteCharges` -- "on" (writes charges to B-factor), "off" (default)
- `nocenter` -- "on" (no centering), "off" (default)
- `reorder_residues` -- "on" (default), "off"
- `OldPrmtopFormat` -- "off" (default)
- `likeCharges` -- check for like charges (default: prompt)
- `keep_chainid` -- "off" (default)
- `FlexibleWater` -- "off" (default), "on"
- `PdbReadBioMT` -- "off" (default)
- `hybrid36` -- "on" (default)

## Common Workflows

### Build and solvate a protein

```
source leaprc.protein.ff19SB
source leaprc.water.opc
mol = loadPdb protein.pdb
check mol
addIons mol Na+ 0
solvateOct mol OPCBOX 12.0
check mol
saveAmberParm mol prmtop rst7
savePdb mol system.pdb
quit
```

### Build a peptide sequence from scratch

```
source leaprc.protein.ff19SB
diala = sequence { ACE ALA NME }
check diala
source leaprc.water.opc
solvateOct diala OPCBOX 10.0
saveAmberParm diala parm7 rst7
quit
```

### Create a disulfide bond

```
mol = loadPdb protein.pdb
bond mol.12.SG mol.35.SG     # CYX residues at positions 12 and 35
```

### Build a branched oligosaccharide (GLYCAM)

```
source leaprc.GLYCAM_06j-1
glycan = sequence { ROH 4YB 4YB 3MB 0MA }
impose glycan {3} { {C1 O4 C4 H4 0.0} }
impose glycan {4} { {C1 O4 C4 H4 0.0} }
saveAmberParm glycan glycan.top glycan.crd
savePdb glycan glycan.pdb
```

Note: For GLYCAM, the sequence order in LEaP is the reverse of the standard notation. Terminal OH is a separate residue. The `impose` command sets glycosidic torsion angles.

### Solvate with a non-water solvent (methanol)

```
source leaprc.protein.ff14SB
loadAmberParams frcmod.meoh
peptide = sequence { ACE VAL NME }
solvateBox peptide MEOHBOX 12.0 0.8
saveAmberParm peptide prmtop prmcrd
quit
```

### Add ions to a specific concentration

```
addIonsRand mol Na+ 50 Cl- 50 5.0   # 50 Na+, 50 Cl-, min 5 Angstrom separation
```

## Reference Tables

### Available pre-equilibrated solvent boxes

| Box name | Solvent | Used with |
|----------|---------|-----------|
| `TIP3PBOX` | TIP3P water | `solvateBox`/`solvateOct` |
| `TIP3PFBOX` | TIP3P/F water | `solvateBox`/`solvateOct` |
| `TIP4PBOX` | TIP4P water | `solvateBox`/`solvateOct` |
| `TIP4PEWBOX` | TIP4P-Ew water | `solvateBox`/`solvateOct` |
| `TIP5PBOX` | TIP5P water | `solvateBox`/`solvateOct` |
| `OPCBOX` | OPC water | `solvateBox`/`solvateOct` |
| `OPC3BOX` | OPC3 water | `solvateBox`/`solvateOct` |
| `OPC3POLBOX` | OPC3-pol water | `solvateBox`/`solvateOct` |
| `POL3BOX` | POL3 water | `solvateBox`/`solvateOct` |
| `SPCBOX` | SPC/E water | `solvateBox`/`solvateOct` |
| `SPCFWBOX` | SPC/Fw water | `solvateBox`/`solvateOct` |
| `QSPCFWBOX` | qSPC/Fw water | `solvateBox`/`solvateOct` |
| `MEOHBOX` | Methanol | `solvateBox`/`solvateOct` |

### solvateBox/solvateOct parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `solute` | UNIT | The UNIT to solvate |
| `solvent` | UNIT | Solvent box UNIT (e.g., `OPCBOX`) |
| `distance` | NUMBER or LIST | Minimum buffer (Angstroms) |
| `"iso"` | STRING (optional) | Make box isometric (rotate solute) |
| `closeness` | NUMBER (optional) | VDW scaling for overlap check (default 1.0) |

For `solvateBox`: if `distance` is a LIST of 3 NUMBERs, they apply to x, y, z axes.
For `solvateOct`: if `distance` is a LIST of 4 NUMBERs, the 4th is diagonal clearance.

### addIons vs addIons2 vs addIonsRand

| Feature | addIons | addIons2 | addIonsRand |
|---------|---------|----------|-------------|
| Method | Coulombic grid | Coulombic grid | Random solvent replacement |
| Solvent treatment | Ignored | Same as solute | Required |
| Ion placement | Energetically favorable | Energetically favorable | Random |
| Speed | Slow | Slow | Fast |
| Reproducibility | Deterministic | Deterministic | Run-dependent |

**Neutralization:** Set numIon = 0 to neutralize (ion1 must be opposite charge to unit).

## Worked Example

Building an alanine dipeptide system with ff19SB and OPC water (from Tutorial 0):

```
# Start tleap
tleap

# Load force field
source leaprc.protein.ff19SB

# Build alanine dipeptide
diala = sequence { ACE ALA NME }

# Load water model and solvate
source leaprc.water.opc
solvateOct diala OPCBOX 10.0

# Check and save
check diala
saveAmberParm diala parm7 rst7

quit
```

This produces `parm7` (topology) and `rst7` (coordinates), ready for minimization and MD with sander or pmemd.

## Key Takeaways

1. **Always `check` before `saveAmberParm`.** The `check` command catches long/short bonds, non-integral charge, missing atom types, and close contacts.
2. **`sequence` links, `combine` does not.** Use `sequence` for proteins/DNA where residues connect head-to-tail. Use `combine` for separate molecules.
3. **`copy` makes a deep copy.** Use it to avoid modifying the original when solvating or adding ions to a variant. Simple assignment (`=` without `copy`) shares the same object.
4. **Solvent order matters.** Solvate after adding ions if using `addIons` (to avoid solvent deletion). Use `addIons2` or `addIonsRand` for solvate-first workflows.
5. **PDB TER records start new UNITs.** They prevent unwanted bonding between chains or across gaps. Use them to separate protein chains, water molecules, and discontinuous segments.

## Connects To

- Chapter 2: Force Fields -- which leaprc files to source
- Chapter 4: PDB Preparation -- preparing PDB files before loading into LEaP
- Section 15: Reading and modifying Amber parameter files
- Section 16: Antechamber and GAFF -- ligand parameterization