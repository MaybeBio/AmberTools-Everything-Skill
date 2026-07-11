# Chapter 6: ParmEd -- Topology File Editor

## Core Commands & Syntax

### Launching parmed
```bash
parmed [-h] [-v] [-i <script>] [-p <prmtop>] [-c <inpcrd>] [-O]
       [-l <logfile>] [--prompt <prompt>] [-n] [-e] [-s] [-r]
       [<prmtop>] [<script>]
```

| Flag | Description |
|------|-------------|
| `-p <prmtop>` | Topology file(s) to load (repeatable) |
| `-c <inpcrd>` | Coordinate file(s) to load (paired with `-p` order) |
| `-i <script>` | Script file with ParmEd commands |
| `-O` | Allow overwriting existing files |
| `-l <logfile>` | Log file (default: `parmed.log`) |
| `-e` | Enable Python interpreter commands |
| `-s`, `--strict` | Halt on unrecognized input (default) |
| `-r`, `--relaxed` | Skip over failed actions |
| `-n`, `--no-splash` | Suppress greeting logo |

### Interactive session
```bash
parmed system.prmtop
> printDetails :1-10
> changeRadii mbondi3
> outparm system_mbondi3.parm7
> quit
```

### Running a script non-interactively
```bash
parmed -p system.prmtop -i commands.in
```

## Key Namelists / Input Files

### Action Command Reference (all case-insensitive)

#### Inspection Commands
| Command | Usage | Description |
|---------|-------|-------------|
| `printDetails` | `printDetails <mask>` | Print atom number, residue, atom name, type, VDW radius, well depth, mass, charge |
| `printFlags` | `printFlags` | List all %FLAG sections in the prmtop |
| `printBonds` | `printBonds <mask> [<mask>]` | Print bonds involving atoms in mask |
| `printAngles` | `printAngles <mask> [<mask> [<mask>]]` | Print angles involving atoms in mask |
| `printDihedrals` | `printDihedrals <mask> [<mask> [<mask> [<mask>]]]` | Print dihedrals (M=multi-term, I=improper) |
| `printInfo` | `printInfo <flag>` | Print data from a specific %FLAG section |
| `printLJTypes` | `printLJTypes [mask]` | Print LJ type for each atom |
| `printLJMatrix` | `printLJMatrix <mask>` | Print how every type interacts with types in mask |
| `printPointers` | `printPointers` | Print all topology pointers with descriptions |
| `summary` | `summary` | Print system composition summary (residues, charge, mass, density) |
| `listParms` | `listParms` | List all loaded topology files |
| `history` | `history` | Print command history |

#### Modification Commands
| Command | Usage | Description |
|---------|-------|-------------|
| `change` | `change <property> <mask> <value>` | Change CHARGE, MASS, RADII, SCREEN, ATOM_NAME, or AMBER_ATOM_TYPE |
| `changeRadii` | `changeRadii <set>` | Set all GB radii to `bondi`, `mbondi`, `mbondi2`, `mbondi3`, or `amber6` |
| `changeLJPair` | `changeLJPair <m1> <m2> <Rmin> <eps>` | Set pairwise LJ parameters (NBFIX equivalent) |
| `changeLJSingleType` | `changeLJSingleType <mask> <Rmin> <eps>` | Change radius and well depth of a nonbonded type |
| `changeLJ14Pair` | `changeLJ14Pair <m1> <m2> <Rmin> <eps>` | Alter 1-4 LJ terms (chamber topologies only) |
| `addLJType` | `addLJType <mask> [radius <r>] [epsilon <e>] [radius_14 <r14>] [epsilon_14 <e14>]` | Assign new VDW type to atoms |
| `addAtomicNumber` | `addAtomicNumber` | Add ATOMIC_NUMBER section (matched by mass) |
| `addPDB` | `addPDB <file> [elem] [strict] [allicodes]` | Add PDB info (chain ID, residue number, element) |
| `addDihedral` | `addDihedral <m1> <m2> <m3> <m4> <k> <per> <phase> <scee> <scnb> [type <t>]` | Add dihedral term; type=`normal` or `improper` |
| `addExclusions` | `addExclusions <m1> <m2>` | Add arbitrary exclusions |
| `setBond` | `setBond <m1> <m2> <k> <Req>` | Change or add a bond |
| `setAngle` | `setAngle <m1> <m2> <m3> <k> <THETeq>` | Change or add an angle |
| `setMolecules` | `setMolecules [solute_ions=True|False]` | Reset molecularity and SOLVENT_POINTERS |
| `setOverwrite` | `setOverwrite [True|False]` | Allow overwriting original topology |
| `deleteBond` | `deleteBond <m1> <m2> [verbose]` | Delete bonds and all dependent valence terms |
| `deleteDihedral` | `deleteDihedral <m1> <m2> <m3> <m4>` | Delete dihedral (all multiterm terms) |
| `deletePDB` | `deletePDB` | Remove PDB info added by `addPDB` |
| `strip` | `strip <mask> [nobox]` | Remove atoms; `nobox` removes unit cell info |
| `defineSolvent` | `defineSolvent <residue_list>` | Define custom solvent residues (comma-separated) |
| `scale` | `scale <FLAG> <factor>` | Multiply all values in a %FLAG section |
| `scee` | `scee <value>` | Set 1-4 electrostatic scaling constant |
| `scnb` | `scnb <value>` | Set 1-4 VDW scaling constant |

#### HMassRepartition
```bash
HMassRepartition [<mass>] [dowater]
```
| Parameter | Description |
|-----------|-------------|
| `<mass>` | Target hydrogen mass in daltons (default: 3.024) |
| `dowater` | Also repartition water hydrogens (default: solute only) |

Transfers mass from heavy atoms to bonded hydrogens, keeping total mass constant. Enables 4 fs timesteps with SHAKE.

#### Other Key Commands
| Command | Usage | Description |
|---------|-------|-------------|
| `loadRestrt` | `loadRestrt <file>` | Load coordinates from restart file |
| `loadCoordinates` | `loadCoordinates <file>` | Load coordinates (auto-detects format: rst7, ncrst, mdcrd, nc, PDB, mmCIF) |
| `outparm` | `outparm <file>` | Write modified topology to new file |
| `writeFrcmod` | `writeFrcmod <file>` | Dump all parameters as a frcmod file |
| `save` | (alias for `outparm`) | Save topology |
| `energy` | `energy [cutoff <c>] [[igb <m>] [saltcon <c>] \| [Ewald]] [nodisper]` | Compute single-point energy |
| `minimize` | `minimize [cutoff <c>] [[igb <m>] [saltcon <c>]] [[restrain <mask>] [weight <k>]] [norun] [tol <t>] [maxcyc <n>]` | Run minimization |
| `checkValidity` | `checkValidity` | Thorough topology validation |
| `changeProtState` | `changeProtState <mask> <state>` | Change protonation state (constant pH) |
| `changeRedoxState` | `changeRedoxState <mask> <state>` | Change redox state (constant redox) |
| `netCharge` | `netCharge` | Print total system charge |
| `interpolate` | `interpolate <n> [parm2 <p>] [eleconly] [prefix <p>] [startnum <n>]` | Create interpolated topologies (TI) |
| `tiMerge` | `tiMerge <mol1> <mol2> <sc1> <sc2> [<sc1N>] [<sc2N>] [<tol>]` | Merge topologies for TI |
| `source` | `source <file>` | Execute commands from a file |
| `go` | `go` | Execute any pending `parmout` and quit |
| `quit` | `quit` | Exit without executing `parmout` |

## Common Workflows

### Workflow 1: Hydrogen Mass Repartitioning (HMR)
```bash
# Load topology and repartition
parmed diala.parm7 <<EOF
HMassRepartition
outparm diala_hmass.parm7
quit
EOF
```
This triples hydrogen masses to 3.024 daltons and subtracts from heavy atoms. After HMR, use `dt=0.004` with `ntc=2, ntf=2`.

### Workflow 2: Change GB Radii
```bash
parmed system.prmtop <<EOF
changeRadii mbondi3
outparm system_mbondi3.parm7
quit
EOF
```
Use `mbondi3` for `igb=8` (GBneck2), `mbondi2` for `igb=2,5,7`.

### Workflow 3: Modify Charges for FEP Setup
```bash
parmed system.prmtop <<EOF
change charge :3@N -0.3479
change charge :3@H 0.2747
change charge :3@CA -0.24
outparm system_state1.parm7
quit
EOF
```

### Workflow 4: Strip Solvent and Ions
```bash
parmed solvated.prmtop <<EOF
strip :WAT,Na+,Cl- nobox
outparm solute_only.parm7
quit
EOF
```

### Workflow 5: Inspect Before Modifying
```bash
parmed system.prmtop <<EOF
summary
printDetails :1-20
printFlags
netCharge
quit
EOF
```

### Workflow 6: Set NBFIX (Pairwise LJ Correction)
```bash
parmed system.prmtop <<EOF
changeLJPair :ZN :OW 1.45 0.015
outparm system_nbfix.parm7
quit
EOF
```

## Reference Tables

### GB Radii Sets
| Set | Recommended For |
|-----|-----------------|
| `bondi` | Original Pauling radii |
| `mbondi` | `igb=1` (modified Bondi) |
| `mbondi2` | `igb=2`, `igb=5` (OBC models) |
| `mbondi3` | `igb=7`, `igb=8` (GBn, GBneck2) |
| `amber6` | `igb=1` legacy |

### `change` Property Options
| Property | Units | Notes |
|----------|-------|-------|
| `CHARGE` | Elementary charge | Stored internally as e; multiplied by 18.2223 when writing |
| `MASS` | g/mol | |
| `RADII` | Angstroms | GB radii |
| `SCREEN` | -- | GB screening parameters |
| `ATOM_NAME` | -- | Atom name string |
| `AMBER_ATOM_TYPE` | -- | NOT the VDW type |

### `addDihedral` Parameter Defaults
| Parameter | Amber Default | CHARMM/GLYCAM Default |
|-----------|---------------|----------------------|
| `scee` | 1.2 | 1.0 |
| `scnb` | 2.0 | 1.0 |

## Worked Example

### HMR on Alanine Dipeptide

**1. Build the system in tleap:**
```bash
tleap -f leaprc.protein.ff14SB
> diala = sequence {ALA ALA}
> saveAmberParm diala diala.parm7 diala.rst7
> quit
```

**2. Run HMR in parmed:**
```bash
parmed diala.parm7
> HMassRepartition
Repartitioning hydrogen masses to 3.024 daltons.
Not changing water hydrogen masses.
> outparm diala_hmass.parm7
Outputting Amber topology file diala_hmass.parm7
> quit
```

**3. Verify the change:**
```bash
parmed diala_hmass.parm7
> printDetails :*
# Check that hydrogen masses are ~3.024 and heavy atoms are reduced
> summary
Total mass (amu):      <same as before>
> quit
```

**4. Run MD with 4 fs timestep:**
```
Production MD (4 fs with HMR)
 &cntrl
  imin=0, ntx=5, irest=1,
  nstlim=2500000, dt=0.004,
  ntc=2, ntf=2,
  temp0=300.0, ntt=3, gamma_ln=2.0,
  ntb=2, ntp=1,
  cut=8.0,
  ntpr=1000, ntwx=5000, ntwr=50000,
  ioutfm=1, ntxo=2,
  ig=-1,
 /
```

## Key Takeaways
1. **`HMassRepartition`** is the most common parmed workflow -- it enables 4 fs timesteps by redistributing mass to hydrogens (default 3.024 daltons)
2. Use **`changeRadii mbondi3`** when switching to `igb=8` (GBneck2) without rebuilding the topology; use `mbondi2` for `igb=2,5`
3. **`changeLJPair`** is the Amber equivalent of CHARMM NBFIX for tuning specific pairwise Lennard-Jones interactions
4. Always run **`summary`** and **`netCharge`** after modifications to verify system integrity
5. **`strip :WAT nobox`** removes solvent and box information for GB implicit solvent simulations

## Connects To
- **Chapter 3: LEaP and System Building** -- parmed modifies topology files created by tleap
- **Chapter 7: MD Engines and Input Files** -- modified topologies feed into sander/pmemd
- **Chapter 8: Minimization and Relaxation** -- HMR enables larger timesteps during relaxation
- **Chapter 9: Production MD** -- HMR is commonly used in production runs
- Section 15 of the Amber 2026 Reference Manual (PDF pages 285-316)
- Tutorial: `05_building_systems/01-4_Hydrogen Mass Repartitioning`