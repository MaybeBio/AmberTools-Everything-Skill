# Chapter 5: Antechamber and GAFF

## Core Commands & Syntax

### antechamber -- primary small-molecule parameterization tool
```bash
antechamber -i <input> -fi <fmt> -o <output> -fo <fmt> [-c <method>] [-nc <net_charge>] [-m <mult>] [-j <idx>] [-s <verb>]
```

Required flags: `-i`, `-o`, `-fi`, `-fo`. All others are optional.

### Input/Output Formats (`-fi` / `-fo`)
| Format  | Index | Description |
|---------|-------|-------------|
| `ac`    | 1     | Antechamber format |
| `mol2`  | 2     | Sybyl Mol2 |
| `pdb`   | 3     | PDB |
| `mpdb`  | 4     | Modified PDB |
| `prepi` | 5     | AMBER PREP (internal) |
| `prepc` | 6     | AMBER PREP (Cartesian) |
| `gcrt`  | 8     | Gaussian Cartesian |
| `gzmat` | 7     | Gaussian Z-Matrix |
| `gout`  | 11    | Gaussian Output |
| `mopint`| 9     | Mopac Internal |
| `mopcrt`| 10    | Mopac Cartesian |
| `mopout`| 12    | Mopac Output |
| `mdl`   | 15    | MDL |
| `rst`   | 17    | AMBER Restart (additional file only) |
| `sqmcrt`| 23    | SQM Input |
| `sqmout`| 24    | SQM Output |
| `gesp`  | 26    | Gaussian ESP |

### Charge Methods (`-c`)
| Method | Description |
|--------|-------------|
| `bcc`  | AM1-BCC (fast, recommended default) |
| `resp` | RESP (requires Gaussian) |
| `gas`  | Gasteiger |
| `rc`   | Read in charge from input |
| `wc`   | Write out charge |
| `mul`  | Mulliken |
| `cm1`  | CM1 |
| `cm2`  | CM2 |
| `esp`  | ESP (Kollman) |
| `dc`   | Delete Charge |

### Key Flags
| Flag | Description |
|------|-------------|
| `-nc <int>` | Net molecular charge (overrides file value) |
| `-m <int>` | Multiplicity (2S+1), default 1 |
| `-j <0-5>` | Atom/bond type prediction index. 4=atom+full bond (default); 0=no assignment; 1=atom type only; 2=full bond; 3=part bond; 5=atom+part bond |
| `-s <0-2>` | Status verbosity: 0=brief, 1=default, 2=verbose |
| `-at <type>` | Atom type: `gaff` (default), `gaff2`, `amber`, `bcc`, `sybyl` |
| `-rn <name>` | Residue name (1-3 chars), overrides input |
| `-eq <0-2>` | Charge equalization: 0=none, 1=atomic paths (default), 2=paths+geometry (E/Z config) |
| `-pf` | Remove intermediate files: `y` or `n` (default: keep) |
| `-dr` | AC doctor mode: `y` (default) or `n` |
| `-du` | Fix duplicate atom names: `y` (default) or `n` |
| `-an` | Adjust atom names: `y` for mol2/ac, `n` for others |
| `-ch <name>` | Check file name for Gaussian |
| `-gk <kw>` | mopac or sqm keyword (inside quotes) |
| `-gv <0/1>` | Gaussian version: 1=G09 (default), 0=older |
| `-ge <file>` | Gaussian ESP file name (iop(6/50)=1) |
| `-pl <int>` | Max path length for charge equivalence (set 10-30 for >= 100 atoms) |
| `-a <file>` | Additional file name |
| `-fa <fmt>` | Additional file format |
| `-ao <op>` | Additional file operation: `crd`, `crg`, `radius`, `name`, `type`, `bond` |

### parmchk2 -- missing parameter detection
```bash
parmchk2 -i <input> -f <fmt> -o <frcmod> [-s <ff_set>] [-a Y|N] [-w Y|N]
```

| Flag | Description |
|------|-------------|
| `-i` | Input file name |
| `-o` | frcmod output file name |
| `-f` | Input format: `prepi`, `prepc`, `ac`, `mol2`, `frcmod`, `leaplog` |
| `-s` | FF parameter set: `1` or `gaff` (default), `2` or `gaff2`, `3` or `parm99`, `4` or `parm10`, `5` or `lipid14` |
| `-a` | Print ALL parameters: `Y` or `N` (default: N) |
| `-w` | Print parameters matching improper dihedrals with `X`: `Y` (default) or `N` |
| `-p` | Custom parmfile (suppresses `-s`) |
| `-frc` | Load frcmods: `ff14SB`, `ff99SB`, `bsc1`, `ol15`, `yil`, etc. (e.g. `ff14SB+bsc1+yil`) |
| `-c` | Atom type correspondence score file (default: `PARMCHK.DAT`) |

### Other Antechamber Suite Programs
| Program | Purpose |
|---------|---------|
| `atomtype` | Assigns GAFF/AMBER atom types |
| `bondtype` | Assigns bond types |
| `prepgen` | Generates residue topology (prep) files |
| `respgen` | Generates RESP charge input files |
| `espgen` | Extracts ESP data from Gaussian output |
| `sqm` | Semi-empirical QM program (`sqm -i <input> -o <output>`) |

## Key Namelists / Input Files

### GAFF2 Atom Types
GAFF2 atom types are all lowercase to distinguish from uppercase AMBER macromolecular types. Common types:
| Type | Description |
|------|-------------|
| `c` | sp2 carbon (C=O) |
| `c1` | sp carbon |
| `c2` | sp2 carbon (alkene) |
| `c3` | sp3 carbon |
| `ca` | aromatic carbon |
| `cc` | carboxylate carbon |
| `ce` | sp2 carbon (guanidine) |
| `cf` | sp2 carbon (CF3) |
| `cg` | sp2 carbon (conjugated) |
| `cp` | sp2 carbon (head of chain) |
| `cq` | sp2 carbon (head of chain, conjugated) |
| `cu` | sp2 carbon (urea) |
| `cv` | sp2 carbon (urea conjugated) |
| `cx` | sp3 carbon (triangular) |
| `cy` | sp3 carbon (triangular conjugated) |
| `h1` | H on alkene carbon |
| `h2` | H on alkane carbon |
| `h3` | H on electron-withdrawing carbon |
| `h4` | H on aromatic/aliphatic carbon |
| `h5` | H on amine |
| `ha` | aromatic H |
| `hc` | H on sp3 carbon |
| `hn` | H on amide nitrogen |
| `ho` | hydroxyl H |
| `hp` | H on terminal acetylene |
| `hs` | thiol H |
| `hw` | H on water |
| `hx` | H on triangular carbon |
| `n` | sp2 nitrogen |
| `n1` | sp nitrogen |
| `n2` | sp2 nitrogen (amide) |
| `n3` | sp3 nitrogen |
| `n4` | sp3 nitrogen (ammonium) |
| `na` | sp2 nitrogen (aromatic) |
| `nh` | sp2 nitrogen (amine) |
| `no` | sp2 nitrogen (nitro) |
| `o`  | sp2 oxygen (C=O) |
| `oh` | sp3 oxygen (hydroxyl) |
| `os` | sp3 oxygen (ether, ester) |
| `ow` | sp3 oxygen (water) |
| `p2` | sp2 phosphorus |
| `p3` | sp3 phosphorus |
| `p4` | sp3 phosphorus (phosphate) |
| `p5` | sp3 phosphorus (phosphonate) |
| `s`  | sp2 sulfur |
| `s2` | sp2 sulfur (sulfoxide) |
| `s4` | sp3 sulfur (sulfone) |
| `s6` | sp3 sulfur (sulfonate) |
| `sh` | sp3 sulfur (thiol) |
| `ss` | sp2 sulfur (thioether) |
| `f`  | fluorine |
| `cl` | chlorine |
| `br` | bromine |
| `i`  | iodine |

## Common Workflows

### Step 1: Prepare ligand structure
```bash
# Add hydrogens with reduce
reduce ligand.pdb > ligand_h.pdb

# Or with OpenBabel
obabel ligand.pdb -O ligand_h.pdb -h
```

### Step 2: Run antechamber with AM1-BCC charges
```bash
antechamber -i ligand_h.pdb -fi pdb -o ligand.mol2 -fo mol2 -c bcc -nc 0 -s 2
```

### Step 3: Check for missing parameters
```bash
parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod
```
Review the output frcmod file. Look for `"ATTN: needs revision"` entries -- these require manual parameterization.

### Step 4: Load into tleap
```bash
tleap -f oldff/leaprc.ff14SB
# In tleap:
source leaprc.gaff
loadAmberParams ligand.frcmod
LIG = loadMol2 ligand.mol2
complex = combine {protein LIG}
saveAmberParm complex complex.prmtop complex.inpcrd
quit
```

### Step 5: For RESP charges (requires Gaussian)
```bash
# Generate Gaussian input
antechamber -i ligand.mol2 -fi mol2 -o gcrt.com -fo gcrt -gv 1 -ge ligand.gesp

# Run Gaussian09 with gcrt.com as input, then:
antechamber -i ligand.gesp -fi gesp -o ligand_resp.mol2 -fo mol2 -c resp -eq 2
```

### Step 6: Check sqm calculation success
```bash
tail sqm.out
# Should show: "Calculation Completed"
```

## Reference Tables

### Charge Method Selection Guide
| Scenario | Recommended Method |
|----------|-------------------|
| Fast, automated workflow | `-c bcc` (AM1-BCC) |
| High accuracy, have Gaussian | `-c resp` (RESP) |
| Quick estimate | `-c gas` (Gasteiger) |
| Read existing charges from file | `-c rc` |
| Mulliken population analysis | `-c mul` |

### `-j` Flag Behavior
| `-j` | Atom Type | Bond Type |
|------|-----------|-----------|
| 0 | None | None |
| 1 | Assigned | None |
| 2 | None | Full |
| 3 | None | Part |
| 4 (default) | Assigned | Full |
| 5 | Assigned | Part |

### Input Format Detection
Antechamber auto-detects the following from the input directly:
- Net charge: from `gout`, `mopout`, `sqmout`, `sqmcrt`, `gcrt`
- Calculated from partial charge sum: from `mol2`, `prepi`

## Worked Example

### Simulating Sustiva (efavirenz) with GAFF

**1. Prepare the PDB:**
```bash
reduce sustiva.pdb > sustiva_h.pdb
```

**2. Run antechamber:**
```bash
antechamber -i sustiva_new.pdb -fi pdb -o sustiva.mol2 -fo mol2 -c bcc -s 2
```
This executes internally: `bondtype` -> `atomtype` -> `sqm` -> `am1bcc` -> `atomtype`. The output `sustiva.mol2` contains Tripos format with GAFF atom types (lowercase) and BCC charges.

**3. Check parameters:**
```bash
parmchk2 -i sustiva.mol2 -f mol2 -o sustiva.frcmod
```
A typical frcmod output shows missing angle parameters (e.g., `ca-c3-c1`, `c1-c1-cx`) and improper dihedrals (e.g., `ca-ca-ca-ha`). If any entry shows `"ATTN: NEEDS REVISION"`, the parameter must be manually supplied.

**4. Build in tleap:**
```bash
tleap -f oldff/leaprc.ff99SB
> source leaprc.gaff
> loadAmberParams sustiva.frcmod
> SUS = loadMol2 sustiva.mol2
> complex = combine {protein SUS}
> saveAmberParm complex complex.prmtop complex.inpcrd
> quit
```

**5. Run GB simulation:**
```bash
pmemd.cuda -O -i md.in -o md.out -p complex.prmtop -c complex.inpcrd \
  -r md.rst7 -x md.nc -inf md.info
```

## Key Takeaways
1. **`antechamber -c bcc`** is the standard fast charge method for automated workflows; use **`-c resp`** for publication-quality charges when Gaussian is available
2. Always run **`parmchk2`** after antechamber to identify missing parameters; critically review any `"ATTN: needs revision"` entries
3. GAFF2 atom types are lowercase to keep them independent from uppercase AMBER macromolecular atom types, enabling mixed force field simulations
4. Use **`-j 4`** (default) for full atom and bond type assignment; use **`-j 5`** if the input already has correct connectivity
5. The mol2 output file is the standard input format for LEaP; load it with `loadMol2` after sourcing `leaprc.gaff`

## Connects To
- **Chapter 3: LEaP and System Building** -- loading GAFF molecules into tleap, combining with protein force fields
- **Chapter 7: MD Engines and Input Files** -- running simulations with the generated prmtop/inpcrd files
- **Chapter 8: Minimization and Relaxation** -- relaxing ligand-protein complexes before production MD
- **Chapter 9: Production MD** -- production simulations of ligand-bound systems
- Section 16 of the Amber 2026 Reference Manual (PDF pages 322-338)
- Tutorial: `06_Developing_Nonstandard_Parameters/02-1_Simulating a pharmaceutical compound with Antechamber and GAFF`