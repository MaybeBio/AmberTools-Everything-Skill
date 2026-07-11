# Chapter 26: LES -- Locally-Enhanced Sampling

## Core Commands & Syntax

LES preparation is done via the **ADDLES** program, which modifies a standard prmtop to create the LES system:

```bash
addles < inputfile > outputfile
```

### ADDLES Input File Format

```
~ lines beginning with ~ are comments
~ all commands are 4 letters
~ max line length is 80 characters; trailing "-" is continuation

file rprm name=(prmtop) read
file rcrd name=(coords) read
~ or rcvd (coords+velocities), rcvb (coords+velo+box), rcbd (coords+box)
file wprm name=(lesparm) wovr
file wcrd name=(lescrd) wovr

action
omas
~ omas: leave masses at original values (otherwise scaled by 1/N)

spac numc=5 pick #mon 1 2 done
spac numc=3 pick #sid 3 3 done

*EOD
```

### PICK Commands

| Syntax | Description |
|--------|-------------|
| `#prt A B` | Atom range A to B by atom number |
| `#mon A B` | Residue range A to B by residue number |
| `#cca A B` | Residue range A to B, dividing between CA and C of residue B |
| `#sid A B` | Sidechain atoms of residue A to B (excludes C,O,CA,HA,N,H,HN) |
| `chem prtc A` | All atoms named A (case-sensitive) |
| `chem mono A` | All residues named A (case-sensitive) |

Wildcards: `H*` matches H, HA, etc.

Boolean operators: `|` (or), `&` (and), `!=` (not). Evaluated in order; no parentheses.

### ADDLES Options

| Command | Description |
|---------|-------------|
| `file` | Open a file |
| `rprm` | Read prmtop |
| `rcrd` | Read coordinates |
| `rcvd` | Read coordinates + velocities |
| `rcvb` | Read coordinates + velocities + box |
| `rcbd` | Read coordinates + box (pack=n option for multiple coord sets) |
| `wprm` | Write new topology file |
| `wcrd` | Write coordinates |
| `wovr` | Overwrite existing file |
| `writ` | Do not overwrite existing file |
| `action` | Start processing (all commands after this) |
| `omas` | Leave masses at original values (default: scale by 1/N) |
| `nomodv` | Do NOT randomize copy velocities (default: slightly randomize) |
| `spac` | Define a subspace: `numc=N pick ... done` |
| `pimd` | Write prmtop for PIMD (no exclusion between copies) |
| `*EOD` | End of input |

### Pack option for different initial coordinates
```bash
file rcrd name=(input.inpcrd) pack=4 read
```
Distributes 4 coord sets among copies: e.g., 20 copies with 4 sets -> copies 1-5 get set 1, 6-10 get set 2, etc.

## SANDER LES Variables

The LES version of sander (`sander.LES`) handles the special scaling. Key modifications:

- All force constants, charges, and VDW epsilon values are scaled by 1/N in the prmtop file.
- LES copies of the same subspace are placed in exclusion lists (no self-interaction).
- Intra-copy interactions are corrected: for N copies in a subspace, diagonal interactions are scaled up by N.
- PME requires separate correction calculations for excluded atom pairs and intra-copy interactions.

**&cntrl variable for LES:**

| Variable | Description |
|----------|-------------|
| `temp0les` | Target temperature for LES particles. If `<0`, single temperature bath for all atoms; otherwise separate thermostats for LES and non-LES atoms. |
| `les` | Flag to enable LES (used in sander.LES executable) |

## LES with Generalized Born Solvation

GB models supported: `igb = 1, 5, or 7`. Surface area calculations not supported with LES. Only a single LES region permitted for GB+LES.

| Variable | Description | Default |
|----------|-------------|---------|
| `rdt` | Effective radii deviation threshold. When multiple Born radii for non-LES atom differ by less than RDT, a single radius is used. | 0.01 |

The GB+LES approach solves the key explicit-water LES problem: with explicit solvent, water interacts in an average way with all copies, causing copies to resist moving apart (cavity creation penalty) and preventing individual solvation. GB provides independent solvation for each copy.

## Design Considerations

### What to copy
- Flexible regions of interest: loops, side chains, entire molecules
- Do NOT copy water molecules
- Copies should not span different molecules for pressure coupling

### How many copies
- 3-10 copies are reasonable; 5 is a good starting point
- More copies: flatter energy surface, wider sampling, but more approximate
- Barrier heights reduced roughly proportional to number of copies

### How many regions
- Larger regions: more independent copies, wider conformational variety, less surface smoothing
- Smaller regions: more inter-copy averaging, better barrier reduction
- At least 2 residues per region when copying backbone
- Rule of thumb: several LES regions (unless copying a very small region)

### Order of operations
1. Build system in LEaP (add solvent, ions)
2. Equilibrate without LES
3. Run ADDLES (LAST step before sander)
4. Run sander.LES

## Common Workflows

### Example 1: Multiple hydroxyl groups (glucose)
```
file rprm name=(parm.solv.top) read
file rcvb name=(glucose.solv.equ.crd) read
file wprm name=(les.prmtop) wovr
file wcrd name=(glucose.les.crd) wovr
action
omas
spac numc=5 pick chem prtc HO1 done
spac numc=5 pick chem prtc HO2 done
spac numc=5 pick chem prtc HO3 done
spac numc=5 pick chem prtc HO4 done
spac numc=5 pick #prt 20 24 done
*EOD
```

### Example 2: RNA tetraloop
```
file rprm name=(prm.top) read
file rcvb name=(rna.crd) read
file wprm name=(les.parm) wovr
file wcrd name=(les.crd) wovr
action
omas
spac numc=5 pick #prt 131 255 done
*EOD
```

### Example 3: Small peptide (AVPA) with 2 regions
```
file rprm name=(prmtop) read
file rcvb name=(md.solv.crd) read
file wprm name=(LES.prmtop) wovr
file wcrd name=(LES.crd) wovr
action
omas
spac numc=5 pick #cca 1 3 | #mon 1 1 done
spac numc=5 pick #cca 4 6 | #mon 6 6 done
*EOD
```

## Key Takeaways

1. **ADDLES** modifies a standard prmtop to create an LES system by copying atoms into subspaces with scaled force constants.
2. **PICK commands** select atoms by residue range (`#mon`), atom range (`#prt`), C-alpha divisions (`#cca`), sidechain (`#sid`), or chemical name (`chem`).
3. **Design choices matter**: what to copy, how many copies (3-10), and how many regions determine effectiveness.
4. **GB+LES** solves explicit-water limitations by providing independent solvation for each copy (igb=1,5,7 only).
5. **Analysis is challenging**: use LES-compatible tools (MOIL-View) or extract single copies from trajectories.

## Connects To

- Chapter 2: LEaP (initial system building)
- Chapter 4: sander (LES uses sander.LES executable)
- Chapter 5: GBNSR6 (GB solvation with LES)
- Chapter 34: gem.pmemd (AMOEBA force field)
- Chapter 35: AI/ML (KMMD for enhanced sampling)