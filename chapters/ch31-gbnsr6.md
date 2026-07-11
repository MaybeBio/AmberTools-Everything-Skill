# Chapter 30: GBNSR6 -- GB Solvation Model

## Core Commands & Syntax

```bash
gbnsr6 -i mdin -o mdout -p prmtop -c inpcrd
```

| Argument | Description |
|----------|-------------|
| `-i mdin` | Input control data file |
| `-o mdout` | Output file; use `-o stdout` to send to terminal |
| `-p prmtop` | Input molecular topology file |
| `-c inpcrd` | Input initial coordinate file |

## Key Namelists / Input Files

The input file uses two namelists: `&cntrl` and `&gb`.

### `&cntrl` Namelist

| Variable | Description | Default |
|----------|-------------|---------|
| `inp` | = 0: do not compute nonpolar, = 1: compute nonpolar solvation | 0 |

### `&gb` Namelist

| Variable | Description | Default |
|----------|-------------|---------|
| `alpb` | = 0: Canonical GB (Still), = 1: ALPB correction | 1 |
| `epsin` | Solute dielectric constant | 1.0 |
| `epsout` | Solvent dielectric constant | 78.5 |
| `istrng` | Ionic strength (mM); 145 = physiological | 0 |
| `dprob` | Solvent probe radius (Angstrom) | 1.4 |
| `space` | Grid spacing for molecular surface (Angstrom) | 0.5 |
| `arcres` | Arc resolution for numerical integration (Angstrom) | 0.2 |
| `B` | Uniform offset to inverse effective radii (Angstrom^-1) | 0.028 |
| `cavity_surften` | Surface tension for nonpolar (kcal/mol/A^2) | 0.005 |
| `rbornstat` | = 0: no print, = 1: print inverse Born radii | 0 |
| `dgij` | = 0: no pairwise, = 1: print polar pairwise DGij | 0 |
| `chagb` | = 0: no CHA, = 1: use Charge Hydration Asymmetry | 0 |
| `Rs` | Dielectric boundary shift (only with chagb=1) | 0.52 |
| `radiopt` | = 0: hardcoded intrinsic radii, = 1: read from topology | 0 |
| `ROH` | RzOH for CHAGB model (defines water model) | 0.586 |
| `tau` | Tau parameter for CHAGB (effective range of neighboring charges) | 1.47 |

## Available GB Equations

The GBNSR6 program is a standalone tool for GB solvation free energy calculations. It uses a grid-based numerical implementation of the R6 integral for computing effective Born radii.

### GB equation variants:
- **Canonical GB** (Still equation): original formulation
- **ALPB** (Analytical Linearized Poisson-Boltzmann): inexpensive correction that restores correct dielectric constant dependence. Recommended except for small non-spherical molecules or non-singly-connected topologies.
- **CHAGB** (Charge Hydration Asymmetry): incorporates asymmetric response to solvated charge based on explicit water model (e.g., TIP3P). Uses a special intrinsic radii set. Tested on neutral small molecules, charged/uncharged amino acid analogs, and small proteins.

### CHAGB water model parameters (ROH):
| Water Model | ROH (Angstrom) |
|-------------|----------------|
| TIP3P / SPC/E | 0.586 (default) |
| OPC | 0.699 |
| TIP4P | 0.734 |
| TIP5P/E | 0.183 |
| Perfect tetrahedral | 0.0 |

## Numerical Implementation

The R6 integral is performed over a grid-based molecular surface using the field-view method. A uniform Cartesian grid discretizes a rectangular box containing the solute. Surface elements traverse the same solid angle as spherical Lee-Richards molecular surface elements.

**Memory considerations:** For large structures (e.g., nucleosome, ~25,000 atoms), the default grid spacing (0.5 A) requires ~2 GB RAM. Increase `space` to reduce memory footprint at the cost of accuracy.

## Common Workflows

### Electrostatic energy only (default parameters)
```
&cntrl
 inp=0
/
```

### Full solvation (electrostatic + nonpolar) with ALPB
```
&cntrl
 inp=1
/
&gb
 epsin=1.0, epsout=78.5, istrng=0, dprob=1.4, space=0.5,
 arcres=0.2, B=0.028, alpb=1, rbornstat=1, cavity_surften=0.005
/
```

### CHAGB with ALPB correction
```
&cntrl
 inp=1
/
&gb
 alpb=1, chagb=1
/
```

### CHAGB with OPC water model and pairwise decomposition
```
&cntrl
 inp=1
/
&gb
 alpb=1, chagb=1, ROH=0.699, dgij=1
/
```

### Large structure with reduced memory
```
&cntrl
 inp=0
/
&gb
 space=1.0, arcres=0.5, alpb=1
/
```

## Key Takeaways

1. **gbnsr6** is a standalone tool for GB solvation free energy calculations using grid-based numerical R6 integrals.
2. **ALPB correction** (`alpb=1`, default) is recommended for all cases except small non-spherical molecules.
3. **CHAGB** (`chagb=1`) adds charge hydration asymmetry correction; requires matching `ROH` to your water model (TIP3P=0.586, OPC=0.699, TIP4P=0.734).
4. **Memory scales with grid resolution**: adjust `space` (grid spacing) for large systems; default 0.5 A is accurate but memory-intensive.
5. **Pairwise decomposition** (`dgij=1`) provides DGij values without the computational expense of PB-based decomposition.

## Connects To

- Chapter 4: GB models in sander (igb settings)
- Chapter 5: GB models in pmemd
- Chapter 6: PBSA (Poisson-Boltzmann)
- Chapter 29: Torch PBSA (GPU-accelerated PB)
- Chapter 26: LES (GB+LES with igb=1,5,7)