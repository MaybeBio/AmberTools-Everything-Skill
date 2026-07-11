# Chapter 27: NAB -- Nucleic Acid Builder and libsff

## Core Commands & Syntax

### NAB Compilation Pipeline

```bash
# The NAB compiler nab2c converts NAB source to C
nab2c -i input.nab -o output.c

# Compile the generated C code with libsff
gcc -o output output.c -I$AMBERHOME/include -L$AMBERHOME/lib -lsff -lm

# Run the compiled program
./output
```

The modern `nabc` provides a C-language interface to libsff:

```bash
# nabc is the C-language interface to libsff
# Available from https://github.com/dacase/nabc
# Requires AmberTools to be installed
```

### NAB Language Basics

NAB (Nucleic Acid Builder) was designed as a "molecular awk" -- a scripting language for manipulation of macromolecules. It was created by Tom Macke in the 1990s. The NAB compiler (`nab2c`) converts NAB source code to C for subsequent compilation.

The SFF ("simple force field") routines were added by David Case, providing Amber-compatible molecular mechanics. Over time, sff evolved to include generalized Born, Poisson-Boltzmann, and RISM approaches, concentrating on implicit solvation. The sff routines were parallelized using both OpenMP and MPI, and second derivatives of the generalized Born model were added by Russ Brown.

### Key NAB Data Types

| Type | Description |
|------|-------------|
| `molecule` | A molecular system (atoms, residues, bonds, parameters) |
| `residue` | A single residue within a molecule |
| `atom` | An individual atom with properties |
| `bond` | Covalent bond between two atoms |
| `angle` | Bond angle between three atoms |
| `dihedral` | Dihedral (torsion) angle between four atoms |
| `improper` | Improper dihedral for planarity |
| `vdw` | Van der Waals parameters |
| `charge` | Partial atomic charge |
| `mass` | Atomic mass |
| `element` | Chemical element |
| `type` | Amber atom type |
| `label` | Atom label |

### NAB Commands for Molecule Building

```
# Building DNA/RNA
molecule m;
m = bdna("gcgcaattcgcg");          # Build B-DNA duplex
m = arna("gcgcaattcgcg");          # Build A-RNA duplex
m = link_na("gcg", "dna", "bform", "5p", "3p");  # Link nucleic acid

# Reading/writing structures
m = getpdb("input.pdb");           # Read PDB file
putpdb("output.pdb", m);           # Write PDB file
savepdb("output.pdb", m);          # Alternative PDB write
loadpdb("input.pdb", m);           # Alternative PDB read

# Parameter manipulation
putbnd("CT", "CT", 310.0, 1.53);   # Set bond: types, k, r0
putang("CT", "CT", "CT", 50.0, 109.5);  # Set angle
putdih("CT", "CT", "CT", "CT", 1.0, 1, 180.0);  # Set dihedral
putimp("C", "O", "N", "CT", 10.5, 180.0, 0.0);  # Set improper
putvdw("CT", 1.9080, 0.1094);      # Set VDW: type, R*, epsilon
putels("CT", -0.1824);             # Set electrostatic charge

# System setup
addions(m, "Na+", 10);             # Add 10 Na+ counterions
addions(m, "Cl-", 9);              # Add 9 Cl- counterions
setbox(m, "rect", 10.0, 10.0, 10.0);  # Set rectangular periodic box
setsugar(m, "dna");                # Set sugar pucker to DNA-like
setpdb(m);                         # Set PDB coordinates
```

### NAB Force Field Functions (libsff)

```c
// Force field lifecycle
mm_options("mdin");                // Read mdin-style options file
mm_init(m);                        // Initialize force field for molecule
mm_create(m);                      // Create force field context
mm_run(m);                         // Execute the simulation
mm_free(m);                        // Free force field resources

// Energy/force evaluation
float e = mm_energy(m);            // Calculate total potential energy
vec3 *f = mm_force(m);             // Calculate forces on all atoms
vec3 *g = mm_gradient(m);          // Calculate gradient
mat *h = mm_hessian(m);            // Calculate Hessian matrix

// Minimization and dynamics
mm_minimize(m, "min.in");          // Energy minimization
mm_dynamics(m, "md.in");           // Molecular dynamics

// Ensemble control
// Available ensembles: NVE, NVT, NPT
// Thermostats: Langevin, Berendsen, Andersen, Nose-Hoover, Nose-Hoover chain
// Barostats: Berendsen, Monte Carlo
```

### mdin Options for NAB/sff (from the manual)

The sff routines use Amber-style mdin files. Key options for PB dynamics in NAB (Section 6.6.2 of the manual):

```
&PBSA input for PBMD with NAB
&cntrl
  imin=0, ntx=1, irest=0,
  ipb=2, ntb=0,
  ntc=2, ntf=2,
  tempi=100, temp0=100, ntt=3, gamma_ln=1,
  nstlim=100000, dt=0.002,
  ntpr=100, ntwr=100000, ntwx=100,
/
&pb
  npbgrid=500, nsnba=5,
/
```

For PB minimization/dynamics, the recommended PB options are:
```
space=0.25
arcres=0.125
fscale=4
eneopt=2
bcopt=6
frcopt=2
```

IPB is explicitly set to 2 to enable PB dynamics. NPBGRID=500 means the finite-difference grid is regenerated every 500 dynamics steps. NSNBA=5 means the atom-based pairlist is generated every 5 steps.

## Key Namelists / Input Files

### NAB mdin Options (cntrl namelist subset)

The sff library reads Amber-style mdin files. Key control variables available in NAB/sff:

| Variable | Description | Default |
|----------|-------------|---------|
| imin | 0=MD, 1=minimization | 0 |
| ntx | Coordinate format option | 1 |
| irest | Restart flag | 0 |
| ntb | Periodic box: 0=none, 1=constant vol, 2=constant P | 0 |
| ipb | PBSA: 0=off, 1=PB energy, 2=PB forces | 0 |
| ntc | SHAKE: 1=none, 2=bonds to H, 3=all bonds | 2 |
| ntf | Force evaluation: 1=complete, 2=no H-bond | 2 |
| tempi | Initial temperature (K) | 0 |
| temp0 | Target temperature (K) | 300 |
| ntt | Thermostat: 0=none, 1=Berendsen, 2=Andersen, 3=Langevin | 0 |
| gamma_ln | Langevin collision frequency (ps^-1) | 0 |
| nstlim | Number of MD steps | 1000 |
| dt | Time step (ps) | 0.001 |
| ntpr | Print frequency (steps) | 50 |
| ntwr | Restart write frequency | 500 |
| ntwx | Trajectory write frequency | 0 |
| cut | Nonbonded cutoff (A) | 8.0 |
| igb | GB model: 0=none, 1=GB^OBC, 2=GB^OBC2, 5=GB^OBC2, 7=GB-neck2, 8=GB-neck2 | 0 |
| saltcon | Salt concentration (M) | 0.0 |

### PB Namelist (for NAB/sff)

| Variable | Description | Default |
|----------|-------------|---------|
| npbgrid | Steps between grid regeneration | 100 |
| nsnba | Steps between pairlist rebuild | 11 |
| space | Grid spacing (A) | 0.5 |
| arcres | Arc resolution (A) | 0.25 |
| fscale | Scale factor for grid | 2 |
| eneopt | Energy calculation method | 1 |
| bcopt | Boundary condition option | 5 |
| frcopt | Force calculation method | 1 |
| epsin | Solute dielectric constant | 1.0 |
| epsout | Solvent dielectric constant | 80.0 |
| istrng | Salt concentration (mM) | 0 |
| radiopt | Radius option | 0 |
| radiscale | Scaling factor for radii | 1.0 |
| protscale | Protein radii scaling factor | 1.0 |

## Common Workflows

### Building a DNA Duplex and Minimizing

```c
// dna_min.nab
molecule m;
float e;

// Step 1: Build B-DNA
m = bdna("gcgcaattcgcg");

// Step 2: Write initial PDB
putpdb("dna_initial.pdb", m);

// Step 3: Minimize with sff (GB implicit solvent)
m = minimize(m, "min.in");

// Step 4: Write minimized PDB
putpdb("dna_min.pdb", m);
```

min.in:
```
&cntrl
  imin=1, maxcyc=500, ncyc=250,
  ntb=0, igb=1, saltcon=0.1,
  cut=999.0,
/
```

### PB Dynamics with NAB (from Section 6.6.2)

```bash
# The example input listed in the manual for PBMD in sander/NAB:
# PB visualization input
cat > pbmd.in << 'EOF'
&cntrl
  imin=0, ntx=1, irest=0,
  ipb=2, ntb=0,
  ntc=2, ntf=2,
  tempi=100, temp0=100, ntt=3, gamma_ln=1,
  nstlim=100000, dt=0.002,
  ntpr=100, ntwr=100000, ntwx=100,
/
&pb
  npbgrid=500, nsnba=5,
/
EOF
```

This input enables PB dynamics (ipb=2), with the grid regenerated every 500 steps and the pairlist every 5 steps. The simulation runs at 100 K with Langevin thermostat, no periodic box, and a 2 fs timestep.

### Electrostatic Force Calculation in PB

The finite-difference Poisson-Boltzmann method in sff/NAB computes three force components:
1. **Reaction field force** -- exists where atomic charges are present; straightforward to map onto atoms
2. **Dielectric boundary force** -- exists on the molecular surface where dielectric constant changes; harder to map
3. **Ionic force** -- much smaller in magnitude, not included in the release

The charge view method (ENEOPT=2, FRCOPT=2) is recommended for stable MD simulations. BCOPT=6 removes charge singularity for stability.

## Reference Tables

### NAB File Types

| Extension | Type | Description |
|-----------|------|-------------|
| .nab | NAB source | NAB scripting language source code |
| .c | C source | Generated C code from nab2c |
| .pdb | PDB file | Protein Data Bank format |
| .crd | Coordinate | Amber coordinate file |
| .rst / .rst7 | Restart | Amber restart file |
| .prmtop / .parm7 | Topology | Amber topology/parameter file |
| .inpcrd | Coordinate | Amber input coordinates |
| .off | OFF library | Amber object file format library |
| .frcmod | Force mod | Amber force field modification |
| .prep / .prepi | Prep | Amber prep/residue template |
| .lib | Library | Amber library file |
| .dat | Data | Generic data file |
| .leaprc | LEaP resource | LEaP configuration file |

### NAB Molecule Building Functions

| Function | Purpose |
|----------|---------|
| `bdna(seq)` | Build canonical B-form DNA duplex from sequence |
| `arna(seq)` | Build canonical A-form RNA duplex from sequence |
| `link_na(seq, type, form, term5, term3)` | Build nucleic acid with specified parameters |
| `getpdb(file)` | Read molecule from PDB file |
| `putpdb(file, mol)` | Write molecule to PDB file |
| `putbnd(t1, t2, k, r0)` | Set bond force constant and equilibrium |
| `putang(t1, t2, t3, k, theta0)` | Set angle force constant and equilibrium |
| `putdih(t1, t2, t3, t4, k, n, delta)` | Set dihedral parameters |
| `putimp(t1, t2, t3, t4, k, delta, phase)` | Set improper dihedral parameters |
| `putvdw(type, R, eps)` | Set van der Waals radius and well depth |
| `putels(type, charge)` | Set atomic partial charge |
| `addions(mol, ion, n)` | Add n counterions of given type |
| `setbox(mol, shape, x, y, z)` | Set periodic box dimensions |
| `setsugar(mol, type)` | Set sugar pucker (dna/rna) |
| `setpdb(mol)` | Set PDB coordinates to current positions |

### sff Force Field Functions

| Function | Purpose |
|----------|---------|
| `mm_options(file)` | Read mdin-style options |
| `mm_init(mol)` | Initialize force field |
| `mm_create(mol)` | Create force field context |
| `mm_run(mol)` | Run simulation |
| `mm_free(mol)` | Free resources |
| `mm_energy(mol)` | Compute total potential energy |
| `mm_force(mol)` | Compute forces |
| `mm_gradient(mol)` | Compute gradient |
| `mm_hessian(mol)` | Compute Hessian |
| `mm_minimize(mol, file)` | Run energy minimization |
| `mm_dynamics(mol, file)` | Run molecular dynamics |

### NAB/sff Capabilities

| Feature | Support |
|---------|---------|
| GB implicit solvent | Full (igb=1,2,5,7,8) |
| PB implicit solvent | Full (ipb=1,2) |
| RISM solvation | Supported |
| Explicit solvent | Not supported (use sander/pmemd) |
| Periodic boundary | Not supported |
| GPU acceleration | Not available |
| OpenMP parallelization | Supported |
| MPI parallelization | Supported |
| Second derivatives (GB) | Supported (Brown) |
| Hierarchical charge partition | Supported |
| Advanced nonbonded list | Supported |
| LMOD (low-mode search) | Supported |

## Worked Example

### Build, Solvate, and Minimize a DNA Hairpin

From the manual's conceptual framework, a complete NAB workflow:

```c
// hairpin.nab -- Build and minimize a DNA hairpin
molecule m, wat;
int i;
float e;

// Step 1: Build the DNA hairpin 
// (sequence includes loop region)
m = link_na("gcgc tttt gcgc", "dna", "bform", "5p", "3p");

// Step 2: Add counterions to neutralize
addions(m, "Na+", 22);

// Step 3: Write initial structure
putpdb("hairpin_initial.pdb", m);

// Step 4: Load force field options
mm_options("min_gb.in");

// Step 5: Initialize and minimize
mm_init(m);
mm_create(m);
mm_run(m);

// Step 6: Get final energy and write structure
e = mm_energy(m);
printf("Final energy: %f kcal/mol\n", e);
putpdb("hairpin_min.pdb", m);

mm_free(m);
```

min_gb.in:
```
&cntrl
  imin=1, maxcyc=1000, ncyc=500,
  ntpr=50, ntb=0,
  igb=5, saltcon=0.1,
  cut=999.0,
/
```

## Key Takeaways

1. **NAB is a legacy scripting language** for molecular manipulation, originally designed as a "molecular awk." The compiler `nab2c` converts NAB source to C. The code has been moved to https://github.com/dacase/nabc and requires AmberTools.

2. **libsff is the Amber-compatible force field library** used by NAB. It supports implicit solvent models (GB, PB, RISM), OpenMP/MPI parallelization, and advanced features like hierarchical charge partitioning and second derivatives. It does NOT support explicit solvent, periodic boundaries, or GPU acceleration.

3. **NAB excels at nucleic acid building** with functions like `bdna()`, `arna()`, and `link_na()` that create canonical A/B-form helices from sequence strings. Parameter manipulation functions (`putbnd()`, `putang()`, etc.) allow custom force field modifications.

4. **PB dynamics in NAB** uses the same input format as sander. For stable PBMD, set `ipb=2`, and use the recommended PB options: `space=0.25`, `arcres=0.125`, `fscale=4`, `eneopt=2`, `bcopt=6`, `frcopt=2`. The charge view method is preferred for stability.

5. **The dielectric boundary force** in PB cannot be easily mapped to atoms because the molecular surface derivatives do not exist for the solvent excluded surface (SES) definition. A partial solution using the Gilson et al. mapping method is implemented.

## Connects To

- [ch03](ch03-leap-system-building.md) -- LEaP system building (alternative to NAB for system preparation)
- [ch07](ch07-md-engines-input.md) -- MD engines and mdin reference (cntrl namelist shared with sff)
- [ch08](ch08-minimization-relaxation.md) -- Minimization protocols
- [ch20](ch20-implicit-solvent.md) -- Implicit solvent models (GB, PB used by NAB/sff)
- [ch22](ch22-nmr-cryoem-nab.md) -- NMR/CryoEM/NAB (existing NAB coverage)
- [ch28](ch28-bar-pbsa.md) -- BAR/PBSA (PB-based free energy analysis)