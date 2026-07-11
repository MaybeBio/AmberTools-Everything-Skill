# Chapter 29: RISM -- 1D and 3D Reference Interaction Site Model

## Core Commands & Syntax

RISM provides an implicit solvent framework based on the Ornstein-Zernike integral equation. It comes in two variants: 1D-RISM for bulk solvent properties and 3D-RISM for solute-solvent distributions.

### 1D-RISM: Bulk Solvent Calculations

```bash
rism1d inputfile
```

The input file must have a `.inp` suffix. Parameters are specified in Fortran namelist format.

### 3D-RISM: Single-Point (standalone)

```bash
rism3d.snglpnt --pdb solute.pdb --prmtop solute.prmtop \
    --xvv solvent.xvv --closure "KH" --buffer 14 --grdspc 0.5,0.5,0.5

# MPI version
mpirun -np 4 rism3d.snglpnt.MPI --pdb solute.pdb --prmtop solute.prmtop \
    --xvv solvent.xvv --closure "KH" --buffer 14 --grdspc 0.5,0.5,0.5
```

### 3D-RISM: sander interface

```bash
sander -O -i mdin -o mdout -p prmtop -c inpcrd \
    -xvv solvent.xvv -guv guv -huv huv -cuv cuv
```

In `&cntrl`, set `irism=1`. All 3D-RISM keywords go in the `&rism` namelist.

### 3D-RISM: Next-Generation (rism3d_ng)

```bash
rism3d_ng --pdb solute.pdb --prmtop solute.prmtop \
    --xvv solvent.xvv [other options]

# MPI version
mpirun -np N rism3d_ng.MPI [options]
```

---

## Theoretical Background

### The Ornstein-Zernike Equation

RISM approximates the OZ equation:

```
h(r12,Omega1,Omega2) = c(r12,Omega1,Omega2) + rho * integral[ dr3 dOmega3 * c(r13,Omega1,Omega3) * h(r32,Omega3,Omega2) ]
```

For infinite dilution, this splits into three equations: solvent-solvent (VV), solute-solvent (UV), and solute-solute (UU). The 1D-RISM equation becomes:

```
h_alpha_gamma(r) = sum_{mu,nu} integral[ dr' dr'' * omega_alpha_mu(|r-r'|) * c_mu_nu(|r'-r''|) * (omega_nu_gamma(r'') + rho_nu * h_nu_gamma(r'')) ]
```

### Closure Relations

The general closure relation is:

```
g(r12,Omega1,Omega2) = exp[-beta*u(r12,Omega1,Omega2) + h(r12,Omega1,Omega2) - c(r12,Omega1,Omega2) + b(r12,Omega1,Omega2)]
```

Three closures are implemented in 3D-RISM:

| Closure | Formula (g > 1 region) | Characteristics |
|---------|----------------------|-----------------|
| **HNC** | `exp(t*)` | Good for charged particles; may not converge with asymmetric sizes |
| **KH** (default) | `1 + t*` | Numerically robust; reliably converges; overestimates non-Coulombic |
| **PSE-n** | `sum_{i=0}^n (t*)^i / i!` | Interpolates KH (n=1) to HNC (n->infinity); balances stability and accuracy |

where `t* = -beta*u + h - c` is the renormalized indirect correlation function.

### Solvation Free Energy

**HNC closure:**
```
mu_ex,HNC = k_B T * sum_alpha rho_alpha^V * integral[ dr * (1/2 * (h_alpha^UV(r))^2 - c_alpha^UV(r) - 1/2 * h_alpha^UV(r) * c_alpha^UV(r)) ]
```

**KH closure:**
```
mu_ex,KH = k_B T * sum_alpha rho_alpha^V * integral[ dr * (1/2 * (h_alpha^UV(r))^2 * Theta(-h_alpha^UV(r)) - c_alpha^UV(r) - 1/2 * h_alpha^UV(r) * c_alpha^UV(r)) ]
```

**Gaussian fluctuation (GF) approximation:**
```
mu_ex,GF = k_B T * sum_alpha rho_alpha^V * integral[ dr * (-c_alpha^UV(r) - 1/2 * h_alpha^UV(r) * c_alpha^UV(r)) ]
```

The solvation free energy decomposes as: `Delta_G_sol = Delta_G_pol + Delta_G_dis + Delta_G_cav`

---

## Key Namelists / Input Files

### 1D-RISM: `&parameters` Namelist

| Variable | Default | Description |
|----------|---------|-------------|
| `theory` | `'DRISM'` | Theory: `DRISM` (dielectrically consistent) or `XRISM` (extended) |
| `closure` | `'KH'` | Closure: `KH`, `PSEn` (e.g. `PSE3`), `HNC`, `PY` |
| `entropicDecomp` | 1 | Solve temperature derivative equations: 0=no, 1=yes |
| `dr` | 0.025 | Grid spacing in real space (Angstrom) |
| `nr` | 16384 | Number of grid points (product of small primes 2,3,5) |
| `temperature` | 298.15 | Temperature in Kelvin |
| `dieps` | (required) | Dielectric constant of solvent |
| `nsp` | (required) | Number of species (molecules) |
| `smear` | 1.0 | Charge smear parameter (Angstrom) for long-range asymptotics |
| `adbcor` | 0.5 | Numeric parameter for DRISM |
| `mdiis_nvec` | 20 | Number of MDIIS vectors |
| `mdiis_del` | 0.3 | MDIIS step size |
| `mdiis_restart` | 10 | Restart threshold multiplier |
| `tolerance` | 1e-12 | Target residual tolerance |
| `maxstep` | 10000 | Maximum iterations |
| `extra_precision` | 1 | Use extra precision routines: 0=no, 1=yes |
| `progress` | 1 | Write residue to stdout every `progress` iterations |
| `ksave` | -1 | Save intermediate solutions every `ksave` steps; <=0 disables |
| `selftest` | 0 | Perform self-consistency check: 0=no, 1=yes |
| `rout` | 0 | Largest real-space separation in output files (Angstrom); 0=all |
| `kout` | 0 | Largest reciprocal-space separation in output; 0=all |

### 1D-RISM: `outlist` Characters

| Char | Output File | Description |
|------|-------------|-------------|
| `U` | `.uvv` | Solvent site-site potential U^VV(r) |
| `X` | `.xvv` | Site-site susceptibility chi^VV(k) -- required for 3D-RISM |
| `G` | `.gvv` | Pair distribution function G^VV(r) |
| `B` | `.bvv` | Bridge correction B^VV(r) |
| `T` | `.therm` | Thermodynamic properties |
| `E` | `.exnvv`, `.n00` | Running and total excess coordination numbers |
| `N` | `.nvv` | Running coordination numbers |
| `Q` | `.q00` | Total excess charge |
| `S` | `.svv` | Structure factor S^VV(k) |

### 1D-RISM: `&species` Namelist (per species)

| Variable | Default | Description |
|----------|---------|-------------|
| `density` | (required) | Density of the species |
| `units` | `'M'` | Units: `M`, `mM`, `'1/A^3'`, `'g/cm^3'`, `'kg/m^3'` |
| `model` | (required) | Path to `.mdl` file with solvent parameters |

### 3D-RISM in sander: `&rism` Namelist

| Variable | Default | Description |
|----------|---------|-------------|
| `irism` | 0 | Enable 3D-RISM (in `&cntrl`): 0=off, 1=on |
| `closure` | `'KH'` | Comma-separated list: `KH`, `HNC`, `PSEn` |
| `tolerance` | `1e-5` | List of max residual values (one/two/n tolerances) |
| `buffer` | 14 | Minimum distance (Angstrom) solute to box edge; <0 for fixed box |
| `grdspc` | `0.5,0.5,0.5` | Linear grid spacing (Angstrom) for variable box |
| `ng3` | (none) | Grid points `nx,ny,nz` for fixed box (used if buffer<0) |
| `solvbox` | (none) | Box size `lx,ly,lz` (Angstrom) for fixed box |
| `ljTolerance` | -1 | LJ cutoff accuracy; -1 = tolerance/10; 0 = no cutoff; >0 = value |
| `asympKSpaceTolerance` | -1 | Reciprocal space asymptotics cutoff accuracy |
| `asympcorr` | `.true.` | Use long-range asymptotic corrections for thermodynamics |
| `periodic` | (none) | Use periodic boundaries: `pme` or `ewald` |
| `solvcut` | (none) | LJ cutoff for periodic calculations |
| `centering` | 1/0 | Solute centering in box (see table below) |
| `zerofrc` | 1 | Zero net solvation force: 0=no, 1=yes |
| `apply_rism_force` | 1 | Calculate solvation forces: 0=no, 1=yes |
| `mdiis_nvec` | 5 | Number of MDIIS vectors |
| `mdiis_del` | 0.7 | MDIIS step size |
| `mdiis_restart` | 10 | MDIIS restart threshold |
| `mdiis_method` | 2 | MDIIS implementation: 0=original, 1=BLAS, 2=BLAS+memory |
| `maxstep` | 10000 | Maximum iterations per solution |
| `npropagate` | 5 | Previous solutions for initial guess (0-5) |
| `verbose` | 0 | Diagnostic detail: 0=none, 1=iterations, 2=per-iteration |
| `write_thermo` | 1 | Print solvation thermodynamics |
| `ntwrism` | 0 | Write solvent density grids every `ntwrism` steps |
| `molReconstruction` | 0 | Output molecular reconstruction (water only) |
| `volfmt` | `'mrc'` | Volumetric data format: `mrc`, `ccp4`, `dx`, `xyzv` |
| `progress` | 1 | Display progress every `progress` iterations |
| `polarDecomp` | 0 | Polar/non-polar decomposition of solvation free energy |
| `entropicDecomp` | 0 | Energy/entropy decomposition (requires .xvv v1.000+) |
| `gfCorrection` | 0 | Gaussian fluctuation excess chemical potential |
| `pcpluscorrection` | 0 | PC+/3D-RISM correction |
| `uccoeff` | `0,0,0,0` | UC correction coefficients: `a,b[,a1,b1]` |
| `rismnrespa` | 1 | r-RESPA MTS: `rismnrespa * dt` = solvation time step |
| `fcestride` | 0 | FCE MTS: full solutions every `fcestride * rismnrespa` steps |
| `fcenbasis` | 20 | Number of previous full solutions stored |
| `fcenbase` | 20 | Number of previous solutions used for extrapolation |
| `fcesort` | 0 | Sort basis vectors by distance: 0=no, 1=yes |
| `fcecrd` | 0 | Coordinate type for FCE: 0=absolute, 1=distance, 2=relative |
| `fceweigh` | 0 | Weighted coordinates: 0=no, 1=yes |
| `fceenormsw` | 0 | Balancing parameter for least squares (fcetrans=2) |
| `fcetrans` | 0 | Extrapolation method: 0=basic, 1=rotated, 2=ASFE, 4-6=GSFE |
| `fceifreq` | 1 | GSFE mapping list update frequency |
| `fcentfrcor` | 0 | Net force correction for GSFE |
| `treeDCF` | `.true.` | Treecode for DCF long-range asymptotics |
| `treeTCF` | `.true.` | Treecode for TCF long-range asymptotics |
| `treeCoulomb` | `.false.` | Treecode for Coulomb potential |
| `treeDCFMAC` | 0.1 | Treecode MAC for DCF |
| `treeTCFMAC` | 0.1 | Treecode MAC for TCF |
| `treeCoulombMAC` | 0.1 | Treecode MAC for Coulomb |
| `treeDCFOrder` | 2 | Taylor series order for DCF |
| `treeTCFOrder` | 2 | Taylor series order for TCF |
| `treeCoulombOrder` | 2 | Taylor series order for Coulomb |
| `treeDCFN0` | 500 | Max targets per leaf for DCF |
| `treeTCFN0` | 500 | Max targets per leaf for TCF |
| `treeCoulombN0` | 500 | Max targets per leaf for Coulomb |

### Centering Options

| Value | Description |
|-------|-------------|
| -4 | Center-of-geometry with grid rounding, first step only |
| -3 | Center-of-mass with grid rounding, first step only |
| -2 | Center-of-geometry, first step only |
| -1 | Center-of-mass, first step only |
| 0 | No centering |
| 1 | Center-of-mass, every step (default for MD, open boundaries) |
| 2 | Center-of-geometry, every step (default for minimization) |
| 3 | Center-of-mass with grid rounding |
| 4 | Center-of-geometry with grid rounding |

### rism3d.snglpnt Command-Line Options

| Option | Description |
|--------|-------------|
| `--pdb` | PDB file (required) |
| `--prmtop` | Topology file (required) |
| `--rst` | Restart file for coordinates |
| `--y`/`--traj` | Trajectory file (NetCDF or ASCII) |
| `--xvv` | Bulk solvent susceptibility file (required) |
| `--guv` | Root name for G^UV output |
| `--cuv` | Root name for C^UV output |
| `--huv` | Root name for H^UV output |
| `--uuv` | Root name for U^UV output |
| `--asymp` | Root name for long-range asymptotics |
| `--quv` | Root name for charge density distribution |
| `--chgdist` | Root name for charge distribution |
| `--exchem` | Root name for excess chemical potential distribution |
| `--solvene` | Root name for solvation energy distribution |
| `--entropy` | Root name for solvation entropy distribution |
| `--potUV` | Root name for solute-solvent potential energy distribution |
| `--closure` | Whitespace-separated closure list: `KH`, `HNC`, `PSEn` |
| `--buffer` | Minimum buffer distance (Angstrom) |
| `--grdspc` | Grid spacing `x,y,z` |
| `--ng` | Grid points `nx,ny,nz` (fixed box) |
| `--solvbox` | Box length `lx,ly,lz` (fixed box) |
| `--tolerance` | Residual tolerance(s) |
| `--ljTolerance` | LJ accuracy |
| `--asympKSpaceTolerance` | Reciprocal space asymptotics accuracy |
| `--periodic` | Periodic potential: `pme` or `ewald` |
| `--solvcut` | LJ cutoff for periodic |
| `--noasympcorr` | Turn off long-range asymptotics for thermo |
| `--centering` | Centering method (-4 to 4) |
| `--mdiis_del` | MDIIS step size |
| `--mdiis_nvec` | Number of MDIIS vectors |
| `--mdiis_restart` | MDIIS restart threshold |
| `--maxstep` | Maximum iterations |
| `--npropagate` | Number of previous solutions |
| `--polarDecomp` | Polar/non-polar decomposition |
| `--entropicDecomp` | Energy/entropy decomposition |
| `--gf` | Gaussian fluctuation correction |
| `--pc+` | PC+/3D-RISM correction |
| `--uccoeff` | UC correction coefficients |
| `--molReconstruct` | Molecular reconstruction (water) |
| `--volfmt` | Volume format: `mrc`, `ccp4`, `dx`, `xyzv` |
| `--verbose` | Verbosity level: 0, 1, 2 |
| `--treeDCF` | Treecode for DCF: 0=no, 1=yes |
| `--treeTCF` | Treecode for TCF: 0=no, 1=yes |
| `--treeCoulomb` | Treecode for Coulomb: 0=no, 1=yes |
| `--treeDCFMAC` | Treecode MAC for DCF |
| `--treeTCFMAC` | Treecode MAC for TCF |
| `--treeCoulombMAC` | Treecode MAC for Coulomb |
| `--treeDCFOrder` | Taylor order for DCF |
| `--treeTCFOrder` | Taylor order for TCF |
| `--treeCoulombOrder` | Taylor order for Coulomb |
| `--treeDCFN0` | Leaf size for DCF |
| `--treeTCFN0` | Leaf size for TCF |
| `--treeCoulombN0` | Leaf size for Coulomb |

---

## Common Workflows

### Workflow 1: Generate bulk solvent .xvv file

```bash
# 1. Create rism1d input file (e.g., water.inp)
cat > water.inp << 'EOF'
&PARAMETERS
    THEORY='DRISM', CLOSURE='KH',
    NR=16384, DR=0.025,
    OUTLIST='X', ROUT=384, KOUT=0,
    MDIIS_NVEC=20, MDIIS_DEL=0.3, TOLERANCE=1.e-12,
    KSAVE=-1, PROGRESS=1, MAXSTEP=10000,
    SMEAR=1, ADBCOR=0.5,
    TEMPERATURE=298.15, DIEPS=78.497, NSP=1
/
&SPECIES
    DENSITY=55.296d0,
    MODEL="path/to/SPC.mdl"
/
EOF

# 2. Run 1D-RISM
rism1d water

# 3. Output: water.xvv (used by 3D-RISM)
```

### Workflow 2: Single-point 3D-RISM solvation analysis

```bash
# Using rism3d.snglpnt
rism3d.snglpnt --pdb solute.pdb --prmtop solute.prmtop \
    --xvv water.xvv --closure "KH" --buffer 14 --grdspc 0.5,0.5,0.5 \
    --guv guv --cuv cuv --huv huv --entropicDecomp --polarDecomp
```

### Workflow 3: 3D-RISM MD with multiple time stepping

```bash
# mdin for sander
cat > md_rism.in << 'EOF'
&cntrl
    imin=0, nstlim=10000, dt=0.001,
    ntt=3, temp0=300, gamma_ln=20,
    ntb=0, cut=999.,
    ntpr=100, ntwx=1000, ntwr=10000,
    irism=1,
/
&rism
    rismnrespa=5,
    fcenbasis=10, fcestride=2, fcecrd=2,
    closure='KH', buffer=14, grdspc=0.5,0.5,0.5,
/
EOF

sander -O -i md_rism.in -o md.out -p prmtop -c inpcrd \
    -xvv water.xvv -r md.rst
```

### Workflow 4: Trajectory post-processing

```bash
cat > post.in << 'EOF'
&cntrl
    imin=5, maxcyc=1,
    ntb=0, cut=9999.,
    ntx=1, ntpr=1, ntwx=1,
    irism=1,
/
&rism
    tolerance=1e-4,
    apply_rism_force=0,
    npropagate=1,
/
EOF

sander -O -i post.in -o post.out -p prmtop -c inpcrd \
    -xvv water.xvv -y trajectory.nc
```

### Workflow 5: Closure bootstrapping for difficult convergence

```bash
# In sander &rism:
closure='KH','PSE3','HNC',
tolerance=1,1,1e-5,

# Or in rism3d.snglpnt:
--closure KH PSE3 HNC --tolerance 1 1 1e-5
```

---

## Thermodynamic Output Quantities

| Quantity | Description | Units |
|----------|-------------|-------|
| `rism_excessChemicalPotential` | Solvation free energy (closure) | kcal/mol |
| `rism_excessChemicalPotentialGF` | Solvation free energy (GF functional) | kcal/mol |
| `rism_excessChemicalPotentialPCPLUS` | Solvation free energy (PC+ correction) | kcal/mol |
| `rism_excessChemicalPotentialUC` | Solvation free energy (UC correction) | kcal/mol |
| `rism_solventPotentialEnergy` | Solute-solvent interaction energy | kcal/mol |
| `rism_excessParticlesCorrected` | Excess solvent particles vs uniform | # |
| `rism_excessChargeCorrected` | Excess solvent charge vs uniform | e |
| `rism_KirkwoodBuff` | All-space integral of h(r) | Angstrom^3 |
| `rism_DCFintegral` | All-space integral of c(r) | Angstrom^3 |
| `solutePotentialEnergy` | Total solute potential energy | kcal/mol |

---

## Solvation Free Energy Corrections

| Correction | Key Property | Notes |
|------------|-------------|-------|
| **UC** (Universal Correction) | `uccoeff=a,b[,a1,b1]` | Parameterized on PMV; a,b for closure; a1,b1 for temperature dependence |
| **PC+/3D-RISM** | `pcpluscorrection=1` | Pressure correction plus; based on PMV |
| **GF** (Gaussian Fluctuation) | `gfCorrection=1` | Not closure-specific; improved absolute solvation free energies |
| **NgB** (Ng Bridge) | Manual calculation | Results from polarDecomp + standard thermo output |

---

## Numerical Accuracy Guidelines

| Quantity | Recommended Value |
|----------|-------------------|
| Grid spacing | 0.3-0.5 Angstrom (0.5 is minimum for usable results) |
| Buffer | >= 14 Angstrom for water; larger for ionic solutions |
| MD tolerance | 1e-5 |
| Minimization tolerance | 1e-11 or lower |
| Trajectory post-processing tolerance | 1e-4 |
| Approximate error in Delta_G_solv | ~10 * tolerance |
| LJ tolerance | u_LJ(r_cut) <= tolerance/10 |
| 1D-RISM grid spacing | 0.025 Angstrom |

---

## Solvent Box Size: Interaction of buffer, ljTolerance, and tolerance

| buffer | ljTolerance < 0 | ljTolerance = 0 | ljTolerance > 0 |
|--------|-----------------|-----------------|-----------------|
| **< 0** (fixed) | Fixed box; LJ cutoff fit to box; correction applied | Fixed box; no LJ cutoff | Fixed box; LJ cutoff with ljTolerance; correction if box large enough |
| **= 0** | ljTolerance=tolerance/10; box fits cutoff | Error | Box fits cutoff; correction if large enough |
| **> 0** (variable) | Box from buffer; LJ cutoff fit to box; correction applied | Box from buffer; no LJ cutoff | Box from buffer; correction if large enough |

---

## Worked Example: Pure Water 1D-RISM for 3D-RISM

Generate an .xvv file for SPC/E water at 310 K:

```
&PARAMETERS
    THEORY='DRISM', CLOSURE='KH',
    NR=16384, DR=0.025,
    OUTLIST='X', ROUT=384, KOUT=0,
    MDIIS_NVEC=20, MDIIS_DEL=0.3, TOLERANCE=1.e-12,
    KSAVE=-1, PROGRESS=1, MAXSTEP=10000,
    SMEAR=1, ADBCOR=0.5,
    TEMPERATURE=310, DIEPS=78.497, NSP=1
/
&SPECIES
    DENSITY=55.296d0,
    MODEL="../../../dat/rism1d/model/SPC.mdl"
/
```

Run: `rism1d water` -> produces `water.xvv`

---

## Key Takeaways

1. **1D-RISM first**: Always run `rism1d` to generate the `.xvv` file before any 3D-RISM calculation. The `.xvv` file contains the bulk solvent susceptibility and is reusable for any solute.
2. **KH closure is default and most robust**: Use KH for molecular mechanics; use HNC or PSE-n for higher accuracy (with closure bootstrapping from KH).
3. **Grid spacing matters**: 0.5 Angstrom is the largest usable; 0.3-0.5 Angstrom is recommended. Buffer should be >= 14 Angstrom for water.
4. **Multiple time stepping saves computation**: Combine r-RESPA (`rismnrespa`) with FCE (`fcestride`) to compute 3D-RISM solutions only every 4 ps.
5. **Temperature derivatives require .xvv v1.000+**: Use `entropicDecomp` in rism1d to include temperature derivative information in the .xvv file.

## Connects To

- Chapter 23: sander molecular dynamics (MTS features)
- Chapter 6: PBSA implicit solvent
- Chapter 4: Generalized Born implicit solvent
- Chapter 37: ambpdb (PDB generation for rism3d.snglpnt)