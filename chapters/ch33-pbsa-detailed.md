# Chapter 32: PBSA -- Poisson-Boltzmann Surface Area

## Core Commands & Syntax

```bash
pbsa [-O] -i mdin -o mdout [-p prmtop -c inpcrd]/[-pqr pqr]
```

PBSA performs single-structure implicit solvent calculations using the Poisson-Boltzmann equation. It supports both Amber topology files and free-format PQR files.

### GPU version

```bash
pbsa.cuda -O -i mdin -o mdout -p prmtop -c inpcrd
```

### In sander (for MD/minimization)

PBSA functionality is available in sander by setting `ipb` to nonzero. All `&pb` namelist options are identical to standalone pbsa.

---

## The Poisson-Boltzmann Equation

```
nabla * [epsilon(r) * nabla * phi(r)] = -4*pi*rho(r) - 4*pi*lambda(r) * sum_i[ z_i * c_i * exp(-z_i * phi(r) / k_B * T) ]
```

where:
- `epsilon(r)` = dielectric constant
- `phi(r)` = electrostatic potential
- `rho(r)` = solute charge distribution
- `lambda(r)` = Stern layer masking function
- `z_i` = charge of ion type i
- `c_i` = bulk number density of ion type i

### Numerical Methods

- **Finite-difference (FD) method**: Default. Maps charges to grid, assigns boundary conditions, applies dielectric model.
- **Immersed Interface Method (IIM)**: Second-order accurate for linear PB. `ipb=4`.
- **Analytical IIM**: `ipb=6` (recommended for GPU). Uses analytical surface routines.
- **X-factor harmonic average**: `ipb=7`.
- **Second-order harmonic average**: `ipb=8`.

### Linear Solvers

| solvopt | Solver | Notes |
|---------|--------|-------|
| 1 | Modified ICCG (or PICCG for PBC) | Default |
| 2 | Geometric multigrid | 4-level v-cycle; fastest for large systems |
| 3 | Conjugate gradient | Requires large `maxitn` |
| 4 | SOR | Requires large `maxitn` |
| 5 | Adaptive SOR | Nonlinear only |
| 6 | Damped SOR | Nonlinear only |

### Nonlinear PB

Set `npbopt=1` for the full nonlinear PB equation. Default is `npbopt=0` (linear).

---

## Key Namelists / Input Files

### `&cntrl` Namelist (basic controls)

| Variable | Default | Description |
|----------|---------|-------------|
| `imin` | 1 | Minimization flag: 0=no, 1=single point (default in pbsa) |
| `ntx` | 1 | Coordinate read format: 1=formatted, 2=unformatted |
| `ipb` | 2 | Dielectric model: 0=off, 1=geometric, 2=level set (default), 4=IIM, 6=analytical IIM, 7=X-factor HA, 8=2nd-order HA |
| `inp` | 2 | Non-polar model: 0=off, 1=total NP (SASA), 2=cavity+dispersion (default) |

### `&pb` Namelist: Physical Constants

| Variable | Default | Description |
|----------|---------|-------------|
| `epsin` | 1.0 | Solute dielectric constant |
| `epsout` | 80.0 | Solvent dielectric constant |
| `epsmem` | 1.0 | Membrane dielectric constant (used if `membraneopt>0`) |
| `smoothopt` | 1 | Dielectric boundary smoothing: 0=equal-weight harmonic, 1=weighted harmonic (default), 2=midpoint assignment |
| `istrng` | 0 | Ionic strength in mM (note: different from GB which uses M) |
| `pbtemp` | 300 | Temperature in K for PB equation |
| `radiopt` | 1 | Atomic radii: 0=from prmtop, 1=Tan-Luo optimized (default) |
| `dprob` | 1.4 | Solvent probe radius for molecular surface (Angstrom) |
| `iprob` | 2.0 | Mobile ion probe radius for Stern layer (Angstrom) |
| `sasopt` | 0 | Molecular surface type: 0=SES, 1=SAS, 2=smooth density function, 3=MLSES |
| `saopt` | 0 | Surface area computation: 0=no, 1=field-view method |
| `triopt` | 1 | Trimer arc dots: 0=off, 1=on (default) |
| `arcres` | 0.25 | Arc dot resolution (Angstrom) |

### `&pb` Namelist: Implicit Membrane

| Variable | Default | Description |
|----------|---------|-------------|
| `membraneopt` | 0 | Membrane model: 0=off, 1=uniform slab, 2=heterogeneous (PCHIP), 3=heterogeneous (Spline) |
| `mprob` | 2.70 | Membrane probe radius (Angstrom) |
| `mthick` | 40.0 | Membrane thickness (Angstrom) |
| `mctrdz` | 0 | Membrane center in z-direction (Angstrom) |
| `poretype` | 0 | Pore searching: 0=off, 1=on |
| `poreradius` | (none) | Cylindrical exclusion radius (Angstrom) |

### `&pb` Namelist: Numerical Procedures

| Variable | Default | Description |
|----------|---------|-------------|
| `npbopt` | 0 | PB equation: 0=linear, 1=nonlinear |
| `solvopt` | 1 | Iterative solver: 1=ICCG, 2=multigrid, 3=CG, 4=SOR, 5=adaptive SOR, 6=damped SOR |
| `accept` | 0.001 | Iteration convergence criterion (relative to initial residue) |
| `maxitn` | 100 | Maximum iterations for FD solver |
| `fillratio` | 2.0 | Ratio of grid dimension to solute dimension |
| `space` | 0.5 | Grid spacing for FD solver (Angstrom) |
| `nbuffer` | 0 | Grid units between solute surface and FD grid boundary |
| `nfocus` | 2 | Number of successive FD calculations for focusing |
| `fscale` | 8 | Ratio between coarse and fine grid spacings |
| `npbgrid` | 1 | Grid regeneration frequency (for MD: >= 100) |

### `&pb` Namelist: Boundary Conditions

| Variable | Default | Description |
|----------|---------|-------------|
| `bcopt` | 5 | Boundary condition: 1=conductor, 5=all grid charges (default), 6=charge singularity free, 10=periodic |
| `eneopt` | 2 | Electrostatic energy method: 1=P3M (EPB=0, EEL=total), 2=surface charge (default, EPB=reaction field), 3=P3M for large systems, 4=P3M for nonlinear |
| `frcopt` | 0 | Force output: 0=no, 1=trilinear interpolation, 2=surface polarized charges, 3=surface charge + boundary field |
| `scalec` | 0 | Scale dielectric boundary charges: 0=no, 1=Gauss's law |
| `cutfd` | 5 | Atom-based cutoff for short-range FD interactions (Angstrom) |
| `cutnb` | 0 | Cutoff for vdW and pairwise Coulombic (Angstrom); 0=infinite |
| `nsnba` | 1 | Pairlist regeneration frequency (for MD: 5) |

### `&pb` Namelist: Non-Polar Solvation

| Variable | Default | Description |
|----------|---------|-------------|
| `decompopt` | 2 | Decomposition scheme: 1=6/12, 2=sigma (default), 3=WCA |
| `use_rmin` | 1 | vdW radii: 0=sigma, 1=rmin (default) |
| `sprob` | 0.557 | Solvent probe radius for SASA in dispersion (Angstrom) |
| `vprob` | 1.300 | Solvent probe radius for molecular volume (Angstrom) |
| `rhow_effect` | 1.129 | Effective water density for dispersion |
| `use_sav` | 1 | Cavity term: 0=SASA, 1=volume enclosed by SASA (default) |
| `cavity_surften` | (optimized) | Regression coefficient for cavity/SASA relation |
| `cavity_offset` | (optimized) | Regression offset for cavity/SASA relation |
| `maxsph` | 400 | Approximate dots for max atomic SASA |

### `&pb` Namelist: Visualization

| Variable | Default | Description |
|----------|---------|-------------|
| `phiout` | 0 | Write electrostatic potential: 0=no, 1=yes |
| `phiform` | 0 | Potential format: 0=Delphi binary, 1=Amber ASCII, 2=DX |
| `outlvlset` | `false` | Write total level set to DX file |
| `outmlvlset` | `false` | Write membrane level set to DX file |
| `npbverb` | 0 | Verbose output: 0=no, 1=yes |

### `&pb` Namelist: Active Site Focusing

| Variable | Default | Description |
|----------|---------|-------------|
| `xmin` | 0 | Lower x boundary of local region |
| `xmax` | 0 | Upper x boundary of local region |
| `ymin` | 0 | Lower y boundary |
| `ymax` | 0 | Upper y boundary |
| `zmin` | 0 | Lower z boundary |
| `zmax` | 0 | Upper z boundary |

### `&pb` Namelist: MLSES (Machine-Learned SES)

| Variable | Default | Description |
|----------|---------|-------------|
| `mlses_opt` | 0 | MLSES runtime: 0=Fortran/CUDA, 1=Torch GENIUSES, 2=Torch Con2SES-2D, 3=Torch Con2SES-3D |

---

## Energy Components in Output

| Component | Description |
|-----------|-------------|
| `EEL` | Total electrostatic energy (Coulombic + reaction field, or Coulombic only depending on `eneopt`) |
| `EPB` | Reaction field energy (when `eneopt=2`); zero when `eneopt=1` |
| `ECAVITY` | Cavity solvation free energy (or total NP when `inp=1`) |
| `EDISPER` | Dispersion solvation free energy (zero when `inp=1`) |
| `VDWAALS` | van der Waals energy |

The total solvation free energy decomposition: `Delta_G_sol = Delta_G_pol + Delta_G_dis + Delta_G_cav`
where `Delta_G_pol = EPB` (or `EEL` minus Coulombic), `Delta_G_dis = EDISPER`, `Delta_G_cav = ECAVITY`.

---

## PQR File Format

```
Tag AtomNumber AtomName ResidueName ChainID ResidueNumber X Y Z Charge Radius
```

- `Tag`: `ATOM` or `HETATM`
- Fields are space-delimited
- Column order must be: Charge then Radius (P.Q.R. order)

---

## Common Workflows

### Workflow 1: Single-point PB calculation with defaults

```
&cntrl
/
&pb
    npbverb=1, istrng=150, fillratio=1.5, saopt=1,
/
```

```bash
pbsa -O -i mdin -o mdout -p prmtop -c inpcrd
```

### Workflow 2: Implicit membrane calculation

```
&cntrl
    ipb=1, inp=0
/
&pb
    radiopt=0, nfocus=1, maxitn=200,
    bcopt=10, eneopt=1, solvopt=1,
    sasopt=0, membraneopt=1, epsmem=4.0,
    outlvlset=true, outmlvlset=true,
/
```

### Workflow 3: Force computation

```
&cntrl
    inp=0
/
&pb
    npbverb=1, radiopt=0, frcopt=2,
/
```

### Workflow 4: PB Molecular Dynamics in sander

```
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

### Workflow 5: Delphi comparison

```
&cntrl
    ipb=1, inp=0
/
&pb
    istrng=150, ivalence=1, iprob=2.0, dprob=1.5,
    radiopt=0, bcopt=5, smoothopt=2, nfocus=1,
/
```

### Workflow 6: GPU-accelerated with multigrid

```
&cntrl
    ntx=1, imin=1, ipb=2, inp=0
/
&pb
    npbverb=1, istrng=0, epsout=80.0, epsin=1.0, space=.5,
    accept=0.0001, dprob=1.4, radiopt=1, fillratio=1.5,
    smoothopt=0, arcres=0.0625, nfocus=1,
    bcopt=1, solvopt=2, maxitn=3000
/
```

### Workflow 7: GPU with analytical IIM

```
&cntrl
    ntx=1, imin=1, ipb=6, inp=0
/
&pb
    npbverb=0, istrng=0, epsout=80.0, epsin=1.0, space=.5,
    accept=0.0001, dprob=1.4, radiopt=0, fillratio=1.5,
    smoothopt=0, nfocus=1, sasopt=2,
    bcopt=2, maxitn=3000, cutnb=15, cutsa=8, cutfd=7
/
```

### Workflow 8: MLSES surface

```
&cntrl
    ntx=1, imin=1, ipb=2, inp=0
/
&pb
    npbverb=1, istrng=0,
    epsout=80.0, epsin=1.0, dprob=1.4, radiopt=0, sasopt=3,
    fillratio=1.25, nfocus=1, space=0.5,
    accept=0.000001, maxitn=10000, solvopt=3,
    npbopt=0, bcopt=6,
    eneopt=1, frcopt=0, cutnb=15, cutsa=8, cutfd=7
/
```

### Workflow 9: Visualization with PyMol

```
&cntrl
    inp=0
/
&pb
    npbverb=1, space=1.,
    phiout=1, phiform=0
/
```

In PyMol:
```python
load molecule.top
load molecule.rst
# Display surface via GUI
load pbsa.phi
ramp_new e_lvl, pbsa, [-7, 0, 7]
set surface_color, e_lvl, molecule
```

### Workflow 10: Visualization with VMD (DX format)

```
&cntrl
    inp=0
/
&pb
    npbverb=1, space=1., sasopt=2,
    phiout=1, phiform=2
/
```

---

## GPU Implementation Details

| Aspect | Detail |
|--------|--------|
| Precision | Hybrid: single for linear solve, double for setup/post-processing |
| Boundary conditions | NBC (`bcopt=1` or `5`) and PBC (`bcopt=10`) |
| MG solver memory | ~75 * N_grid bytes |
| MG-Jacobi-PCG hybrid memory | ~135 * N_grid bytes |
| Jacobi-PCG memory | ~92 * N_grid bytes |
| Supported GPUs | NVIDIA CUDA (single precision) |

---

## Accuracy Recommendations

| Parameter | Recommended Value |
|-----------|-------------------|
| Grid spacing (single point) | 0.5 Angstrom |
| Grid spacing (MD/minimization) | 0.25 Angstrom |
| Arc resolution (MD) | 0.125 Angstrom |
| Focus scale (MD) | 4 |
| Reaction field convergence | ~2% at 0.5 Angstrom with `eneopt=2` |
| `fillratio` for small molecules | 4.0 (default 2.0 may be too small) |
| `maxitn` for CG/SOR | 10,000 (default 100 is too small) |
| `arcres` with TRIOPT | max(0.125, 0.5*h) |
| `arcres` without TRIOPT | max(0.0625, 0.25*h) |

---

## Worked Example: Single-Point Solvation Free Energy

Compute the solvation free energy of a protein at physiological salt concentration:

```
Sample single point PB calculation
&cntrl
/
&pb
    npbverb=1, istrng=150, fillratio=1.5, saopt=1,
/
```

Output quantities:
- `EPB` = reaction field (polar) free energy
- `ECAVITY` = cavity free energy (from SASA or volume)
- `EDISPER` = dispersion free energy
- `EEL` = Coulombic energy (or total electrostatic depending on `eneopt`)

The total solvation free energy = `EPB + ECAVITY + EDISPER`.

---

## Key Takeaways

1. **Use `eneopt=2` for better convergence**: The surface charge method gives reaction field energies that converge to ~2% at 0.5 Angstrom grid spacing.
2. **Increase `fillratio` for small molecules**: The default of 2.0 may cause solute to lie outside the grid; use 4.0 for ligands.
3. **GPU multigrid is fastest for large systems**: Use `solvopt=2` with `bcopt=1` (conductor boundary) for best performance.
4. **Use `inp=2` (default) for accurate NP**: Separates cavity and dispersion terms; the sigma decomposition scheme (`decompopt=2`) is recommended.
5. **For membrane proteins**: Use `membraneopt=1`, `sasopt=0` (SES), periodic boundary conditions (`bcopt=10`), and `eneopt=1`.

## Connects To

- Chapter 4: Generalized Born implicit solvent
- Chapter 29: RISM (3D-RISM implicit solvent)
- Chapter 23: sander MD (PB dynamics)
- Chapter 15: parmed (topology manipulation)
- Chapter 37: ambpdb (PDB generation)