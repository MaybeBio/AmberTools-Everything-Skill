# Chapter 34: Free Energies -- Complete Methodology

## Core Commands & Syntax

### Thermodynamic Integration (sander -- multisander mode)
```bash
mpirun -np 4 sander.MPI -ng 2 -groupfile groupfile
```
Groupfile format:
```
-O -i mdin -p prmtop.0 -c eq1.x -o md1.o -r md1.x -inf mdinfo
-O -i mdin -p prmtop.1 -c eq1.x -o md1b.o -r md1b.x -inf mdinfob
```

### Thermodynamic Integration (pmemd -- single prmtop)
```bash
pmemd.cuda -O -i ti.in -o ti.out -p ti.prmtop -c eq.rst7 -r ti.rst7 -x ti.nc
```

## Key Namelist Variables: &cntrl (sander TI)

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `icfe` | int | 0 | Enable TI: 0=off, 1=on. Uses mixing rule Eq.27.3 |
| `clambda` | float | 0.0 | λ value (0=V0, 1=V1) |
| `klambda` | int | 1 | Exponent k in mixing rule Eq.27.4. k=1 linear, k≥4 handles dummy atoms |
| `ifsc` | int | 0 | Softcore potential: 0=off, 1=on. Recommended over klambda>1 |
| `scalpha` | float | 0.5 | Softcore α parameter (Å²). Purely empirical |
| `scbeta` | float | 12.0 | Softcore β parameter. Purely empirical |
| `tishake` | int | 0/1 | SHAKE in TI: 0=keep constraints (sander default), 1=remove between common+softcore (pmemd default) |
| `idecomp` | int | 0 | Per-residue decomposition: 0=off, 1=decomp (1-4→internal), 2=decomp (1-4→EEL/VDW). Values 3,4 NOT compatible with icfe=1 |
| `dvdl_nout` | int | -- | Frequency to write ∂V/∂λ data |
| `logdvdl` | int | 0 | Write ∂V/∂λ to dvdl_dat every `dvdl_nout` steps |

### sander-specific TI flags

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `icfe` | int | 0 | Free energy calculation on/off (sander) |
| `timask1` | str | '' | Atoms unique to V0 (ambmask format) |
| `timask2` | str | '' | Atoms unique to V1 (ambmask format) |
| `scmask1` | str | '' | Softcore atoms for V0 |
| `scmask2` | str | '' | Softcore atoms for V1 |
| `crdgrow` | int | 0 | Growing atoms during TI (1=on) |
| `ifmols` | int | 0 | Allow molecule count to change during TI |
| `ntlambda` | int | -- | Number of λ windows |

### pmemd-specific TI flags

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `timask1` | str | '' | Atoms unique to V0 in ambmask format |
| `timask2` | str | '' | Atoms unique to V1 |
| `scmask1` | str | '' | Softcore region for V0 |
| `scmask2` | str | '' | Softcore region for V1 |
| `sccomp` | float | -- | Softcore composition parameter for pmemd |
| `gti_cpu_output` | int | 0 | CPU-compatible TI output format |
| `gti_add_sc` | int | 0 | Add softcore to output |

### pmemd.cuda-specific TI flags

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `icfe` | int | 0 | TI on/off |
| `ifsc` | int | 0 | Softcore: 0=off, 1=on, 2=smoothstep |
| `scalpha` | float | 0.5 | Softcore α |
| `scbeta` | float | 12.0 | Softcore β |
| `logdvdl` | int | 0 | Write ∂V/∂λ values |
| `dvdl_nout` | int | -- | Write frequency |
| `dvdl_dat` | str | '' | dvdl data file |
| `dvdl_flt` | str | '' | Float dvdl output file |
| `ntlambda` | int | -- | Number of λ values |
| `clambda` | float | -- | Current λ value |

## Background Theory

### TI Integral
```
ΔA = A(λ=1) − A(λ=0) = ∫₀¹ ⟨∂V/∂λ⟩_λ dλ
```

### Numerical Integration (Gaussian Quadrature)
```
ΔA = Σ wᵢ ⟨∂V/∂λ⟩_i
```

| Points | λ values | Weights |
|--------|----------|---------|
| 2 | 0.21132, 0.78868 | 0.5, 0.5 |
| 3 | 0.11270, 0.5, 0.88730 | 0.27778, 0.44444, 0.27778 |
| 5 | 0.04691, 0.23076, 0.5, 0.76923, 0.95308 | 0.11846, 0.23931, 0.28444, 0.23931, 0.11846 |
| 7 | 0.02545, 0.12923, 0.29708, 0.5, 0.70292, 0.87077, 0.97455 | 0.06468, 0.14048, 0.19000, 0.20938, 0.19000, 0.14048, 0.06468 |

### Mixing Rules

**Linear (klambda=1):**
```
V(λ) = (1−λ)V₀ + λV₁
⟨∂V/∂λ⟩ = V₁ − V₀
```

**Nonlinear (klambda=k):**
```
V(λ) = (1−λ)^k V₀ + [1−(1−λ)^k] V₁
```
k≥4 keeps integrand finite as λ→1 with dummy atoms.

## Softcore Potentials

### Traditional Softcore (ifsc=1, sander/pmemd)
Eliminates endpoint singularities when creating/annihilating atoms:
```
V_sc(r,λ) = 4ε(1−λ) [1/(αλ+(r/σ)^6)² − 1/(αλ+(r/σ)^6)]
```
where α=`scalpha`, the softcore α parameter (default 0.5).

### Smoothstep Softcore (ifsc=2, pmemd.cuda only)
Uses a smoothstep function for λ-scheduling, providing better numerical stability:
```
f(λ) = 1 − smoothstep(λ, sccomp, ...)
```
- `sccomp` controls the composition of softcore potential
- `sccompr` controls the radial composition

### Common TI Workflow
```
1. Charge removal (λ: 0→1, remove charges on disappearing atoms)
2. vdW removal    (λ: 0→1, remove vdW on disappearing atoms)  
3. vdW addition   (λ: 0→1, add vdW on appearing atoms)
4. Charge addition(λ: 0→1, add charges on appearing atoms)
```

## Performance Notes

| Engine | Performance vs MD | Notes |
|--------|-------------------|-------|
| sander (multisander) | ~50% (2 groups) | Two force calls per step |
| pmemd (CPU) PME TI | ~75% | Single reciprocal+dual direct |
| pmemd (CPU) GB TI | ~50% | GB radii not pairwise decomposable |
| pmemd.cuda TI | ~70% | GPU-accelerated |
| vdW-only softcore | ~100% | Charges identical, recip sum done once |

## GaMD + TI Integration

### GaMD Boosted TI (igmad=14)
Combines Gaussian accelerated MD with thermodynamic integration for enhanced sampling of TI windows.

```fortran
&cntrl
  icfe = 1, ifsc = 1, gti_cpu_output = 0, gti_add_sc = 1,
  timask1 = ':1-3', scmask1 = ':1-3',
  timask2 = '', scmask2 = '',
  igamd = 14, iE = 1, iEP = 1, iED = 1, irest_gamd = 0,
  ntcmd = 1000000, nteb = 1000000, ntave = 50000,
  ntcmdprep = 200000, ntebprep = 200000,
  sigma0P = 6.0, sigma0D = 6.0,
/
```

### Ligand GaMD TI (igmad=15)
For ligand-only boosting in protein-ligand TI.

### PPI-GaMD TI (igmad=17)
For protein-protein interaction boosting:
```fortran
bgpro2atm=869, edpro2atm=1736,
```

## Post-Processing

### alchemical_analysis.py
```bash
alchemical_analysis.py -d dvdl_data/ -o results/ -t 300 -u kcal
```

### BAR/MBAR Analysis
```bash
# Using pymbar
python -c "import pymbar; ..."
```

## Key Takeaways

1. **pmemd.cuda TI is ~70% of standard MD performance** -- use GPU for production TI
2. **Prefer softcore (ifsc=1) over klambda tricks** -- better convergence for atom creation/annihilation
3. **Gaussian quadrature with 5-7 windows is standard** -- more is better for rough transformations
4. **vdW-only softcore runs at full MD speed** -- charges identical, skip extra recip sum
5. **Always validate**: check restart files match (sander) and run 50-step test with ntpr=1
6. **pmemd uses single prmtop** with both V0 and V1 -- sander uses two prmtops + multisander

## Connects To

- **ch11**: Basic TI workflow with softcore potentials
- **ch12**: MMPBSA.py for endpoint free energies
- **ch14**: Enhanced sampling (GaMD, REMD) that can accelerate TI
- **ch21**: FEW (Free Energy Workflow) for automated TI setup
- **ch28**: BAR/PBSA post-processing analysis