# AmberTools Cheatsheet

## Force Field Selection

| System | Force Field | leaprc |
|--------|------------|--------|
| Protein (current) | FF19SB | `leaprc.protein.ff19SB` |
| Protein (legacy) | FF14SB | `leaprc.protein.ff14SB` |
| DNA | OL15/OL21 | `leaprc.DNA.OL15` |
| RNA | OL3 | `leaprc.RNA.OL3` |
| Small molecule | GAFF2 | via antechamber |
| Lipid | Lipid21 | `leaprc.lipid21` |
| Carbohydrate | GLYCAM_06j | `leaprc.GLYCAM_06j` |
| Unnatural AA | ff15ipq-m | `leaprc.ff15ipq` |
| Ionic liquid | GAFF2 | `leaprc.gaff2` |

## Water Model Selection

| Model | Type | Best With | Speed |
|-------|------|-----------|-------|
| OPC | 4-point | FF19SB | Moderate |
| TIP3P | 3-point | FF14SB, GAFF | Fastest |
| TIP4P-Ew | 4-point | FF14SB, OL15 | Moderate |
| OPC3 | 3-point | FF19SB | Fast |
| SPC/E | 3-point | General | Fast |
| TIP4P-D | 4-point | Disordered proteins | Moderate |

## MD Engine Selection

| Engine | When to Use |
|--------|------------|
| `pmemd.cuda` | Production MD on GPU (recommended) |
| `pmemd.cuda.MPI` | Multi-GPU production |
| `pmemd` | CPU-only, subset of sander features |
| `pmemd.MPI` | Multi-node CPU |
| `sander` | Full features: NMR, QM/MM, TI, LES |
| `sander.MPI` | Parallel sander for QM/MM, REMD |

## GB Model Selection (igb)

| igb | Model | Radii Set | Use When |
|-----|-------|-----------|----------|
| 1 | GB^OBC I | mbondi | Standard GB |
| 2 | GB^OBC II | mbondi2 | Proteins |
| 5 | GB^OBC II | mbondi2 | Standard (FF14SB) |
| 7 | GB^neck2 | bondi | Nucleic acids |
| 8 | GB^neck2 | mbondi3 | Recommended (FF19SB) |

## Charge Method Selection

| Method | Command | Accuracy | Speed |
|--------|---------|----------|-------|
| AM1-BCC | `antechamber -c bcc` | Good | Fast |
| RESP | `antechamber -c resp` | Best | Slow (needs Gaussian) |
| Gasteiger | `antechamber -c gas` | Rough | Instant |
| Mulliken | `antechamber -c mul` | Poor | Via sqm |

## Minimization Protocol Defaults

| Parameter | Stage 1 (backbone restraint) | Stage 2 (all atom) |
|-----------|------------------------------|---------------------|
| imin | 1 | 1 |
| ncyc | 500 | 500 |
| maxcyc | 1000 | 1000 |
| ntr | 1 | 0 |
| restraint_wt | 10.0 | -- |
| restraintmask | `:1-300@CA,C,N` | -- |
| ntc | 1 | 2 |
| ntf | 1 | 2 |

## Standard MD Timestep

| dt | SHAKE | HMR Required | Use Case |
|----|-------|-------------|----------|
| 0.001 (1 fs) | Optional | No | QM/MM, high-T, flexible water |
| 0.002 (2 fs) | ntc=2, ntf=2 | No | **Standard production** |
| 0.004 (4 fs) | ntc=2, ntf=2 | Yes | Fast production with HMR |

## Cutoff Defaults

| Method | cut | Ewald/PME |
|--------|-----|-----------|
| PME (explicit water) | 8.0 A | Auto |
| GB (igb>0) | 999.0 (infinite) | N/A |
| Gas phase (igb=0, ntb=0) | 999.0 | N/A |
| PBSA | 999.0 | N/A |

## Temperature Control

| ntt | Method | gamma_ln | When |
|-----|--------|----------|------|
| 1 | Berendsen | -- | Equilibration only |
| 2 | Andersen | -- | NVE checks |
| 3 | Langevin | 1.0-5.0 | **Production (default)** |
| 9 | Bussi | -- | Modern alternative |

## SHAKE Defaults

```fortran
ntc=2, ntf=2, tol=0.00001
```

## Box Types

| Command | Shape | Water Savings |
|---------|-------|--------------|
| `solvateoct` | Truncated octahedron | ~25% |
| `solvatebox` | Cubic | 0% (baseline) |

## Common Ion Parameters

| Ion | Residue Name | LEaP Command |
|-----|-------------|-------------|
| Na+ | Na+ | `addions2 mol Na+ 0` |
| Cl- | Cl- | `addions2 mol Cl- 0` |
| K+ | K+ | `addions2 mol K+ 0` |
| Mg2+ | Mg2+ | `addions2 mol Mg2+ 0` |
| Ca2+ | Ca2+ | `addions2 mol Ca2+ 0` |

## Trajectory Output

| Parameter | Value | Notes |
|-----------|-------|-------|
| ioutfm | 1 | NetCDF (recommended) |
| ioutfm | 0 | ASCII mdcrd (legacy) |
| ntwx | 5000 | Write every 5000 steps = 10 ps at dt=0.002 |
| ntwr | 5000 | Restart file frequency |
| ntpr | 500 | Print to mdout frequency |
| ntwv | -1 | Velocity output (off for production) |
| ntwe | 5000 | Energy output to mden |

## Common Restraint Setup

| Type | ntr | restraint_wt | restraintmask |
|------|-----|-------------|---------------|
| Backbone | 1 | 10.0 | `:1-300@CA,C,N` |
| Heavy atom | 1 | 5.0 | `:1-300&!@H=` |
| Distance | nmropt=1 | via DISANG | rk2=20.0 |

## Typical nstlim Values

| Stage | nstlim | Time at dt=0.002 |
|-------|--------|------------------|
| Minimization | 500-2500 | N/A |
| Heating | 10000-50000 | 20-100 ps |
| Equilibration | 50000-250000 | 100-500 ps |
| Production | 500000-5000000 | 1-10 ns |

## CPPTRAJ Atom Mask Quick Reference

| Mask | Selects |
|------|---------|
| `:1-10@CA` | Residues 1-10, C-alpha atoms |
| `@/H` | Strip all hydrogens |
| `:WAT` | All water molecules |
| `:Na+,Cl-` | Sodium and chloride ions |
| `:1-300@CA,C,N` | Protein backbone |
| `:LIG` | Residue named LIG |
| `!@H=` | All non-hydrogen (heavy) atoms |
| `:POPC` | POPC lipid residues |
| `@CA` | All C-alpha atoms |

## TI Window Defaults

| Parameter | Value |
|-----------|-------|
| clambda | 0.0 to 1.0 |
| nlambda | 12-24 windows |
| icfe | 1 |
| ifsc | 1 (softcore) |
| scalpha | 0.5 |
| scbeta | 12.0 |

## MMPBSA.py igb Defaults

| Parameter | Default |
|-----------|---------|
| igb | 5 |
| saltcon | 0.15 |
| startframe | 1 |
| endframe | 9999999 |
| interval | 1 |

## REMD Temperature Distribution

```python
# Geometric distribution: T_i = T_min * (T_max/T_min)^(i/(N-1))
# Typical: 32 replicas, 300 K to 450 K
```

## antechamber -c Method Comparison

| -c value | Method | Requires |
|----------|--------|----------|
| bcc | AM1-BCC | sqm |
| resp | RESP | Gaussian03+ |
| gas | Gasteiger | Nothing |
| mul | Mulliken | sqm |
| rc | Read from mol2 | mol2 with charges |
| dc | Delete charges | Nothing |

## Input File Naming Conventions

| Extension | File Type | Program |
|-----------|-----------|---------|
| .prmtop / .parm7 | Topology | LEaP output |
| .inpcrd / .rst7 | Coordinates | LEaP output |
| .mdin | MD input | sander/pmemd input |
| .mdout | MD output log | sander/pmemd output |
| .mdcrd | ASCII trajectory | MD output |
| .nc | NetCDF trajectory | MD output (ioutfm=1) |
| .mden | Energy file | MD output |
| .mdvel | Velocity file | MD output |
| .mdinfo | MD info | MD output |
| .mol2 | Small molecule | antechamber |
| .frcmod | Force field mod | parmchk2 |
| .off | OFF library | LEaP |
| .prep / .prepi | PREP residue | antechamber |
| .lib | Library file | LEaP |
| .leaprc | LEaP resource | LEaP |
| .pdb | Protein structure | pdb4amber |
| .f | DISANG restraints | sander (NMR) |
| .dx | Volumetric map | 3D-RISM, GIST |
| .ccp4 / .mrc | EM density map | EMAP |
| .xvv | Solvent susceptibility | rism1d |
| .yaml | YAML config | BAR/PBSA |
| .json | JSON config | KMMD, FRETrest |
| .xml | XML config | edgembar |