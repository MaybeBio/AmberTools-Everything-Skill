# Chapter 8: Minimization, Relaxation & Equilibration

## Core Commands & Syntax

### Minimization Run
```bash
pmemd.cuda -O -i min.in -o min.out -p system.prmtop -c system.inpcrd \
    -r min.rst7 -ref system.inpcrd -inf min.info
```
The `-ref` flag is required when `ntr=1` (positional restraints). It provides the reference coordinates for the harmonic restraint potential.

### Heating Run (NVT, temperature ramp)
```bash
pmemd.cuda -O -i heat.in -o heat.out -p system.prmtop -c min.rst7 \
    -ref min.rst7 -r heat.rst7 -x heat.nc -inf heat.info
```

### Equilibration Run (NPT, restart from heating)
```bash
pmemd.cuda -O -i eq.in -o eq.out -p system.prmtop -c heat.rst7 \
    -r eq.rst7 -x eq.nc -inf eq.info
```

### Post-Processing Trajectory Energies
```bash
sander -O -i post.in -o post.out -p system.prmtop -c system.inpcrd \
    -y trajectory.nc -r post.rst7
```
Set `imin=5, maxcyc=1` for single-point energy evaluation of each trajectory frame.

## Key Namelist Variables

### Minimization Control (imin=1)
| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `imin` | int | 0 | 0=MD, 1=minimization, 5=trajectory analysis, 6=MD analysis per frame |
| `maxcyc` | int | 1 | Maximum number of minimization cycles |
| `ncyc` | int | 10 | Steepest descent cycles before switching to conjugate gradient (when ntmin=1) |
| `ntmin` | int | 1 | 0=CG only (4 SD at start), 1=SD then CG (default), 2=SD only, 3=XMIN, 4=LMOD, 5=DL-Find |
| `drms` | float | 1e-4 | RMS gradient convergence criterion (kcal/mol/Angstrom). Minimization halts when RMS gradient < drms |
| `dx0` | float | 0.01 | Initial step length. The minimizer auto-adjusts if too large |

### Restraints and Frozen Atoms
| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ntr` | int | 0 | 1=harmonic positional restraints on atoms in `restraintmask` |
| `restraintmask` | str | '' | Amber mask selecting restrained atoms (e.g., `':1-81'`, `'@CA,N,C'`) |
| `restraint_wt` | float | 0.0 | Force constant in kcal/mol/Ang^2. Potential: k*(delta_x)^2 per Cartesian dimension |
| `ibelly` | int | 0 | 1=belly dynamics (freeze atoms not in `bellymask`). Not available with igb>0 |
| `bellymask` | str | '' | Amber mask for moving atoms when ibelly=1 |
| `nmropt` | int | 0 | 1=read `&wt` namelist blocks for time-varying parameters (temperature, restraints) |

### Temperature Control (for Heating/Equilibration)
| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ntt` | int | 0 | 0=NVE, 1=Berendsen, 2=Andersen, 3=Langevin (recommended), 9=OIN, 10=SINR |
| `temp0` | float | 300.0 | Target temperature (K) |
| `tempi` | float | 0.0 | Initial temperature (K). Velocities sampled from Maxwellian at tempi if ntx=1 and tempi>0 |
| `gamma_ln` | float | 0.0 | Langevin collision frequency (ps^{-1}). Typical: 1.0-5.0 for explicit, 1.0 for GB |
| `ig` | int | -1 | Random seed. -1 = from system clock (recommended for Langevin) |
| `tautp` | float | 1.0 | Heat bath coupling time constant (ps) for ntt=1. Range 0.5-5.0 |

### Pressure Control and Periodic Boundaries
| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ntb` | int | auto | 0=no PBC (GB), 1=constant volume (NVT), 2=constant pressure (NPT). Auto-detected from igb/ntp |
| `ntp` | int | 0 | 0=no pressure scaling, 1=isotropic, 2=anisotropic, 3=semi-isotropic |
| `pres0` | float | 1.0 | Target pressure (bar). 1 bar ~ 0.987 atm |
| `taup` | float | 1.0 | Pressure relaxation time (ps). Recommended 1.0-5.0 |
| `barostat` | int | 1 | 1=Berendsen, 2=Monte Carlo (recommended for NPT sampling) |
| `mcbarint` | int | 100 | Steps between MC barostat volume change attempts |

### SHAKE and Timestep
| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ntc` | int | 1 | 1=no SHAKE, 2=bonds with H constrained, 3=all bonds |
| `ntf` | int | 1 | 1=all forces, 2=omit bond H forces, 3=omit all bond forces. Typically ntc=ntf |
| `dt` | float | 0.001 | Timestep (ps). 0.001=1 fs, 0.002=2 fs with SHAKE, 0.004=4 fs with HMR |
| `tol` | float | 1e-5 | SHAKE tolerance (Angstrom). Max recommended: 0.00005 |

### &wt Namelist (Time-Varying Parameters)
Used when `nmropt=1`. Each `&wt` block specifies a linear ramp between two step values.
```
&wt type='TEMP0', istep1=0, istep2=9000, value1=0.0, value2=300.0 /
&wt type='TEMP0', istep1=9001, istep2=10000, value1=300.0, value2=300.0 /
&wt type='END' /
```
| Parameter | Description |
|-----------|-------------|
| `type` | Variable to change: `TEMP0` (temperature), `RESTRAINT` (restraint weight) |
| `istep1` | Start step for the linear ramp |
| `istep2` | End step for the linear ramp |
| `value1` | Parameter value at `istep1` |
| `value2` | Parameter value at `istep2` |

## Common Workflows - Step-by-step Relaxation Protocol

### Explicit Solvent Relaxation Protocol (Standard 9-Step)

This protocol, from the official Amber tutorial, gradually relaxes the system before production MD. Each step is 1 ns at 1 fs timestep with restraints on the solute.

**Step 1: Restrained Minimization (solute heavy)**
```
Step 1 - Minimize with solute restraints
 &cntrl
  imin=1, maxcyc=1000, ncyc=20, ntx=1,
  ntpr=50, ntwr=500, ntc=2, ntf=2,
  ntb=1, ntp=0, cut=10.0,
  ntr=1, restraintmask=':1-81', restraint_wt=100.0,
  ioutfm=1, ntxo=2,
 /
```
Key: 100 kcal/mol/Ang^2 restraints on solute. Steepest descent for 20 cycles, then conjugate gradient.

**Step 2: Heat (NVT, 100K to 298K, with solute restraints)**
```
Step 2 - Heat from 100K to 298K
 &cntrl
  imin=0, nstlim=1000000, dt=0.001,
  irest=0, ntx=1, ig=-1,
  tempi=100.0, temp0=298.0,
  ntc=2, ntf=2, tol=0.00001,
  ntwx=10000, ntwr=1000, ntpr=1000,
  cut=8.0, iwrap=0,
  ntt=3, gamma_ln=1.0, ntb=1, ntp=0,
  nscm=0,
  ntr=1, restraintmask=':1-81', restraint_wt=100.0,
  nmropt=1,
  ioutfm=1, ntxo=2,
 /
&wt TYPE='TEMP0', istep1=0, istep2=1000000, value1=100.0, value2=298.0 /
&wt TYPE='END' /
```
Key: Start at 100K (not 0K), ramp linearly to 298K over 1M steps. 1 fs timestep during heating.

**Step 3: NPT Equilibration (solute heavy restraint)**
```
Step 3 - NPT with solute restraints
 &cntrl
  imin=0, nstlim=1000000, dt=0.001,
  irest=1, ntx=5, ig=-1,
  temp0=298.0,
  ntc=2, ntf=2, tol=0.00001,
  ntwx=10000, ntwr=1000, ntpr=1000,
  cut=8.0, iwrap=0,
  ntt=3, gamma_ln=1.0, ntb=2, ntp=1, barostat=2,
  nscm=0,
  ntr=1, restraintmask=':1-81', restraint_wt=100.0,
  ioutfm=1, ntxo=2,
 /
```
Key: Switch to NPT (ntb=2, ntp=1). Monte Carlo barostat (barostat=2). Restart from Step 2.

**Step 4: NPT, reduced restraint (10 kcal/mol)**
Same as Step 3 but `restraint_wt=10.0`.

**Step 5: Minimization with backbone restraints only**
```
Step 5 - Minimize with backbone restraints
 &cntrl
  imin=1, maxcyc=1000, ncyc=30, ntx=1,
  ntwr=500, ntpr=50,
  ntc=2, ntf=2, ntb=1, ntp=0,
  cut=8.0,
  ntr=1, restraintmask='@CA,N,C', restraint_wt=10.0,
  ioutfm=1, ntxo=2,
 /
```
Key: Restraint mask changes from `':1-81'` (all solute atoms) to `'@CA,N,C'` (backbone only).

**Step 6: NPT, backbone restraint (10 kcal/mol)**
Same as Step 3/4 but `restraintmask='@CA,N,C'`, `restraint_wt=10.0`.

**Step 7: NPT, backbone restraint reduced (1 kcal/mol)**
`restraintmask='@CA,N,C'`, `restraint_wt=1.0`.

**Step 8: NPT, backbone restraint further reduced (0.1 kcal/mol)**
`restraintmask='@CA,N,C'`, `restraint_wt=0.1`.

**Step 9: NPT, no restraints**
```
Step 9 - NPT, unrestrained
 &cntrl
  imin=0, nstlim=1000000, dt=0.001,
  irest=1, ntx=5, ig=-1,
  temp0=298.0,
  ntc=2, ntf=2, tol=0.00001,
  ntwx=10000, ntwr=1000, ntpr=1000,
  cut=8.0, iwrap=0,
  ntt=3, gamma_ln=1.0, ntb=2, ntp=1,
  nscm=1000, barostat=2,
  ioutfm=1, ntxo=2,
 /
```
Key: `nscm=1000` (remove COM motion). No `ntr` or `restraintmask` -- system is now fully relaxed.

### GB Implicit Solvent Relaxation Protocol (Compact 3-Step)

Implicit solvent simulations use `ntb=0`, `igb=8`, `gbsa=3`, `cut=1000.0` (effectively infinite), and no pressure control.

**Step 1: Minimization (no restraints, no SHAKE)**
```
GB minimization
 &cntrl
  imin=1, maxcyc=100,
  ntx=1, ntwr=100, ntpr=10,
  ioutfm=0, ntxo=1,
  cut=1000.0, ntb=0,
  igb=8, gbsa=3, surften=0.007, saltcon=0.0,
 /
```
Key: `cut=1000.0` (effectively infinite). No SHAKE (ntc=1 default). No PBC (ntb=0).

**Step 2: Heating (NVT, 100K to 300K, backbone restraints)**
```
GB heating with backbone restraints
 &cntrl
  imin=0, nstlim=500000, dt=0.002,
  ntx=1, irest=0, ig=-1,
  ntt=3, gamma_ln=1.0, temp0=100.0,
  ntc=2, ntf=2, nscm=1000,
  ntwx=5000, ntwr=500, ntpr=5000,
  cut=1000.0, igb=8, gbsa=3, surften=0.007,
  ntb=0, saltcon=0.0,
  ntr=1, restraintmask='@CA,N,C,O', restraint_wt=10.0,
  nmropt=1,
 /
&wt TYPE='TEMP0', istep1=0, istep2=500000, value1=100.0, value2=300.0 /
&wt TYPE='END' /
```
Key: `dt=0.002` (2 fs with HMR). Backbone restraints at 10 kcal/mol. `nscm=1000` for COM removal.

**Step 3: Equilibration (two stages with decreasing restraints)**
```
GB eq1 - backbone restraint 1.0 kcal/mol
 &cntrl
  imin=0, nstlim=125000, dt=0.002,
  ntx=5, irest=1, ig=-1,
  ntt=3, gamma_ln=1.0, temp0=300.0,
  ntc=2, ntf=2, nscm=1000,
  ntwx=5000, ntwr=500, ntpr=5000,
  cut=1000.0, igb=8, gbsa=3, surften=0.007,
  ntb=0, saltcon=0.0,
  ntr=1, restraintmask='@CA,N,C,O', restraint_wt=1.0,
 /
```
```
GB eq2 - backbone restraint 0.1 kcal/mol
 &cntrl
  ... (same as eq1 but restraint_wt=0.1)
 /
```

## Reference Tables

### Standard Explicit Solvent Relaxation Protocol
| Step | Type | Ensemble | Duration | dt | Restraint Mask | restraint_wt | Purpose |
|------|------|----------|----------|----|---------------|-------------|---------|
| 1 | Minimize | NVT | 1000 cyc | -- | :1-81 (solute) | 100 | Remove bad contacts |
| 2 | Heat | NVT | 1 ns | 1 fs | :1-81 (solute) | 100 | Ramp 100K to 298K |
| 3 | MD | NPT | 1 ns | 1 fs | :1-81 (solute) | 100 | Relax density |
| 4 | MD | NPT | 1 ns | 1 fs | :1-81 (solute) | 10 | Reduce solute restraint |
| 5 | Minimize | NVT | 1000 cyc | -- | @CA,N,C (backbone) | 10 | Re-minimize |
| 6 | MD | NPT | 1 ns | 1 fs | @CA,N,C (backbone) | 10 | Relax with backbone |
| 7 | MD | NPT | 1 ns | 1 fs | @CA,N,C (backbone) | 1 | Reduce backbone restraint |
| 8 | MD | NPT | 1 ns | 1 fs | @CA,N,C (backbone) | 0.1 | Further reduce restraint |
| 9 | MD | NPT | 1 ns | 1 fs | none | 0 | Final unrestrained relaxation |

### Restraint Reduction Schedule
| Stage | Restraint Target | Weight (kcal/mol/Ang^2) |
|-------|-----------------|------------------------|
| Heavy | All solute atoms | 100.0 |
| Medium | All solute atoms | 10.0 |
| Backbone heavy | CA, N, C (backbone) | 10.0 |
| Backbone medium | CA, N, C (backbone) | 1.0 |
| Backbone light | CA, N, C (backbone) | 0.1 |
| None | -- | 0.0 |

### Explicit vs. GB Implicit Solvent Comparison
| Feature | Explicit Solvent | GB Implicit Solvent |
|---------|-----------------|---------------------|
| PBC | ntb=1 (NVT) or ntb=2 (NPT) | ntb=0 (no PBC) |
| Pressure control | ntp=1, barostat=2 | Not applicable |
| Nonbonded cutoff | cut=8.0 (PME) | cut=1000.0 (effectively infinite) |
| Thermostat | ntt=3, gamma_ln=1.0-5.0 | ntt=3, gamma_ln=1.0 |
| COM removal | nscm=0 (heating) to nscm=1000 | nscm=1000 always |
| Heating timestep | dt=0.001 (1 fs) | dt=0.002 (2 fs with HMR) |
| Minimization SHAKE | ntc=2, ntf=2 (explicit) | ntc=1, ntf=1 (no SHAKE) |
| Solvent atoms | ~10K-100K water molecules | None (GB replaces solvent) |
| Speed | Slow (many particles) | Fast (few particles) |

### Minimizer Methods (ntmin)
| ntmin | Method | Notes |
|-------|--------|-------|
| 0 | Full conjugate gradient | 4 steepest descent cycles at start and after pairlist update |
| 1 | SD then CG (default) | `ncyc` cycles of steepest descent, then conjugate gradient |
| 2 | Steepest descent only | Slower convergence but robust for bad initial geometries |
| 3 | XMIN | Limited-memory BFGS (Section 26.8.4) |
| 4 | LMOD | Low-mode conformational search (Section 26.8.5) |
| 5 | DL-Find | Interface to external DL-Find library (Section 26.9) |

## Worked Example

### Explicit Solvent: Lysozyme Relaxation (Two-Stage Minimization)

**Input files** for a protein-ligand system in explicit water:

```bash
# Stage 1: Restrained minimization (solute frozen)
pmemd -O -i min1.in -p system.prmtop -c system.inpcrd \
    -o min1.out -r min1.rst -ref system.inpcrd
```

**min1.in:**
```
Restrained minimization - solute fixed
&cntrl
    imin=1, maxcyc=1000, ncyc=500,
    ntb=1, ntr=1, cut=10.0,
    restraintmask=':1-161', restraint_wt=500.0,
/
Hold solute fixed
500.0
RES 1 161
END
END
```

```bash
# Stage 2: Unrestrained minimization (whole system)
pmemd -O -i min2.in -p system.prmtop -c min1.rst \
    -o min2.out -r min2.rst
```

**min2.in:**
```
Unrestrained minimization - whole system
&cntrl
    imin=1, maxcyc=2500, ncyc=1000,
    ntb=1, ntr=0, cut=10.0,
/
```

```bash
# Stage 3: Heating (200 ps, 0K to 300K, NVT, weak restraints)
pmemd.cuda -O -i heat.in -p system.prmtop -c min2.rst \
    -ref min2.rst -o heat.out -r heat.rst -x heat.nc
```

**heat.in:**
```
200 ps heating with weak restraint on solute
&cntrl
    imin=0, irest=0, ntx=1,
    ntb=1, cut=10.0,
    ntr=1, restraintmask=':1-161', restraint_wt=10.0,
    ntc=2, ntf=2,
    tempi=0.0, temp0=300.0,
    ntt=3, gamma_ln=1.0,
    nstlim=100000, dt=0.002,
    ntpr=1000, ntwx=1000, ntwr=10000,
/
```

```bash
# Stage 4: NPT Equilibration (100 ns, 300K, 1 atm)
pmemd.cuda -O -i equil.in -p system.prmtop -c heat.rst \
    -o equil.out -r equil.rst -x equil.nc
```

**equil.in:**
```
100 ns NPT equilibration
&cntrl
    imin=0, irest=1, ntx=7,
    ntb=2, pres0=1.0, ntp=1, taup=2.0,
    cut=10.0, ntr=0,
    ntc=2, ntf=2,
    tempi=300.0, temp0=300.0,
    ntt=3, gamma_ln=1.0,
    nstlim=50000000, dt=0.002,
    ntpr=100000, ntwx=100000, ntwr=1000000,
/
```

**Verification checklist after each stage:**
1. Energy decreases monotonically during minimization
2. RMS gradient (GMAX, RMS) decreases to acceptable levels
3. Temperature reaches and maintains target during heating
4. Density stabilizes to ~1.0 g/mL during NPT (explicit solvent)
5. No SHAKE failures or vlimit warnings
6. Visual inspection: structure intact, no unnatural distortions

## Key Takeaways

1. **Always minimize before dynamics.** Start with heavy restraints on solute (100-500 kcal/mol), then reduce or remove. Steepest descent (ncyc) handles bad contacts; conjugate gradient (maxcyc-ncyc) refines the minimum.

2. **Heat gradually, not instantly.** Start at 100K (never 0K), ramp to target temperature over hundreds of ps using `nmropt=1` and `&wt` blocks. Keep solute restrained during heating. Use NVT (ntb=1) before switching to NPT.

3. **Use NPT to relax density.** After NVT heating, switch to NPT (ntb=2, ntp=1) with Monte Carlo barostat (barostat=2) to let the box adjust to proper density. Keep restraints on during density equilibration.

4. **Reduce restraints progressively.** Follow the staged reduction: 100 -> 10 -> 1 -> 0.1 kcal/mol, then shift from all-atom solute to backbone-only, then remove entirely. This prevents sudden structural drift.

5. **GB systems differ fundamentally.** No periodic boundaries (ntb=0), no pressure control, cutoff=1000 (effectively infinite), no SHAKE during minimization. Use fewer but longer relaxation stages. Always set `nscm` to remove COM motion.

## Connects To
- **Chapter 7: MD Engines and mdin Reference** -- complete namelist variable reference for &cntrl, &ewald, &wt
- **Chapter 9: Production MD** -- transitioning from relaxation to multi-ns production runs
- **Chapter 20: Implicit Solvent (GB, PBSA, RISM)** -- details on igb, gbsa, saltcon, GBneck2 models
- **Chapter 3: LEaP System Building** -- creates the prmtop/inpcrd files that are the starting point for relaxation
- **Chapter 6: ParmEd** -- HMR (hydrogen mass repartitioning) for 4 fs timesteps, used in GB protocols
- Sections 23.6 of the Amber 2026 Reference Manual (PDF pages 425-455) -- complete variable documentation
- Tutorials: `03-1_Relaxation of Explicit Water Systems`, `03-2_Relaxation of Implicit Solvent System GB`, `07_Creating_Stable_Systems_and_Running_MD`