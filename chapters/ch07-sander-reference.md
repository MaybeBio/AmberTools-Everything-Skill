# Chapter 7: sander Full Reference

## Core Commands & Syntax

```bash
sander [-help] [-O] [-A] -i mdin -o mdout -p prmtop -c inpcrd -r restrt \
       -ref refc -mtmd mtmd -x mdcrd -y inptraj -v mdvel -frc mdfrc -e mden \
       -inf mdinfo -radii radii -cpin cpin -cpout cpout -cprestrt cprestrt \
       -cein cein -ceout ceout -cerestrt cerestrt -evbin evbin -suffix suffix
```
- `-O`: Overwrite output files if they exist.
- `-A`: Append output files if they exist (used mainly for replica exchange).

### Key File Arguments
| File | Type | Description |
|------|------|-------------|
| `mdin` | input | Control data for the min/md run |
| `mdout` | output | User-readable state info and diagnostics |
| `prmtop` | input | Molecular topology, force field, periodic box type |
| `inpcrd` | input | Initial coordinates and optionally velocities/box |
| `restrt` | output | Final coordinates, velocities, box for restart |
| `refc` | input | Reference coords for position restraints / targeted MD |
| `mdcrd` | output | Coordinate trajectory |
| `mdvel` | output | Velocity trajectory |
| `mdfrc` | output | Force trajectory |
| `mden` | output | Energy data over trajectory |
| `inptraj` | input | Trajectory for analysis (imin=5,6) |
| `cpin` | input | Protonation state definitions (constant pH) |
| `cpout` | output | Protonation state data over trajectory |
| `evbin` | input | EVB potentials |

## Namelist Input Syntax

```
&name
 var1=value, var2=value, var3(sub)=value,
 var4=1,2,3,4,
 var5=25*3.1415,
/
```

- Variables are case-insensitive; string constants use single quotes.
- Repeat counts: `25*3.1415` = 25 copies of 3.1415.
- Array values go into successive locations.
- The ending `/` is standard; `&end` is non-standard but accepted.

## Complete Namelist Variable Reference

### &cntrl -- General Flags

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `imin` | int | 0 | 0=MD, 1=minimization, 5=trajectory analysis (min), 6=trajectory analysis (MD), 7=i-PI socket |
| `nmropt` | int | 0 | 0=no NMR, 1=NMR restraints+weight changes, 2=+NOESY/chem shifts/dipolar |
| `ntx` | int | 1 | 1=coords only, 5=coords+velocities |
| `irest` | int | 0 | 0=new simulation, 1=restart (requires ntx=5) |
| `ntxo` | int | 2 | 1=ASCII formatted restrt, 2=NetCDF restrt |
| `ntpr` | int | 50 | Print energy info every ntpr steps to mdout/mdinfo |
| `ntave` | int | 0 | Steps for running averages; 0=disabled |
| `ntwr` | int | nstlim | Write restart every ntwr steps; <0 = unique restrt_N |
| `iwrap` | int | 0 | 1=wrap coords into primary box |
| `ntwx` | int | 0 | Write coords to mdcrd every ntwx steps; 0=disabled |
| `ntwv` | int | 0 | Write velocities to mdvel; -1=combined mdcrd (NetCDF only) |
| `ntwf` | int | 0 | Write forces to mdfrc; -1=combined mdcrd (NetCDF only) |
| `ntwe` | int | 0 | Write energies to mden; 0=disabled |
| `ioutfm` | int | 1 | 0=ASCII trajectory, 1=NetCDF binary |
| `ntwprt` | int | 0 | 0=all atoms; >0=first N atoms in trajectory |
| `idecomp` | int | 0 | 0=off; 1-4=per-residue energy decomposition |
| `ionstepvelocities` | int | 0 | 0=half-step-ahead velocities; 1=on-step velocities |

### &cntrl -- Frozen/Restrained Atoms

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ibelly` | int | 0 | 1=belly dynamics (freeze non-selected atoms) |
| `bellymask` | str | '' | Mask for moving atoms when ibelly=1 |
| `ntr` | int | 0 | 1=Cartesian restraints |
| `restraint_wt` | float | 0 | Force constant (kcal/mol/A^2) |
| `restraintmask` | str | '' | Mask for restrained atoms |

### &cntrl -- Minimization

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `maxcyc` | int | 1 | Max minimization cycles |
| `ncyc` | int | 10 | Steepest descent cycles before conjugate gradient |
| `ntmin` | int | 1 | 0=CG, 1=SD then CG, 2=SD only, 3=XMIN, 4=LMOD, 5=DL-Find |
| `dx0` | float | 0.01 | Initial step length |
| `drms` | float | 1e-4 | Convergence criterion (kcal/mol/A) |

### &cntrl -- Molecular Dynamics

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `nstlim` | int | 1 | Number of MD steps |
| `nscm` | int | 1000 | COM motion removal frequency |
| `t` | float | 0.0 | Start time (ps) |
| `dt` | float | 0.001 | Time step (ps); max 0.002 w/ SHAKE |
| `nrespa` | int | 1 | Multiple time-stepping; >1 evaluates slow forces less frequently |

### &cntrl -- Temperature Regulation

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ntt` | int | 0 | 0=NVE, 1=Berendsen, 2=Andersen, 3=Langevin, 9=OIN (Nose-Hoover chain), 10=SINR, 11=Bussi |
| `temp0` | float | 300.0 | Target temperature (K) |
| `tempi` | float | 0.0 | Initial temp; 0.0=from forces; Maxwellian at tempi if ntx=1 |
| `temp0les` | float | -1 | Target temp for LES particles; <0=single bath |
| `ig` | int | -1 | Random seed; -1=use date/time |
| `tautp` | float | 1.0 | Heat bath coupling time (ps) for ntt=1 |
| `gamma_ln` | float | 0.0 | Collision frequency (ps^-1) for ntt=3; must be >0 for ntt=9,10 |
| `vrand` | int | 1000 | Velocity randomization frequency for ntt=2 |
| `vlimit` | float | 20.0 | Max velocity component; GPU default=-1 (disabled) |
| `nkija` | int | 1 | Thermostat substeps for ntt=9,10 |
| `idistr` | int | 0 | Thermostat velocity distribution accumulation freq (ntt=9) |
| `sinrtau` | float | 1.0 | Time scale for SINR thermostat masses (ntt=10) |

### &cntrl -- Pressure Regulation

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ntp` | int | 0 | 0=no scaling, 1=isotropic, 2=anisotropic, 3=semiisotropic, 4=targeted volume |
| `barostat` | int | 1 | 1=Berendsen, 2=Monte Carlo |
| `mcbarint` | int | 100 | Steps between MC volume attempts |
| `pres0` | float | 1.0 | Reference pressure (bar) |
| `comp` | float | 44.6 | Compressibility (10^-6 bar^-1) |
| `taup` | float | 1.0 | Pressure relaxation time (ps) |
| `baroscalingdir` | int | 0 | 0=random, 1=x, 2=y, 3=z (MC barostat+anisotropic) |
| `csurften` | int | 0 | 0=off; 1=yz, 2=xz, 3=xy interfaces |
| `gamma_ten` | float | 0.0 | Surface tension (dyne/cm) |
| `ninterface` | int | 2 | Number of interfaces |

### &cntrl -- "Middle" Scheme (Extended Thermostat/Barostat)

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ischeme` | int | 0 | 0=conventional, 1=LFMiddle leapfrog |
| `ithermostat` | int | 1 | 1=Langevin, 2=Andersen (middle scheme) |
| `baro_stochastic` | int | 0 | 1=stochastic cell rescaling (with barostat=1) |
| `therm_par` | float | 5.0 | Thermostat parameter (ps^-1), must be >0 |

### &cntrl -- SHAKE Constraints

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ntc` | int | 1 | 1=no SHAKE, 2=H-bonds only, 3=all bonds |
| `ntf` | int | 1 | 1=all forces, 2=omit H-bond forces, 3-8=omit more |
| `tol` | float | 1e-5 | SHAKE tolerance (A) |
| `jfastw` | int | 0 | 0=fast water SHAKE; 4=disable fast water |
| `noshakemask` | str | '' | Atoms NOT to SHAKE (forces ntf=1) |
| `watnam` | str | 'WAT ' | Water residue name |
| `owtnm` | str | 'O ' | Water oxygen atom name |
| `hwtnam1` | str | 'H1 ' | Water H1 atom name |
| `hwtnam2` | str | 'H2 ' | Water H2 atom name |

### &cntrl -- Potential Function

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ntb` | int | auto | 0=no PBC, 1=const vol, 2=const pressure |
| `dielc` | float | 1.0 | Dielectric constant for electrostatics |
| `cut` | float | 8.0 | Nonbonded cutoff (A); GB default=9999 |
| `fswitch` | float | -1 | Force switching distance; <=0=off |
| `nsnb` | int | 25 | Nonbonded list update freq (nbflag=0) |
| `ipol` | int | 0 | 1=polarizable force field |
| `ipgm` | int | 0 | 1=pGM polarizable Gaussian Multipole |
| `ifqnt` | int | 0 | 1=QM/MM (requires &qmmm) |
| `igb` | int | 0 | GB implicit solvent model |
| `ipb` | int | 0 | Poisson-Boltzmann implicit solvent |
| `irism` | int | 0 | 3D-RISM solvation |
| `ievb` | int | 0 | 1=Empirical Valence Bond |
| `iamoeba` | int | 0 | 1=Amoeba polarizable potential |
| `lj1264` | int | 0 | 1=12-6-4 Lennard-Jones |
| `plj1264` | int | 0 | 1=pairwise 12-6-4 LJ |
| `ips` | int | 0 | 0=off; 1-6=Isotropic Periodic Sum |
| `infe` | int | 0 | 1=enable non-equilibrium free energy (&smd, &abmd, &bbmd) |

### &cntrl -- Water Cap, EMAP, NMR

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ivcap` | int | 0 | 0=from prmtop, 1=cap from box, 2=inactivate, 5=shell |
| `fcap` | float | 0 | Cap restraint force constant |
| `cutcap` | float | 0 | Cap radius |
| `xcap, ycap, zcap` | float | 0 | Cap center |
| `iemap` | int | 0 | 1=EMAP restraints on |
| `gammamap` | float | 1.0 | EMAP friction constant (1/ps) |
| `iscale` | int | 0 | Additional scaling parameters (NMR) |
| `noeskp` | int | 1 | NOESY eval frequency |
| `ipnlty` | int | 1 | Penalty: 1=abs, 2=square, 3=1/6 power |
| `mxsub` | int | 1 | Max submolecules for NOESY |
| `scalm` | float | 100.0 | Mass for scaling parameters |
| `pencut` | float | 0.1 | Penalty cutoff for summaries |
| `tausw` | float | 0.1 | Mixing time threshold (s) |

### &cntrl -- SGLD (Self-Guided Langevin Dynamics)

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `isgld` | int | 0 | 0=off; 1=uniform; 2=atom-specific; 3=balanced |
| `tsgavg` | float | 0.2 | Local averaging time (ps) |
| `sgft` | float | 0.0 | Momentum guiding factor (-1 to 1) |
| `sgff` | float | 0.0 | Force guiding factor (-0.32 to 0.32) |
| `sgfg` | float | 0.0 | SGLD-GLE momentum guiding factor (-1 to 1) |
| `tempsg` | float | 0.0 | Effective guiding temperature |
| `isgsta` | int | 1 | First atom of SGLD region |
| `isgend` | int | natom | Last atom of SGLD region |
| `sgmask` | str | ':*' | SGLD atom mask |
| `sgtype` | int | 1 | Spatial average: 1=none, 2=bonded, 3=bond+angle+dihedral, 4=cutoff |
| `sgsize` | float | 3.0 | Cutoff for local spatial average (A) |

### &cntrl -- Targeted MD / MTMD

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `itgtmd` | int | 0 | 0=off, 1=targeted MD, 2=MTMD |
| `tgtrmsd` | float | 0.0 | Target RMSD |
| `tgtmdfrc` | float | 0.0 | Force constant |
| `tgtfitmask` | str | '' | Atoms for RMS fit |
| `tgtrmsmask` | str | '' | Atoms for RMSD/restraint force |

### &cntrl -- Electric Field (pmemd only)

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `efx, efy, efz` | float | 0.0 | Electric field components (kcal/(mol*A*e)) |
| `efn` | int | 0 | 1=normalize to box size |
| `efphase` | float | 0 | Phase offset (degrees) |
| `effreq` | float | 0 | Frequency (timestep units) |

### &cntrl -- MC Water / RAMD / Reweight (pmemd only)

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `mcwat` | int | 0 | 1=MC water equilibration |
| `nmd` | int | 1000 | MD steps per MC/MD cycle |
| `nmc` | int | 100000 | MC steps per cycle |
| `mcwatmask` | str | '' | Center of MC region |
| `mcligshift` | float | auto | MC region half-length |
| `ramdboost` | float | 1 | RAMD acceleration boost |
| `ramdboostfreq` | int | 0 | RAMD boost increase frequency |
| `ramdboostrate` | float | 0 | RAMD boost increase rate |
| `ramdint` | int | 0 | RAMD boost interval |
| `ramdmaxdist` | float | 0 | RAMD termination distance |
| `ramdligmask` | str | '' | RAMD ligand mask |
| `ramdprotmask` | str | '' | RAMD protein mask |
| `reweight` | int | 0 | 1=trajectory re-evaluation |
| `midpoint` | int | 0 | 1=3D spatial decomposition (pmemd.MPI only) |

### &cntrl -- TI Decomposition

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ntwd` | int | 0 | 1=TI free energy decomposition |
| `decompmask` | str | '' | Atoms for dV/dL output |
| `ligmask` | str | '' | Ligand atoms for region decomposition |
| `proteinmask` | str | '' | Protein atoms for region decomposition |
| `cofactormask` | str | '' | Cofactor atoms for region decomposition |

### &cntrl -- pmemd-specific

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `mdout_flush_interval` | int | 300 | Min seconds between mdout flushes (0-3600) |
| `mdinfo_flush_interval` | int | 60 | Min seconds between mdinfo flushes (0-3600) |
| `es_cutoff` | float | auto | Electrostatic cutoff (pmemd CPU only) |
| `vdw_cutoff` | float | auto | VDW cutoff; must >= es_cutoff |
| `no_intermolecular_bonds` | int | 1 | 1=fuse bonded molecules; 0=prmtop definition |
| `ene_avg_sampling` | int | ntpr | Steps between energy samples for averages |

### &ewald -- Particle Mesh Ewald

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `nfft1, nfft2, nfft3` | int | auto | Charge grid dimensions; ~1A spacing recommended |
| `order` | int | 4 | B-spline interpolation order (min 3) |
| `verbose` | int | 0 | 0-3, higher=more output |
| `ew_type` | int | 0 | 0=PME, 1=regular Ewald |
| `dsum_tol` | float | 1e-5 | Direct sum tolerance; RMS force error ~10-50x this |
| `rsum_tol` | float | 5e-5 | Reciprocal sum tolerance |
| `mlimit(1,2,3)` | int | auto | Reciprocal vectors for regular Ewald |
| `ew_coeff` | float | auto | Ewald coefficient (A^-1); determined from dsum_tol and cut |
| `nbflag` | int | 1 | 0=nsnb-based list, 1=skin-based list |
| `skinnb` | float | 2.0 | Nonbonded skin width (A) |
| `skin_permit` | float | 0.5 | pmemd.cuda: pair list rebuild threshold (0.5-1.0) |
| `nbtell` | int | 0 | 1=print list update messages |
| `netfrc` | int | 1 | 1=remove net force; set 0 when ntr>0 |
| `vdwmeth` | int | 1 | 0=no correction, 1=continuum VDW correction |
| `eedmeth` | int | 1 | Coulomb switch: 1=cubic spline, 2=linear table, 3=exact |
| `eedtbdns` | float | 500 | Spline/table density (points/unit) |
| `column_fft` | int | 0 | 1=column-mode FFT (parallel) |
| `use_axis_opt` | int | auto | pmemd: 1=force axis optimization; 0=disable |
| `fft_grids_per_ang` | float | 1.0 | pmemd: reciprocal grid density |
| `block_fft` | int | auto | pmemd: block/pencil FFT |
| `fft_blk_y_divisor` | int | auto | pmemd: block FFT division |
| `excl_recip` | int | auto | pmemd: exclude reciprocal from direct |
| `excl_master` | int | auto | pmemd: reserve master for I/O |
| `atm_redist_freq` | int | auto | pmemd: atom redistribution frequency |

### &ewald -- Extra Points / Polarizable

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `frameon` | int | 1 | 1=remove EP bonds/angles, transfer torques |
| `chngmask` | int | 1 | 1=rebuild 1-2,1-3,1-4 masks for EP |
| `indmeth` | int | 3 | 0-2=iterative; 3=Car-Parrinello (default) |
| `diptol` | float | 1e-4 | Dipole convergence (Debye) |
| `maxiter` | int | 20 | Max iterations for indmeth<3 |
| `dipmass` | float | 0.33 | Fictitious dipole mass |
| `diptau` | float | 11.0 | Dipole temp coupling; >10=off |
| `irstdip` | int | 0 | 1=read dipole restart |
| `scaldip` | int | 1 | 1=scale 1-4 charge-dipole/dipole-dipole |

### &ewald -- IPS (Isotropic Periodic Sum)

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `raips` | float | -1 | Local region radius; -1=longest box side |
| `mipsx, mipsy, mipsz` | int | -1 | DFFT grid counts |
| `mipso` | int | 4 | IPS B-spline order |
| `gridips` | float | 2.0 | DFFT grid size (A) |
| `dvbips` | float | 1e-8 | Volume tolerance for IPS grid update |

### &ewald -- MPI Timing

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `profile_mpi` | int | 0 | 1=detailed per-thread MPI timings |

### &wt -- Varying Conditions (read if nmropt>0)

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `type` | str | (required) | Quantity to vary (see list below) |
| `istep1, istep2` | int | 0 | Step range; istep2=0 means until end |
| `value1, value2` | float | 0 | Values at istep1 and istep2 |
| `iinc` | int | 0 | 0=continuous; >0=step function interval |
| `imult` | int | 0 | 0=linear; 1=multiplicative scaling |

**Valid TYPE values**: `TEMP0`, `TEMP0LES`, `TAUTP`, `CUT`, `BOND`, `ANGLE`, `TORSION`, `IMPROP`, `VDW`, `HB`, `ELEC`, `NB`, `ATTRACT`, `REPULSE`, `RSTAR`, `INTERN`, `ALL`, `REST`, `RESTS`, `RESTL`, `NOESY`, `SHIFTS`, `SHORT`, `TGTRMSD`, `NSTEP0`, `STPMLT`, `DISAVE`, `ANGAVE`, `TORAVE`, `DISAVI`, `ANGAVI`, `TORAVI`, `DUMPFREQ`, `END`.

### &emap -- CryoEM Map Restraints

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `mapfile` | str | '' | Map (.map, .ccp4, .mrc) or PDB file; ''=from coordinates |
| `atmask` | str | ':*' | Atom mask for restrained atoms |
| `fcons` | float | 0.05 | Restraint force constant (kcal/g) |
| `move` | int | 0 | 0=fixed map; >0=allow map to move |
| `resolution` | float | 2.0 | Resolution for structure->map conversion (A); <0=boundary |
| `ifit` | int | 0 | 0=no fit; 1=fit map to coords; 2=fit coords to map |
| `grids` | int(6) | 1,1,1,1,1,1 | Grid points for fitting (x,y,z,phi,psi,theta) |
| `mapfit` | str | '' | Output filename for final map |
| `molfit` | str | '' | Output filename for final coordinates (.pdb) |

### &smd -- Steered Molecular Dynamics

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `output_file` | str | 'nfe-smd.txt' | SMD output filename |
| `output_freq` | int | 50 | Output frequency (MD steps) |
| `cv_file` | str | 'nfe-smd-cv' | Collective variable definition file |

**&colvar entries for SMD** (in cv_file):
| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `cv_type` | str | (required) | Reaction coordinate type |
| `cv_ni` | int | 0 | Number of cv_i integers |
| `cv_nr` | int | 0 | Number of cv_r reals |
| `cv_i` | int[] | - | Atom indices |
| `cv_r` | float[] | - | Coefficients/reference positions |
| `path` | float[] | (required) | Steering path points |
| `npath` | int | 0 | Number of path elements |
| `path_mode` | str | 'SPLINE' | 'SPLINE' or 'LINES' |
| `harm` | float[] | (required) | Harmonic constant(s) |
| `nharm` | int | 0 | Number of harm elements |
| `harm_mode` | str | 'SPLINE' | 'SPLINE' or 'LINES' |

**Valid cv_type values**: `DISTANCE`, `COM_DISTANCE`, `DF_COM_DISTANCE`, `LCOD`, `ANGLE`, `COM_ANGLE`, `TORSION`, `COM_TORSION`, `COS_OF_DIHEDRAL`, `SIN_OF_DIHEDRAL`, `PAIR_DIHEDRAL`, `PATTERN_DIHEDRAL`, `R_OF_GYRATION`, `MULTI_RMSD`, `N_OF_BONDS`, `N_OF_STRUCTURES`, `HANDEDNESS`.

### &abmd -- Adaptively Biased MD

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `mode` | str | (required) | 'ANALYSIS', 'UMBRELLA', or 'FLOODING' |
| `cv_file` | str | (required) | Collective variable definition file |
| `monitor_file` | str | (required) | Monitor output file |
| `monitor_freq` | int | (required) | Monitor output frequency |
| `timescale` | float | (required) | Flooding timescale tau_F (ps) |
| `umbrella_file` | str | (required) | Biasing potential file (for UMBRELLA mode) |
| `snapshots_basename` | str | '' | Snapshot basename for biasing potential |
| `snapshots_freq` | int | 0 | Snapshot frequency (0=disabled) |
| `selection_freq` | int | 0 | Multiple-walker resampling frequency |
| `selection_constant` | float | 0 | Resampling parameter C |
| `selection_epsilon` | float | 1 | Stopping criterion epsilon (0-1) |
| `wt_temperature` | float | inf | Pseudo-temperature T' for well-tempered |
| `wt_umbrella_file` | str | '' | True biasing potential output |

**&colvar entries for ABMD** (in cv_file):
| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `cv_min` | float | auto | Minimum CV value |
| `cv_max` | float | auto | Maximum CV value |
| `resolution` | float | (required) | Spatial resolution |

### &bbmd -- BBMD (same as &abmd, different internal flavor)

Uses identical variables as &abmd. The difference is purely technical. Both are activated by `infe=1` in &cntrl.

### &qmmm -- QM/MM

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `qmmask` | str | '' | Mask for QM region atoms |
| `qmcharge` | int | 0 | Total charge of QM region |
| `spin` | int | 1 | Spin multiplicity |
| `qm_theory` | str | '' | QM method: 'PM3', 'AM1', 'MNDO', 'PDDG/PM3', 'PM3CARB1', 'DFTB', 'DFTB2', 'DFTB3', 'EXTERN' |
| `qmcut` | float | 0.0 | QM/MM cutoff (A) |
| `qmshake` | int | 1 | 0=no SHAKE in QM region |
| `qmqmdx` | str | '' | External QM program path |
| `verbosity` | int | 0 | QM/MM output verbosity |
| `qmmm_switch` | int | 0 | 1=QM/MM switching function on |
| `qmmm_int` | int | 0 | QM/MM interaction type |
| `qmmm_ewald` | int | 0 | 1=QM/MM Ewald for periodic |
| `qmmm_pme` | int | 0 | 1=PME for QM/MM charges |
| `qm_ewald` | int | 0 | 1=QM Ewald sum |
| `qm_pme` | int | 0 | 1=PME for QM region |
| `dftb_maxiter` | int | 70 | DFTB max SCF iterations |
| `dftb_telec` | float | 300.0 | DFTB electronic temperature |
| `dftb_chg` | int | 0 | DFTB initial charge |
| `tight_pmax` | int | auto | Maximum tight-binding P |
| `scfconv` | float | 1e-8 | SCF convergence |
| `dftb_3rd_order` | str | '' | DFTB3 3rd order parameters |
| `dftb_dispersion` | int | 0 | DFTB dispersion correction |
| `dftb_hubbard_derivs` | str | '' | DFTB Hubbard derivatives |
| `vsolv` | int | 0 | Adaptive solvent QM/MM: 0=off, 2=ONIOM-XS, 3=Hot Spot |
| `qxd` | int | 0 | 1=charge exchange decomposition |

### &rism -- 3D-RISM

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `irism` | int | 0 | 1=enable 3D-RISM (in &cntrl) |
| `closure` | str | 'KH' | 'KH', 'HNC', 'PSE-n' |
| `gfCorrection` | int | 0 | 1=GF excess chemical potential |
| `pcpluscorrection` | int | 0 | 1=PC+/3D-RISM correction |
| `uccoeff` | float(4) | 0,0,0,0 | UC correction coefficients |
| `periodic` | str | '' | 'pme' or 'ewald' for periodic |
| `solvcut` | float | auto | LJ cutoff for periodic |
| `grdspc` | float(3) | 0.5,0.5,0.5 | Grid spacing (A) |
| `asympcorr` | bool | .true. | Long-range asymptotics |
| `treeDCF` | bool | .true. | Treecode for DCF asymptotics |
| `treeTCF` | bool | .true. | Treecode for TCF asymptotics |
| `treeCoulomb` | bool | .false. | Treecode for Coulomb |
| `treeDCFMAC` | float | 0.1 | DCF multipole acceptance |
| `treeTCFMAC` | float | 0.1 | TCF multipole acceptance |
| `treeCoulombMAC` | float | 0.1 | Coulomb multipole acceptance |
| `treeDCFOrder` | int | 2 | DCF Taylor order |
| `treeTCFOrder` | int | 2 | TCF Taylor order |
| `treeCoulombOrder` | int | 2 | Coulomb Taylor order |
| `treeDCFN0` | int | 500 | DCF leaf cluster size |
| `treeTCFN0` | int | 500 | TCF leaf cluster size |
| `rismnrespa` | int | 1 | RISM multiple time-stepping |

### &debugf -- Debug Options

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `do_debugf` | int | 0 | 1=enable debug |
| `atomn` | int(25) | 0 | Atom numbers for force test |
| `nranatm` | int | 0 | Random atoms to test |
| `ranseed` | int | 71277 | Random seed for atom selection |
| `neglgdel` | int | 5 | -log10(delta) for numerical diff (default delta=1e-5 A) |
| `chkvir` | int | 0 | 1=test virial numerically |
| `dumpfrc` | int | 0 | 1=dump forces to forcedump.dat |
| `rmsfrc` | int | 0 | 1=compare forces to forcedump.dat |
| `zerochg` | int | 0 | 1=zero all charges |
| `zerovdw` | int | 0 | 1=remove all VDW |
| `zerodip` | int | 0 | 1=remove atomic dipoles |
| `do_dir, do_rec, do_adj, do_self, do_bond, do_cbond, do_angle, do_ephi, do_xconst, do_cap` | int | 1 | Toggle force subroutines on/off |

### &pmd -- Umbrella Sampling

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `output_file` | str | 'nfe-pmd.txt' | Output filename |
| `output_freq` | int | 50 | Output frequency |
| `cv_file` | str | (required) | CV definition file |

**&colvar entries for Umbrella Sampling**:
| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `anchor_position` | float(4) | 0,0,0,0 | Umbrella rectangle (r1,r2,r3,r4) |
| `anchor_strength` | float(2) | 0,0 | Force constants (k1,k2) |

### &pol_gauss -- pGM Model

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `pol_gauss_verbose` | int | 0 | 1=extra printing |
| `pol_gauss_ips` | int | 0 | 1=use IPS for pGM |
| `ee_dsum_cut` | float | 9.0 | Ewald direct sum cutoff (A) |
| `dipole_scf_tol` | float | 1e-2 | Induction convergence tolerance |
| `dipole_solv_opt` | int | 3 | 3=PCG solver; 4=SOR solver |
| `scf_cg_niter` | int | 50 | Max CG iterations |
| `scf_sor_coefficient` | float | 0.65 | SOR relaxation parameter |
| `scf_sor_niter` | int | 100 | Max SOR iterations |

## Common Workflows

### 1. Simple Minimization with Cartesian Restraints
```bash
# mdin:
&cntrl
  imin=1, maxcyc=200,
  ntpr=5,
  ntr=1, restraint_wt=1.0,
  restraintmask=':1-58',
/
sander -O -i mdin -p prmtop -c inpcrd -ref refc -r restrt -o mdout
```

### 2. Standard NPT MD Production Run
```bash
# mdin:
&cntrl
  imin=0, irest=1, ntx=5,
  ntt=3, temp0=300.0, gamma_ln=5.0,
  ntp=1, taup=2.0,
  ntb=2, ntc=2, ntf=2,
  nstlim=500000, dt=0.002,
  ntwx=1000, ntpr=200, ntwe=1000, ntwr=10000,
  ioutfm=1, iwrap=1,
  cut=8.0,
/
sander -O -i mdin -p prmtop -c inpcrd -r restrt -x mdcrd -o mdout
```

### 3. SGLD Enhanced Sampling Run
```bash
# mdin:
&cntrl
  imin=0, irest=0, ntx=1,
  ntt=3, temp0=300.0, gamma_ln=10.0,
  ntc=2, ntf=2,
  nstlim=500000, dt=0.002,
  ntwx=1000, ntpr=200,
  isgld=1, tsgavg=0.2, sgft=1.0,
/
sander -O -i mdin -p prmtop -c inpcrd -r restrt -x mdcrd -o mdout
```

### 4. QM/MM MD
```bash
# mdin:
&cntrl
  imin=0, irest=0, ntx=1,
  ntt=3, temp0=300.0, gamma_ln=2.0,
  ntc=2, ntf=2, nstlim=10000, dt=0.001, cut=8.0,
  ifqnt=1,
/
&qmmm
  qmmask=':ACE,ALA,NME',
  qmcharge=0, qmshake=0,
  qm_theory='PM3', qmcut=8.0,
/
sander -O -i mdin -p prmtop -c inpcrd -r restrt -x mdcrd -o mdout
```

### 5. EMAP-Constrained Simulation
```bash
# mdin:
&cntrl
  imin=0, nstlim=100000, ntc=2, ntf=2, cut=9.0,
  ntt=3, gamma_ln=10.0, dt=0.001, iemap=1,
/
&emap
  mapfile='map.ccp4', atmask=':1-56', fcons=0.1,
  move=1, ifit=1, mapfit='final.ccp4', molfit='final.pdb',
/
sander -O -i mdin -p prmtop -c inpcrd -r restrt -o mdout
```

### 6. ABMD Free Energy Calculation
```bash
# mdin:
&cntrl
  imin=0, irest=0, ntx=1, nstlim=1000000,
  ntt=3, temp0=300.0, gamma_ln=2.0, dt=0.002, cut=8.0,
  infe=1,
/
&abmd
  mode='FLOODING', cv_file='cv.in', monitor_file='monitor.txt',
  monitor_freq=100, timescale=100.0,
  umbrella_file='umbrella.nc', snapshots_basename='snap',
  snapshots_freq=10000,
/
sander -O -i mdin -p prmtop -c inpcrd -r restrt -o mdout
```

### 7. multisander (Replica Parallel)
```bash
mpirun -np 32 sander.MPI -ng 4 -groupfile groupfile
# groupfile:
# -O -p prmtop1 -c inpcrd1 -i replica1.mdin -suffix replica1
# -O -p prmtop2 -c inpcrd2 -i replica2.mdin -suffix replica2
# -O -p prmtop3 -c inpcrd3 -i replica3.mdin -suffix replica3
# -O -p prmtop4 -c inpcrd4 -i replica4.mdin -suffix replica4
```

## Performance Optimization

- **GPU**: Use pmemd.cuda instead of sander for supported features
- **Parallel**: Use sander.MPI with `mpirun -np N`; PME scales well to 32-64 cores
- **Column FFT**: Set `column_fft=1` in &ewald for high processor counts
- **NetCDF**: Always use `ioutfm=1` for smaller, faster I/O
- **Output frequency**: Set ntwx, ntpr, ntwr to >= 1000 for production; avoid frequent writes
- **Cutoff**: 8.0 A is standard for PME; consider 9.0 A for better accuracy
- **Grid spacing**: ~1.0 A (set nfft1/2/3 accordingly); must be powers of 2,3,5
- **HMR**: Hydrogen Mass Repartitioning allows dt=0.004 with SHAKE
- **RESPA**: nrespa>1 can speed up GB but may hurt PME scaling at high processor counts
- **Middle Scheme**: Use `ischeme=1, ithermostat=1, therm_par=5.0` for efficient NVT/NPT sampling

## Worked Example: Standard NPT Equilibration + Production

**Step 1: Minimization**
```bash
cat > min.in << 'EOF'
Initial minimization
&cntrl
  imin=1, maxcyc=1000, ncyc=500,
  ntb=1, ntc=2, ntf=2, cut=8.0,
  ntpr=100,
/
EOF
sander -O -i min.in -p system.prmtop -c system.inpcrd -r min.rst -o min.out
```

**Step 2: NVT Equilibration with restraints**
```bash
cat > heat.in << 'EOF'
NVT heating with restraints
&cntrl
  imin=0, irest=0, ntx=1,
  ntb=1, ntc=2, ntf=2, cut=8.0,
  ntt=3, temp0=300.0, tempi=0.0, gamma_ln=2.0,
  nstlim=50000, dt=0.002,
  ntr=1, restraint_wt=10.0, restraintmask=':1-100@CA',
  ntpr=500, ntwx=500, ntwr=5000, ioutfm=1,
/
EOF
sander -O -i heat.in -p system.prmtop -c min.rst -ref min.rst -r heat.rst -x heat.nc -o heat.out
```

**Step 3: NPT Equilibration**
```bash
cat > npt.in << 'EOF'
NPT density equilibration
&cntrl
  imin=0, irest=1, ntx=5,
  ntb=2, ntp=1, barostat=2, mcbarint=100,
  pres0=1.0, taup=2.0,
  ntc=2, ntf=2, cut=8.0,
  ntt=3, temp0=300.0, gamma_ln=2.0,
  nstlim=100000, dt=0.002,
  ntr=1, restraint_wt=5.0, restraintmask=':1-100@CA',
  ntpr=1000, ntwx=1000, ntwr=10000, ioutfm=1,
/
EOF
sander -O -i npt.in -p system.prmtop -c heat.rst -ref heat.rst -r npt.rst -x npt.nc -o npt.out
```

**Step 4: Production NPT**
```bash
cat > prod.in << 'EOF'
Production NPT
&cntrl
  imin=0, irest=1, ntx=5,
  ntb=2, ntp=1, barostat=2, mcbarint=100,
  pres0=1.0, taup=2.0,
  ntc=2, ntf=2, cut=8.0,
  ntt=3, temp0=300.0, gamma_ln=2.0,
  nstlim=50000000, dt=0.002,
  ntpr=5000, ntwx=5000, ntwr=50000, ntwe=5000, ioutfm=1, iwrap=1,
  ig=-1,
/
EOF
sander -O -i prod.in -p system.prmtop -c npt.rst -r prod.rst -x prod.nc -o prod.out
```

## Key Takeaways

1. **sander is the most feature-rich Amber engine** -- supports QM/MM, polarizable force fields, NMR restraints, LES, EMAP, ABMD, SMD, 3D-RISM, and more. Use it when pmemd lacks the needed feature.
2. **Always set netfrc=0 when using ntr>0** (Cartesian restraints) to avoid artifacts from force removal.
3. **Use ioutfm=1** (NetCDF) for all trajectories -- it is smaller, faster, higher precision, and required for many analysis tools.
4. **The "middle" scheme** (`ischeme=1, ithermostat=1`) provides more accurate configuration sampling and allows larger time steps.
5. **For production, prefer pmemd** over sander when your simulation is within pmemd's feature envelope -- it is significantly faster and scales better.

## Connects To

- **Chapter 1**: Force field parameters and topology preparation
- **Chapter 8**: pmemd CPU reference (faster alternative for standard MD)
- **Chapter 9**: pmemd GPU reference (pmemd.cuda, multi-GPU, performance)
- **Chapter 4**: Generalized Born implicit solvent
- **Chapter 6**: QM/MM and EVB methods
- **Chapter 27**: Free energy calculations (TI, ABMD, umbrella sampling)
- **Chapter 26**: Enhanced sampling (SGLD, accelerated MD, GaMD)
- **Chapter 32**: CryoEM and EMAP restraints
- **Chapter 31**: NMR refinement