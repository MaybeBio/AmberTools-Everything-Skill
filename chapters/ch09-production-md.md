# Chapter 9: Production MD with pmemd

## Core Commands & Syntax

### Single-GPU Production MD
```bash
pmemd.cuda -O -i md.in -o md.out -p system.prmtop -c eq.rst7 \
    -r md.rst7 -x md.nc -inf md.info
```

### Command-line flags
| Flag | Purpose |
|------|---------|
| `-i` | Input file (mdin) |
| `-o` | Output file (mdout) |
| `-p` | Topology file (prmtop) |
| `-c` | Input coordinates (inpcrd/rst7) |
| `-r` | Output restart file |
| `-x` | Output trajectory file |
| `-inf` | MD info file with performance data |
| `-O` | Overwrite existing output files |

### GPU Selection
```bash
export CUDA_VISIBLE_DEVICES=1    # Use GPU 1
unset CUDA_VISIBLE_DEVICES       # Use all GPUs
```

### Multi-GPU (pmemd.cuda.MPI)
```bash
mpirun -np 4 pmemd.cuda.MPI -O -i md.in -o md.out \
    -p system.prmtop -c eq.rst7 -r md.rst7 -x md.nc
```
One MPI process per GPU. Use `CUDA_VISIBLE_DEVICES` to select specific GPUs per node.

### MVAPICH2-GDR (for improved multi-GPU scaling)
```bash
# Build with:
# -DMVAPICH2GDR_GPU_DIRECT_COMM=TRUE
mpirun -np 4 pmemd.cuda.MPI -O -i md.in ...
```
Activated for same-node pmemd.cuda.MPI runs on more than 2 GPUs.

## Key Namelist / Input File

### Production MD input (NPT, 300 K, 1 atm)
```
Production
 &cntrl
  imin=0,        ! MD (not minimization)
  ntx=5,         ! Read coordinates and velocities from NetCDF restart
  irest=1,       ! Restart simulation
  nstlim=5000000,! 5M steps = 10 ns at dt=0.002
  dt=0.002,      ! 2 fs timestep
  ntf=2,         ! Skip SHAKE force calculations
  ntc=2,         ! SHAKE on bonds involving H
  temp0=300.0,   ! Target temperature (K)
  ntpr=1000,     ! Print to mdout every 1000 steps
  ntwx=5000,     ! Write trajectory frame every 5000 steps
  ntwr=50000,    ! Write restart every 50000 steps
  ioutfm=1,      ! NetCDF trajectory (binary, compact)
  ntxo=2,        ! NetCDF restart format
  cut=8.0,       ! Nonbonded cutoff (Angstroms)
  ntb=2,         ! Constant pressure PBC
  ntp=1,         ! Isotropic pressure scaling
  ntt=3,         ! Langevin thermostat
  gamma_ln=2.0,  ! Collision frequency (ps^{-1})
  ig=-1,         ! Random seed from system clock
 /
```

### Key parameter reference
| Variable | Description | Typical Values |
|----------|-------------|----------------|
| `imin` | 0=MD, 1=minimization | 0 |
| `ntx` | 1=coords only, 5=NetCDF coords+v | 1 or 5 |
| `irest` | 0=new run, 1=restart | 0 then 1 |
| `nstlim` | Number of MD steps | 5000000 (10 ns) |
| `dt` | Timestep in ps | 0.002 (2 fs) |
| `ntb` | 0=no PBC, 1=const vol, 2=const press | 1 (NVT) or 2 (NPT) |
| `ntp` | 0=no press ctrl, 1=isotropic, 2=anisotropic | 0 (NVT) or 1 (NPT) |
| `ntt` | Thermostat: 1=Berendsen, 3=Langevin | 3 |
| `ntpr` | Print frequency (steps) | 1000 |
| `ntwx` | Trajectory write frequency | 5000 |
| `ntwr` | Restart write frequency | 50000 |
| `ioutfm` | 0=ASCII, 1=NetCDF trajectory | 1 |
| `ntxo` | 0=ASCII, 1=binary, 2=NetCDF restart | 2 |
| `ntave` | Interval for averaging energies | 0 (off) |
| `cut` | Nonbonded cutoff (Angstroms) | 8.0 |
| `gamma_ln` | Langevin collision frequency | 2.0 |
| `barostat` | 1=Berendsen, 2=Monte Carlo | 1 or 2 |

## Common Workflows

### Workflow 1: Complete production MD pipeline (GPU)
```bash
# 1. Minimization
pmemd.cuda -O -i 01_Min.in -o 01_Min.out -p parm7 -c rst7 \
    -r 01_Min.ncrst -inf 01_Min.info

# 2. Heating (NVT, 0K -> 300K)
pmemd.cuda -O -i 02_Heat.in -o 02_Heat.out -p parm7 \
    -c 01_Min.ncrst -r 02_Heat.ncrst -x 02_Heat.nc -inf 02_Heat.info

# 3. Equilibration (NPT, short)
pmemd.cuda -O -i 03_Equil.in -o 03_Equil.out -p parm7 \
    -c 02_Heat.ncrst -r 03_Equil.ncrst -x 03_Equil.nc

# 4. Production MD (NPT, long)
pmemd.cuda -O -i 04_Prod.in -o 04_Prod.out -p parm7 \
    -c 03_Equil.ncrst -r 04_Prod.ncrst -x 04_Prod.nc
```

### Workflow 2: Independent parallel runs (multiple replicas)
```bash
# Run N independent simulations with different random seeds
for i in $(seq 1 5); do
  sed "s/ig=-1/ig=$RANDOM/" 04_Prod.in > 04_Prod_${i}.in
  pmemd.cuda -O -i 04_Prod_${i}.in -o 04_Prod_${i}.out \
    -p parm7 -c 03_Equil_${i}.ncrst \
    -r 04_Prod_${i}.ncrst -x 04_Prod_${i}.nc &
done
wait
```

### Workflow 3: Multi-GPU single simulation
```bash
# Requires pmemd.cuda.MPI built with NCCL support
export CUDA_VISIBLE_DEVICES=0,1,2,3
mpirun -np 4 pmemd.cuda.MPI -O -i md.in -o md.out \
    -p system.prmtop -c eq.rst7 -r md.rst7 -x md.nc
```

### Monitoring a running simulation
```bash
tail -f 03_Prod.out     # Monitor output
cat 03_Prod.info         # Check performance / ETA
```

### Restarting from a checkpoint
```bash
pmemd.cuda -O -i 04_Prod.in -o 04_Prod_part2.out \
    -p parm7 -c 04_Prod.ncrst -r 04_Prod_part2.ncrst \
    -x 04_Prod_part2.nc
```
Set `irest=1, ntx=5` in the input file. The restart file contains both coordinates and velocities.

## Reference Tables

### NVT vs NPT
| Property | NVT | NPT |
|----------|-----|-----|
| `ntb` | 1 | 2 |
| `ntp` | 0 | 1 (isotropic) |
| Volume | Constant | Fluctuates |
| Pressure | Fluctuates | Constant (target) |
| Use case | Heating, equilibration | Production |

### Trajectory write strategies
| Format | `ioutfm` | `ntxo` | Pros |
|--------|----------|--------|------|
| ASCII | 0 | 0 | Human-readable, portable |
| NetCDF | 1 | 2 | Compact, fast, self-describing |

### Proper equilibration checks
```bash
# Check density convergence
cpptraj -p parm7 <<EOF
trajin 03_Equil.nc
density out density.dat
run
EOF

# Check temperature
cpptraj -p parm7 <<EOF
trajin 03_Equil.nc
temperature out temp.dat
run
EOF
```

## Worked Example

### Alanine dipeptide production MD (NPT, 300 K, 60 ps)
Minimization input:
```
Minimize
 &cntrl
  imin=1, ntx=1, irest=0, maxcyc=2000, ncyc=1000,
  ntpr=100, ntwx=0, cut=8.0,
 /
```

Heating input (NVT, 0K -> 300K over 20 ps):
```
Heat
 &cntrl
  imin=0, ntx=1, irest=0, nstlim=10000, dt=0.002,
  ntf=2, ntc=2, tempi=0.0, temp0=300.0,
  ntpr=100, ntwx=100, cut=8.0,
  ntb=1, ntp=0, ntt=3, gamma_ln=2.0,
  nmropt=1, ig=-1,
 /
&wt type='TEMP0', istep1=0, istep2=9000, value1=0.0, value2=300.0 /
&wt type='TEMP0', istep1=9001, istep2=10000, value1=300.0, value2=300.0 /
&wt type='END' /
```

Production input (NPT, 60 ps):
```
Production
 &cntrl
  imin=0, ntx=5, irest=1, nstlim=30000, dt=0.002,
  ntf=2, ntc=2, temp0=300.0, ntpr=100, ntwx=100,
  cut=8.0, ntb=2, ntp=1, ntt=3, barostat=1,
  gamma_ln=2.0, ig=-1,
 /
```

Run sequence:
```bash
sander -O -i 01_Min.in -o 01_Min.out -p parm7 -c rst7 -r 01_Min.ncrst
sander -O -i 02_Heat.in -o 02_Heat.out -p parm7 -c 01_Min.ncrst -r 02_Heat.ncrst -x 02_Heat.nc
pmemd -O -i 03_Prod.in -o 03_Prod.out -p parm7 -c 02_Heat.ncrst -r 03_Prod.ncrst -x 03_Prod.nc
```

## Key Takeaways

1. **Always use NetCDF format** (`ioutfm=1, ntxo=2`) for trajectories and restarts -- it is compact, fast, and self-describing.
2. **Write restart files frequently** (`ntwr=50000` or so) to enable recovery from crashes. Checkpoints are essential.
3. **Use `pmemd.cuda`** for GPU-accelerated MD; use `pmemd.cuda.MPI` only when you need multi-GPU for a single simulation or replica exchange.
4. **Monitor equilibration** by checking density, temperature, and pressure convergence before starting production.
5. **For independent replicas**, run separate `pmemd.cuda` processes with different random seeds rather than using MPI -- this is simpler and more efficient for embarrassingly parallel workloads.

## Connects To

- **Chapter 8 (Equilibration)**: Minimization and heating steps precede production MD
- **Chapter 10 (CPPTRAJ Analysis)**: All trajectory analysis including RMSD, RMSF, autoimage, clustering
- **Chapter 11 (Free Energy TI)**: TI production runs use similar pmemd.cuda flags but add `icfe=1`, `clambda`, `scalpha`, `scbeta`
- **Chapter 12 (MMPBSA)**: MMPBSA.py post-processes production trajectories


# Chapter 37: pmemd -- Performance-Optimized MD Engine

## Command Line

```bash
pmemd         -O -i mdin -o mdout -p prmtop -c inpcrd -r rst7 -x mdcrd -inf mdinfo
pmemd.cuda    -O -i mdin -o mdout -p prmtop -c inpcrd -r rst7 -x mdcrd -inf mdinfo
pmemd.cuda.MPI -O -i mdin -o mdout -p prmtop -c inpcrd -r rst7 -x mdcrd -inf mdinfo
```

## pmemd vs sander Feature Comparison

| Feature | sander | pmemd | pmemd.cuda |
|---------|--------|-------|------------|
| Minimization | Full (XMIN, LMOD, NEB) | Standard | Standard |
| MD (NVE/NVT/NPT) | Yes | Yes | Yes |
| TI | Yes (multisander) | Yes (single prmtop) | Yes |
| Softcore TI | Yes | Yes | Yes (smoothstep) |
| QM/MM | Yes | No | No |
| NMR restraints | Yes | No | No |
| Targeted MD | Yes | No | No |
| NEB | Yes | Yes | Yes |
| LMOD | Yes | No | No |
| GaMD | No (sander only) | Yes | Yes |
| REMD | Via multisander | Yes | Yes |
| Constant pH | Via sander | Yes | Yes |
| NFE (Umbrella/SMD/ABMD) | Yes | Yes | Yes |
| GB implicit solvent | Yes | Yes | Yes |
| PB implicit solvent | Yes | No | No |

## GPU Acceleration (pmemd.cuda)

### GPU Selection
```bash
# Single GPU
export CUDA_VISIBLE_DEVICES=0
pmemd.cuda -O -i md.in -o md.out -p prmtop -c inpcrd

# Multi-GPU (one process per GPU)
export CUDA_VISIBLE_DEVICES=0,1,2,3
mpirun -np 4 pmemd.cuda.MPI -O -i md.in -o md.out -p prmtop -c inpcrd
```

### AMD GPU Support (HIP/ROCm)
```bash
# Build with HIP support
cmake .. -DHIP=ON
# Supported: MI100, MI210, MI250, MI300, RDNA2 (RX 6900 XT)
# For RDNA2: -DHIP_WARP64=OFF
# Run tests: make test.hip.serial
```

### GPU Performance Options

| Variable | Type | Description |
|----------|------|-------------|
| `CUDA_VISIBLE_DEVICES` | env | Select which GPU(s) to use |
| `gpu_id` | int | GPU device ID |
| `peer_to_peer` | int | Enable GPU P2P communication |
| `gpu_direct` | int | Enable GPU Direct RDMA |
| `NCCL` | build | NVIDIA Collective Communications Library for multi-GPU |
| `MVAPICH2-GDR` | build | GPU-to-GPU direct communication (84% improvement) |

## Performance Optimization Guidelines (from manual)

| # | Rule | Detail |
|---|------|--------|
| 1 | Use GPU-suitable GBSA (igb=8) | Not igb=2 or igb=5 on GPU |
| 2 | Avoid `ntave/=0` | Running averages calculation hurts performance |
| 3 | Avoid NPT (ntb=2) when not needed | Monte Carlo barostat is expensive |
| 4 | Use `gbsa=3` for GB on GPU | Best GPU performance |
| 5 | Use Berendsen (ntt=1) or Andersen (ntt=2) over Langevin | Lower overhead |
| 6 | Set `netfrc=0` in `&ewald` | Legacy force calculation mode |
| 7 | Small GB systems (<~10K atoms) may be faster on CPU | GPU overhead dominates for tiny systems |
| 8 | Use pmemd.cuda.MPI with NCCL for multi-GPU | Proper MPI+GPU build |
| 9 | Multi-GPU per single run rarely beneficial | Benchmark before scaling |

## Key Namelist Variables for pmemd

### Performance-Related

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `skinnb` | float | 2.0 | Nonbonded pairlist skin width (Å) |
| `skinnb_lab` | float | -- | Lab-specific skin width |
| `nbtell` | int | -- | Verbose pairlist reporting |
| `vdwmeth` | int | 0 | vdW method: 0=PME direct, 1=PME via FFT |
| `eedmeth` | int | -- | Electrostatic energy decomposition method |
| `netfrc` | int | 1 | Ewald force calculation: 0=legacy, 1=current |
| `ee_type` | int | 0 | Electrostatic interaction type |
| `dsum_tol` | float | 1e-6 | Direct sum tolerance |
| `rsum_tol` | float | 5e-5 | Reciprocal sum tolerance |
| `maxexp` | float | -- | Maximum exponent in Ewald |
| `nbflag` | int | -- | Nonbonded calculation flag |

### PME Grid

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `nfft1` | int | auto | FFT grid X-dimension |
| `nfft2` | int | auto | FFT grid Y-dimension |
| `nfft3` | int | auto | FFT grid Z-dimension |
| `order` | int | 4 | B-spline interpolation order |
| `use_pme` | int | 1 | PME on/off |
| `use_fft` | int | 1 | FFT on/off |

### HMR (Hydrogen Mass Repartitioning) -- pmemd.cuda
When using HMR (parmed with hmassrepartition), set `dt=0.004` and `ntc=2, ntf=2`.

## NFE Toolkit in pmemd

The Nonequilibrium Free Energy toolkit works with pmemd and pmemd.cuda via `infe=1`.

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `infe` | int | 0 | Enable NFE: 0=off, 1=on |
| `nfe_umbrella` | int | 0 | Umbrella sampling mode |
| `nfe_abmd` | int | 0 | Adaptive biasing MD |
| `nfe_bbmd` | int | 0 | Barrier-bypassing MD |
| `nfe_steered` | int | 0 | Steered MD (SMD) |
| `nfe_jar` | int | 0 | Jarzynski sampling |
| `nfe_pulling` | str | '' | Pulling direction |
| `nfe_target` | str | '' | Target coordinate file |
| `nfe_moving` | str | '' | Moving atoms mask |
| `nfe_mark1` | str | '' | First marker atom |
| `nfe_mark2` | str | '' | Second marker atom |
| `nfe_spring` | float | -- | Spring constant |
| `nfe_irate` | int | -- | Rate control |
| `nfe_igathering` | int | -- | Gathering interval |
| `nfe_nwork` | int | -- | Number of work values |
| `nfe_dwork` | float | -- | Work increment |
| `nfe_nstep` | int | -- | Number of steps |
| `nfe_nsteep` | int | -- | Number of steep steps |

## NEB in pmemd (Partial NEB)

```bash
# pmemd NEB requires multipmemd with groupfile
mpirun -np <N*num_images> pmemd.MPI -ng <num_images> -groupfile groupfile
```

### NEB Variables

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ineb` | int | 0 | NEB flag: 0=off, 1=on |
| `skmin` | float | -- | Spring constant minimum |
| `skmax` | float | -- | Spring constant maximum |
| `dsk` | float | -- | Spring constant delta (Eq. 26.5) |
| `neb_nxyz` | int | -- | Number of data items per image |
| `nebfitmask` | str | '*' | Atoms for RMS fit between images |
| `nebmask` | str | '*' | Atoms to decouple via NEB |
| `tgtfitmask` | str | '*' | Atoms for endpoint RMS fit |

## Key Takeaways

1. **pmemd is the production engine** -- use pmemd.cuda for GPU, sander only for QM/MM or NMR
2. **No QM/MM in pmemd** -- that's sander-only territory
3. **pmemd TI uses single prmtop** -- simpler setup than sander multisander mode
4. **GPU scaling tips**: reduce `ntave`, avoid NPT when possible, use `igb=8` on GPU
5. **pmemd.cuda TI performs at ~70% of standard MD** -- budget accordingly
6. **NFE toolkit works in pmemd** -- umbrella sampling, ABMD, SMD, Jarzynski
7. **AMD GPUs supported** -- use HIP build with ROCm stack

## Connects To

- **ch07**: sander reference and basic mdin structure
- **ch07-sander**: sander-specific features (QM/MM, NMR, LMOD)
- **ch09**: Production MD and parallel runs
- **ch11**: Free energy TI
- **ch13**: NFE umbrella sampling and SMD
- **ch14**: Enhanced sampling (GaMD, REMD)
- **ch34**: Free energy detailed methodology