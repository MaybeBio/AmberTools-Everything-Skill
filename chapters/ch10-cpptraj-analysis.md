# Chapter 10: CPPTRAJ Trajectory Analysis

## Core Commands & Syntax

### Running CPPTRAJ
```bash
# Batch mode (recommended)
cpptraj -p system.prmtop -i analysis.in

# Interactive mode
cpptraj -p system.prmtop
# or simply:
cpptraj
```

### Command-line flags
| Flag | Purpose |
|------|---------|
| `-p <Top0>` | Load topology file(s) |
| `-i <Input0>` | Read input from file(s) |
| `-y <trajin>` | Read trajectory |
| `-x <trajout>` | Write trajectory |
| `-c <reference>` | Reference coordinates |
| `-o <output>` | Redirect STDOUT to file |
| `--log <logfile>` | Log interactive commands |
| `--rng <type>` | RNG: marsaglia, stdlib, mt, pcg32, xo128 |

### Command categories
- **Immediate**: `parm`, `reference`, `readdata`, `writedata`, `list`, `calc`, `clear`, `go`, `run`
- **Queued** (executed at `run`/`go`): `trajin`, `trajout`, `rms`, `rmsd`, `hbond`, `clustering`, `strip`, `autoimage`, `center`, `image`, `distance`, `angle`, `dihedral`, `radgyr`, `radial`, `surf`, `vector`, `matrix`, `atomicfluct`, `nativecontacts`, `secstruct`, `lifetime`, `lipidorder`, `density`, `dipole`, `grid`, `volume`, `xtalsymm`, `watershell`

## Key Input File Structure

### Basic analysis template
```cpptraj
parm system.prmtop
reference ref.pdb
trajin prod.nc 1 last 10        # read every 10th frame

# Center and image the trajectory
center :1-300 mass origin
image origin center familiar
autoimage

# Strip water and ions for solute analysis
strip :WAT,:Na+,:Cl- outprefix solute

# RMSD to reference
rmsd reference out rmsd.dat :1-300@CA

# RMSF
atomicfluct out rmsf.dat :1-300@CA

run
writedata
```

### Atom mask syntax
| Mask | Meaning |
|------|---------|
| `:1-10` | Residues 1 through 10 |
| `@CA` | All CA atoms |
| `@/H` | Exclude all hydrogens |
| `:WAT` | All water residues |
| `:Na+` | Sodium ions |
| `@CA,C,N` | Backbone heavy atoms |
| `:1-10@CA` | CA atoms in residues 1-10 |
| `!@H=` | All non-hydrogen atoms |
| `:LIG<@5.0` | Atoms within 5 A of residue LIG |

### For loops in cpptraj
```cpptraj
for i=1;i<10;i++
    trajin md$i.nc
done
# or with masks:
for maskpairs mask1,mask2 in :1-10@CA
    distance d_$mask1 out dist.dat $mask1 $mask2
done
```

## Common Workflows

### Workflow 1: RMSD and RMSF analysis
```cpptraj
parm tz2.parm7
trajin tz2.nc
autoimage
rms first out rmsd.dat :1-12@CA
atomicfluct out rmsf.dat :1-12@CA
run
```

### Workflow 2: Hydrogen bond analysis
```cpptraj
parm system.prmtop
trajin prod.nc
hbond HB series hbt out hbond.dat \
    avgout hbavg.dat \
    solventdonor :WAT \
    solventacceptor :WAT@O
run
```

### Workflow 3: Clustering (hierarchical, k-means, dbscan, DPEAK)
```cpptraj
parm system.prmtop
trajin prod.nc
rms first :1-300@CA

# Hierarchical agglomerative (clusters 5)
cluster C1 \
    hierarchicalagglo clusters 5 \
    rms :1-300@CA \
    sieve 10 \
    out cnumvtime.dat \
    repout rep repfmt pdb

# k-means
cluster C2 \
    kmeans clusters 5 \
    rms :1-300@CA \
    out cnumvtime_kmeans.dat

# DBSCAN
cluster C3 \
    dbscan minpoints 5 epsilon 2.0 \
    rms :1-300@CA \
    out cnumvtime_dbscan.dat

# DPEAK
cluster C4 \
    dpeak rms :1-300@CA \
    out cnumvtime_dpeak.dat
run
```

### Workflow 4: PCA (Principal Component Analysis)
```cpptraj
parm myparm.parm7
trajin mytraj.nc

# Step 1: Average structure
rms first !@H=
average crdset AVG
run

# Step 2: Covariance matrix
rms ref AVG !@H=
matrix covar name MyMatrix !@H=
createcrd CRD1
run

# Step 3: Diagonalize
runanalysis diagmatrix MyMatrix vecs 2 name MyEvecs

# Step 4: Project along eigenvectors
crdaction CRD1 projection evecs MyEvecs !@H= \
    out project.dat beg 1 end 2
```

### Workflow 5: Dihedral PCA (phi/psi)
```cpptraj
parm ../1rrb_vac.prmtop
trajin ../1rrb_vac.mdcrd
multidihedral BB phi psi resrange 2
run
matrix dihcovar dihedrals BB[*] out dihcovar.dat name DIH
diagmatrix DIH vecs 4 out modes.dat name DIHMODES
run
projection evecs DIHMODES out dih.project.dat beg 1 end 4 dihedrals BB[*]
run
```

### Workflow 6: Post-processing (strip, autoimage, center)
```cpptraj
parm solvated.prmtop
trajin prod.nc

# Autoimage for periodic systems
autoimage

# Strip water and ions, write new topology
strip :WAT,:Na+,:Cl- outprefix dry

# Center on protein
center :1-300 mass origin

# Write processed trajectory
trajout processed.nc
run
```

## Reference Tables

### Key Action Commands
| Command | Purpose | Key Options |
|---------|---------|-------------|
| `rms` / `rmsd` | RMS fit + RMSD | `first`, `ref <ref>`, `perres`, `out` |
| `atomicfluct` / `rmsf` | RMS fluctuations | `out`, `byres` |
| `autoimage` | Auto-center in PBC | `origin`, `familiar` |
| `center` | Center coordinates | `mass`, `origin`, `box` |
| `image` | Wrap molecules in PBC | `origin`, `center` |
| `strip` | Remove atoms | `:WAT`, `outprefix` |
| `hbond` | Hydrogen bonds | `series`, `solventdonor`, `solventacceptor` |
| `distance` | Distance between masks | `out` |
| `angle` | Angle between 3 masks | `out` |
| `dihedral` | Dihedral angle | `out` |
| `radgyr` / `rog` | Radius of gyration | `out` |
| `radial` / `rdf` | Radial distribution function | `out`, `spacing`, `maximum` |
| `surf` | Surface area | `out` |
| `nativecontacts` | Native contacts | `mindist`, `out`, `series` |
| `secstruct` | Secondary structure | `out`, `sumout` |
| `dipole` | Dipole moment | `out` |
| `vector` | Vector between masks | `out`, `magnitude` |
| `matrix` | Pairwise matrix | `covar`, `dist`, `idea`, `dihcovar` |
| `multidihedral` | Multiple dihedrals | `phi`, `psi`, `resrange` |
| `lipidorder` | Lipid order parameters | `out` |
| `molsurf` | Molecular surface | `out` |
| `grid` | Grid data set | `out` |
| `volume` | Volume | `out` |
| `xtalsymm` | Crystal symmetry | `out` |
| `watershell` | Water shell | `out` |
| `lifetime` | Lifetime analysis | `out` |
| `timecorr` | Time correlation | `out` |
| `rotdif` | Rotational diffusion | `out` |
| `average` | Average coordinates | `crdset` |
| `createcrd` | Create COORDS data set | `name` |
| `cluster` | Cluster (action) | `hierarchicalagglo`, `kmeans`, `dbscan`, `dpeak` |
| `clustering` | Cluster (analysis) | `clusters`, `sieve`, `rms` |

### Key Analysis Commands
| Command | Purpose |
|---------|---------|
| `modes` | Normal mode analysis |
| `projection` | Project onto eigenvectors |
| `diagmatrix` | Diagonalize matrix |
| `hist` | Histogram |
| `kde` | Kernel density estimation |
| `lifetime` | Lifetime analysis of data |
| `statistics` | Statistical analysis |
| `corr` / `correlationcoe` | Correlation coefficient |
| `crosscorr` | Cross-correlation |
| `curvefit` | Curve fitting |
| `fft` | Fast Fourier Transform |
| `integrate` | Numerical integration |
| `multicurve` | Multi-curve fitting |
| `runningavg` | Running average |
| `spline` | Spline interpolation |
| `wavelet` | Wavelet analysis |
| `ti` | Thermodynamic integration analysis |
| `tica` | Time-lagged ICA |

### Data file commands
| Command | Purpose |
|---------|---------|
| `readdata` | Read data from file |
| `writedata` / `write` | Write data to file |
| `datafile` | Declare data file |
| `dataset` | Modify data set |
| `list` | List data sets |

## Worked Example

### Complete RMSD + clustering analysis of Tz2
```cpptraj
parm tz2.parm7
trajin tz2.nc

# RMSD to first frame
rms first :1-12 out rmsd.dat

# RMSD to average structure
rms first :1-12
average crdset AVG
run
rms ref AVG :1-12 out rmsd_avg.dat

# RMSF
atomicfluct out rmsf.dat :1-12@CA

# Hierarchical clustering
cluster C1 \
    hierarchicalagglo clusters 5 \
    rms :1-12@CA \
    sieve 10 \
    out cnumvtime.dat \
    repout rep repfmt pdb

run
```

## Key Takeaways

1. **Always autoimage** before analysis of periodic systems to fix molecules that cross the periodic boundary.
2. **Use `strip` to remove water and ions** before solute analysis; this reduces file sizes and speeds up calculations.
3. **The `sieve` option** in clustering drastically reduces memory usage by using every Nth frame for initial clustering.
4. **For PCA**, you must align to a reference structure first to remove global translation/rotation before computing the covariance matrix.
5. **Use `outprefix` with `strip`** to write a new topology file matched to the stripped atoms, enabling further analysis.

## Connects To

- **Chapter 9 (Production MD)**: CPPTRAJ processes trajectories produced by pmemd.cuda
- **Chapter 11 (Free Energy TI)**: `dv/dlambda` analysis via `ti` command and `alchemical_analysis`
- **Chapter 12 (MMPBSA)**: MMPBSA.py uses cpptraj/molsurf internally for surface area calculations
- **Chapter 8 (Equilibration)**: Pre-analysis checks of density, temperature, and pressure evolution