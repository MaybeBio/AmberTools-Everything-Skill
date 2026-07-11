# Chapter 1: Installation & Quickstart

## Core Commands & Syntax

Amber 26 uses CMake as its build system. AmberTools and Amber (pmemd) are separate distributions with independent source trees and install directories.

### Extract the source

```bash
cd /home/xxxx
tar xvfj ambertools26.tar.bz2    # extracts into ambertools26_src/
tar xvfj pmemd26.tar.bz2         # extracts into pmemd26_src/
```

### Apply patches

```bash
cd ambertools26_src
./update_amber --help
./update_amber --update          # applies any available patches

cd ../pmemd26_src
./update_pmemd --help
./update_pmemd --update
```

### Configure and build (serial)

```bash
cd ambertools26_src/build
./run_cmake                      # edit run_cmake first to adjust options
make install
```

### Source the environment

```bash
source /home/xxxx/ambertools26/amber.sh   # bash, zsh, ksh
# or
source /home/xxxx/pmemd26/amber.sh
```

This sets `AMBERHOME` (or `PMEMDHOME`) and adds `$AMBERHOME/bin` to your `PATH`. Add it to `~/.bashrc` or `~/.zshrc` for persistence.

### Conda install (pre-compiled AmberTools only)

```bash
conda install ambertools
```

Users choosing this route only need to download and compile `pmemd26.tar.bz2`. Note: remove any existing miniconda from your PATH while building Amber from source (`conda deactivate`).

### Run tests

```bash
cd $AMBERHOME
make test.serial
```

Check `*.dif` files in `$AMBERHOME/AmberTools/test` or `$AMBERHOME/test` for failures. Differences limited to round-off in the last printed digit are acceptable. All failure messages are collected in `$AMBERHOME/logs`.

## Key Namelists / Input Files

### run_cmake script options

The `run_cmake` script in `ambertools26_src/build/` can be edited to set CMake flags:

```
-DMPI=TRUE                  # Enable MPI parallelism
-DCUDA=TRUE                 # Enable CUDA GPU acceleration
-DOPENMP=TRUE               # Enable OpenMP threading
-DHIP=ON                    # Enable AMD GPU (ROCm) support
-DNCCL=TRUE                 # Enable NVIDIA NCCL multi-GPU comms
-DINSTALL_TESTS=TRUE        # Install test suite
-DCMAKE_INSTALL_PREFIX=/path/to/install  # Install destination
```

### CMake variables for CUDA

```bash
-DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-11.0
```

### CMake variables for AMD HIP

```bash
-DHIP_TOOLKIT_ROOT_DIR=/opt/rocm/
-DHIP_WARP64=OFF            # Required for RDNA2 GPUs (RX 6900 XT, etc.)
```

## Common Workflows

### Full serial build (AmberTools)

```bash
cd ambertools26_src/build
./run_cmake
make install
source /home/xxxx/ambertools26/amber.sh
make test.serial
```

### Full serial build (pmemd)

```bash
cd pmemd26_src/build
./run_cmake
make install
source /home/xxxx/pmemd26/amber.sh
make test.serial
```

### MPI parallel build

```bash
cd ambertools26_src/build
# Edit run_cmake: set -DMPI=TRUE
./run_cmake
make install
export DO_PARALLEL="mpirun -np 2"
source /home/xxxx/ambertools26/amber.sh
make test.parallel

# Test with more threads
export DO_PARALLEL="mpirun -np 4"
make test.parallel
```

MPI builds create `.MPI` suffix executables alongside the serial versions. Amber uses CMake's `FindMPI` module (not compiler wrappers). To specify a particular MPI implementation:

```bash
-DMPI_C_COMPILER=/path/to/mpicc
-DMPI_CXX_COMPILER=/path/to/mpicxx
-DMPI_Fortran_COMPILER=/path/to/mpif90
```

### CUDA GPU build

```bash
# Edit run_cmake: set -DCUDA=TRUE
./run_cmake
make install
```

Builds CUDA versions of pmemd, cpptraj, mdgx, pbsa, and QUICK. If MPI is also enabled, MPI+CUDA variants are built as well.

### AMD GPU build (HIP/ROCm)

```bash
# Edit run_cmake: set -DHIP=ON
./run_cmake
make install
# Test:
cd $PMEMDHOME/test
make test.hip.serial
make test.hip.parallel   # requires GPU-aware MPI
```

Requires ROCm 6.0+. Supports MI100, MI210, MI250, MI300, and RDNA2 GPUs. A `compile_with_hip.sh` script in `pmemd26_src/` provides detailed instructions.

## Reference Tables

### Build flags reference

| Flag | Purpose | Example |
|------|---------|---------|
| `-DMPI=TRUE` | MPI parallelization | `-DMPI=TRUE` |
| `-DCUDA=TRUE` | NVIDIA GPU support | `-DCUDA=TRUE` |
| `-DOPENMP=TRUE` | OpenMP threading | `-DOPENMP=TRUE` |
| `-DHIP=ON` | AMD GPU support | `-DHIP=ON` |
| `-DNCCL=TRUE` | Multi-GPU communication | `-DNCCL=TRUE` |
| `-DINSTALL_TESTS=TRUE` | Include test suite | `-DINSTALL_TESTS=TRUE` |
| `-DCUDA_TOOLKIT_ROOT_DIR` | CUDA install path | `-DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-11.0` |
| `-DHIP_TOOLKIT_ROOT_DIR` | ROCm install path | `-DHIP_TOOLKIT_ROOT_DIR=/opt/rocm/` |
| `-DHIP_WARP64=OFF` | RDNA2 GPU support | `-DHIP_WARP64=OFF` |

### Build output suffixes

| Configuration | Suffix |
|---------------|--------|
| Serial | (none) |
| MPI | `.MPI` |
| OpenMP | `.OMP` |
| CUDA | `.cuda` |
| CUDA + MPI | `.cuda.MPI` |

### Key environment variables

| Variable | Set by | Purpose |
|----------|--------|---------|
| `AMBERHOME` | `source amber.sh` (from ambertools) | AmberTools install root |
| `PMEMDHOME` | `source amber.sh` (from pmemd) | pmemd install root |
| `DO_PARALLEL` | User | MPI launcher for tests |
| `PATH` | `source amber.sh` | Adds `$AMBERHOME/bin` |

## Worked Example

Building AmberTools 26 with MPI and CUDA support on a Linux cluster:

```bash
# 1. Extract and patch
tar xvfj ambertools26.tar.bz2
cd ambertools26_src
./update_amber --update

# 2. Configure
cd build
# Edit run_cmake: set -DMPI=TRUE -DCUDA=TRUE -DINSTALL_TESTS=TRUE
./run_cmake

# 3. Build
make install

# 4. Set up environment
source /home/xxxx/ambertools26/amber.sh

# 5. Test serial
make test.serial

# 6. Test parallel
export DO_PARALLEL="mpirun -np 4"
make test.parallel

# 7. Verify CUDA executables exist
ls $AMBERHOME/bin/pmemd.cuda*
```

## Key Takeaways

1. **AmberTools and pmemd are separate.** Each has its own tarball, build directory, and `amber.sh` resource file. `AMBERHOME` is for AmberTools, `PMEMDHOME` is for pmemd.
2. **The build directory must differ from the source.** `-DCMAKE_INSTALL_PREFIX` must NOT point to the `ambertools26_src` directory.
3. **Deactivate conda before building.** Run `conda deactivate` to avoid conflicts. You cannot use cmake from a conda distribution.
4. **Test after every build.** Always run `make test.serial` before `make test.parallel`. Round-off differences in the last digit are normal.
5. **Water model must be explicitly loaded.** In Amber 26, `leaprc.water.xxxx` is required for explicit solvent simulations -- TIP3P is no longer loaded by default.

## Connects To

- Chapter 2: Force Fields -- which force fields to load after installation
- Chapter 3: LEaP System Building -- building systems with your installed AmberTools
- Section 24.6.4 (pmemd.cuda) -- GPU-specific compilation
- [ambermd.org/Installation.php](https://ambermd.org/Installation.php) -- OS-specific prerequisites
- [ambermd.org/pmwiki/index.php/Main/CMake-Quick-Start](https://ambermd.org/pmwiki/index.php/Main/CMake-Quick-Start)