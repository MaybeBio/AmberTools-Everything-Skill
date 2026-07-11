# Chapter 29: Torch PBSA

## Overview

Torch PBSA is a dedicated runtime within the PBSA framework that uses **LibTorch** (the C++ distribution of PyTorch) to accelerate performance-critical components of Poisson-Boltzmann calculations. It provides:

1. **Machine-learned solvent-excluded surface (MLSES) models** (GENIUSES, Con2SES) for fast SES construction
2. **AmberTorchPB** -- a tensor-based PB linear system solver with multiple iterative solvers and preconditioners

## Core Commands & Syntax

### PBSA invocation with Torch PBSA
```bash
pbsa -i mdin -o mdout -p prmtop -c inpcrd
pbsa.cuda -i mdin -o mdout -p prmtop -c inpcrd
```

### LibTorch Installation

LibTorch is enabled via CMake built-in mode:

```bash
cmake ... \
  -DLIBTORCH=ON \
  -DLIBTORCH_PLATFORM_OPT=cu126
```

Two PyTorch versions are supported:
- PyTorch 2.10.0 for CUDA 12.6+ (platform IDs: cu126, cu128, cu129, cu130)
- PyTorch 2.5.1 for CUDA 11.8-12.5 (platform IDs: cu118, cu121, cu124)

Requires conda (bundled Miniconda via `-DDOWNLOAD_MINICONDA=TRUE` or system conda).

## Key Namelists / Input Files

### `&pb` Namelist -- MLSES Options

| Variable | Description | Default |
|----------|-------------|---------|
| `sasopt` | Surface area model; set to `3` for MLSES | -- |
| `mlses_opt` | MLSES runtime selection | 0 |
| | = 0: Fortran/CUDA GENIUSES on CPU/GPU | |
| | = 1: Torch PBSA GENIUSES on CPU/GPU | |
| | = 2: Torch PBSA Con2SES-2D on CPU/GPU | |
| | = 3: Torch PBSA Con2SES-3D on CPU/GPU | |

### `&pb` Namelist -- AmberTorchPB Options

| Variable | Description | Default |
|----------|-------------|---------|
| `torchpb` | Enable AmberTorchPB solver | 0 |
| | = 0: Standard PBSA solver | |
| | = 1: Enable AmberTorchPB | |
| `torchpb_solvopt` | Iterative solver selection | 0 |
| | = 0: Conjugate gradient (CG) -- recommended | |
| | = 1: Biconjugate gradient (BiCG) | |
| | = 2: GMRES (fewer iterations, higher per-iteration cost) | |
| | = 3: SGD gradient optimization | |
| | = 4: Adam gradient optimization | |
| | = 5: Neural network solver | |
| | = 6: Successive over-relaxation (SOR) | |
| | = 7: Red-black SOR (RB-SOR) -- more parallel-efficient | |
| `torchpb_cuda_opt` | Compute device | 0 |
| | = 0: CPU with LibTorch native multithreading | |
| | = 1: GPU CUDA acceleration (auto-set in pbsa.cuda) | |
| `torchpb_layout_type` | Sparse matrix layout | 0 |
| | = 0: CSR (compressed sparse row) | |
| | = 1: COO (coordinate) | |
| | = 9: Matrix-free stencil (no explicit storage, lower memory) | |
| `torchpb_precision_type` | Floating-point precision | 1 |
| | = 0: float64 (double) | |
| | = 1: float32 (single) -- default for CPU/GPU | |
| | = 2: float16 (half) -- may cause instability | |
| | = 3: bfloat16 (CG only) | |
| `torchpb_preconditioner_opt` | Preconditioner | 0 |
| | = 0: No preconditioning | |
| | = 1: Block Jacobi (high parallel efficiency) | |
| | = 2: Incomplete Cholesky (ICC) (symmetric positive-definite) | |
| | = 3: Algebraic multigrid (AMG) (strongest convergence) | |
| `torchpb_verbose_opt` | Verbosity | 0 |
| | = 0: No output | |
| | = 1: Print residual at each iteration | |
| | = 2: Print detailed debug information | |

## MLSES Models

### GENIUSES (Grid-robust Efficient Neural Interface for Universal Solvent-Excluded Surface)
- Point cloud-based neural network for SES construction
- Dense matrix multiplications -- efficient on CPU and GPU
- Robust to grid spacing variations
- ~95% fidelity vs. classical AMBER SES
- ~26-33x speedup on GPU
- Custom CUDA implementation for additional overhead reduction

### Con2SES
- Convolutional neural network with learnable 2D and 3D kernels
- Explicitly models grid-context interactions (many-body effects)
- Higher accuracy in complex molecular interiors
- Two variants: Con2SES-2D and Con2SES-3D (3D recommended)
- ~99% accuracy, ~28x speedup vs. classical SES
- Sliding voxel window for scalable inference

## Common Workflows

### Single-point solvation with MLSES
```
&cntrl
   ntx=1, imin=1, ipb=2, inp=0
/
&pb
   npbverb=0, istrng=0,
   epsout=80.0, epsin=1.0, dprob=1.4, radiopt=0,
   sasopt=3,           # Enable MLSES
   mlses_opt=3,        # Con2SES-3D via Torch PBSA
   fillratio=1.5, nfocus=1, space=0.95,
   accept=0.000001, maxitn=100000, solvopt=3,
   npbopt=0, bcopt=6,
   eneopt=1, frcopt=0, cutnb=15, cutsa=8, cutfd=7,
/
```

### Single-point solvation with AmberTorchPB
```
&cntrl
   ntx=1, imin=1, ipb=2, inp=0
/
&pb
   npbverb=0, istrng=0,
   epsout=80.0, epsin=1.0, dprob=1.4, radiopt=0,
   sasopt=0,
   fillratio=1.5, nfocus=1, space=0.95,
   accept=0.000001, maxitn=100000, solvopt=3,
   npbopt=0, bcopt=6,
   eneopt=1, frcopt=0, cutnb=15, cutsa=8, cutfd=7,
   torchpb=1, torchpb_solvopt=0,
   torchpb_layout_type=0, torchpb_precision_type=1,
/
```

### GPU-accelerated CG with AMG preconditioner
```
   torchpb=1, torchpb_solvopt=0,
   torchpb_cuda_opt=1,
   torchpb_preconditioner_opt=3,
   torchpb_layout_type=0, torchpb_precision_type=1,
```

## Key Takeaways

1. **Torch PBSA** requires LibTorch installation via CMake (`-DLIBTORCH=ON`); uses conda to manage PyTorch packages.
2. **MLSES** replaces geometric SES algorithms with neural networks: GENIUSES (point-cloud) or Con2SES (CNN). Set `sasopt=3` and `mlses_opt`.
3. **AmberTorchPB** is a unified tensor-based PB solver (`torchpb=1`): supports CG, BiCG, GMRES, SOR, RB-SOR with CSR, COO, or matrix-free layouts.
4. **Precision control**: float32 is default and recommended; bfloat16 only for CG; float64 for high fidelity.
5. **GPU acceleration** up to 28x speedup via `torchpb_cuda_opt=1` or `pbsa.cuda` executable.

## Connects To

- Chapter 6: PBSA (classical Poisson-Boltzmann)
- Chapter 5: GBNSR6 (GB solvation)
- Chapter 35: AI/ML (KMMD, neural network potentials)
- Chapter 31: External library interface (external computation)