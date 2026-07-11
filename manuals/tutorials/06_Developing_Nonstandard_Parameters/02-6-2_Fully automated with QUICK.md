

# Automated RESP Charge Derivation with PyPE_RESP and QUICK
 This tutorial is compatible with AmberTools25 or later versions.
 ## Learning Outcomes
 - Learn how to derive RESP charges starting from molecular coordinates alone, without pre-computed QM data.
- Understand how PyPE_RESP integrates with the QUICK quantum chemistry program for geometry optimization and ESP calculation.
- Learn about multi-orientation RESP fitting and why it produces more robust, orientation-independent charges.
- Understand the QUICK-related keywords in the PyPE_RESP input file for controlling the QM calculation.
 ## Introduction
 The [previous PyPE_RESP tutorial](https://ambermd.org/tutorials/basic/tutorial22/index.php) showed how to derive RESP charges from pre-computed Gaussian or GAMESS output. In many situations, however, the user may only have molecular coordinates (in XYZ or MOL2 format) and would prefer not to set up and run separate QM calculations. PyPE_RESP can handle the entire workflow — from coordinates to RESP charges — by interfacing with the [QUICK](https://quick-docs.readthedocs.io/) quantum chemistry program that ships with AmberTools.
 When the input format is set to `xyz` or `mol2`, PyPE_RESP automatically invokes QUICK to:
 1. Optionally optimize the molecular geometry.
2. Calculate the electrostatic potential (ESP) on a grid surrounding each molecule.
3. Generate multiple molecular orientations and compute ESP for each, producing orientation-independent RESP charges.
 The result is a fully automated pipeline: coordinates in, RESP charges out.
 For background on the RESP algorithm, see the [PyRESP tutorial](https://ambermd.org/tutorials/basic/tutorial19/index.php). For details on PyPE_RESP input file syntax, constraints, and multi-molecule fits, see the [PyPE_RESP tutorial](https://ambermd.org/tutorials/basic/tutorial22/index.php). The [PyPE_RESP paper](https://doi.org/10.1021/acs.jcim.5c00041) describes the program in detail.
 ## Prerequisites
 You need molecular coordinates in XYZ or MOL2 format. QUICK must be available in your `$PATH` (it is included with AmberTools).
 ## Process
 1. ### Preparing Coordinate Files
 Create a working directory and place your coordinate files in it. For this tutorial, we use ethanol as an example. You can use either XYZ or MOL2 format.
 An XYZ file ( `ethanol.xyz`) looks like:
 ```
9
ethanol
C    -0.7516    0.0198   -0.0346
C     0.7405   -0.0804    0.1245
O     1.1936    1.2032   -0.2420
H    -1.1663   -0.9497    0.2424
H    -1.0507    0.2091   -1.0684
H    -1.1815    0.7988    0.5960
H     1.0044   -0.2801    1.1677
H     1.2219   -0.8456   -0.4922
H     2.1499    1.1868   -0.1438

```
 ```
mkdir ethanol_quick
cp ethanol.xyz ethanol_quick/
cd ethanol_quick

```
 If you have multiple conformers, place all of them in the working directory. PyPE_RESP will detect them by filename pattern matching, just as with Gaussian output files.
2. ### A Minimal QUICK-Based Input File
 Create a file called `ethanol_quick.in`:
 ```
# PyPE_RESP input: automated RESP via QUICK
# Starting from XYZ coordinates

MOL_1 = ethanol

FORMAT    = xyz

```
 Setting `FORMAT = xyz` (or `mol2`) tells PyPE_RESP that the input files contain coordinates only and that QUICK should be used for the QM calculations. As in the [previous tutorial](https://ambermd.org/tutorials/basic/tutorial22/index.php), every keyword that is not set takes its default value — in particular `WORK_DIR` defaults to the directory in which `pype-resp.py` is executed, and `ATOM_TYPE` defaults to GAFF2, so neither needs to appear in a minimal input run from inside the working directory. Run it with:
 ```
pype-resp.py -i ethanol_quick.in

```
 With default settings, PyPE_RESP will:
 - Read the ethanol coordinates from the XYZ file.
- Calculate the ESP at the HF/6-31G* level of theory with a grid spacing of 0.25 Å.
- Generate 6 molecular orientations and compute ESP for each.
- Combine all orientations and run the two-stage RESP fit.
- Write mol2 files with fitted charges and a charge summary.
3. ### Geometry Optimization with QUICK
 If your starting coordinates are approximate (e.g., generated from a SMILES string or taken from an unoptimized structure), you should optimize the geometry before computing the ESP. Enable this with the `GEOM_OPT` keyword:
 ```
# PyPE_RESP input with geometry optimization

MOL_1 = ethanol

FORMAT     = xyz

# Optimize geometry before ESP calculation
GEOM_OPT   = True
GEOM_FUNC  = B3LYP
GEOM_BASIS = 6-31G*

```
 The `GEOM_FUNC` and `GEOM_BASIS` keywords control the level of theory for the optimization step. They default to HF/6-31G* if not specified. After optimization, PyPE_RESP uses the optimized coordinates for the subsequent ESP calculation.
 Both the optimization and the ESP step accept any basis set and functional supported by QUICK:
 - **Basis sets:** STO-3G, 3-21G, 6-31G, 6-31G*, 6-31G**, 6-311G*, 6-311+G*, 6-31++G**, cc-pVDZ, cc-pVTZ, aug-cc-pVDZ, aug-cc-pVTZ, def2-SVP, def2-TZVP, and others.
- **Functionals:** HF, B3LYP, BLYP, PBE, PBE0, BP86, CAM-B3LYP, and many more (including LIBXC functionals).
4. ### Controlling the ESP Calculation
 The quality of RESP charges depends on the ESP calculation. PyPE_RESP provides several keywords to control this step:
 ```
# Fine-grained control over the QUICK ESP calculation

MOL_1 = ethanol

FORMAT    = xyz

GEOM_OPT  = True
GEOM_FUNC = B3LYP
GEOM_BASIS = 6-31G**

# ESP calculation settings
ESP_FUNC       = HF             # Functional (default: HF)
ESP_BASIS      = 6-31G*         # Basis set (default: 6-31G*)
ESP_GRID_SPACE = 0.25           # Grid spacing in Angstroms (default: 0.25)

```
 A smaller `ESP_GRID_SPACE` produces a denser grid with more ESP points, which can improve fitting quality at the cost of longer computation time. The default of 0.25 Å is a good balance for most molecules.
 As noted above, `ESP_FUNC` and `ESP_BASIS` accept the same wide range of QUICK-supported basis sets and functionals as the geometry-optimization step.
5. ### Multi-Orientation RESP Fitting
 A well-known limitation of RESP fitting is that the resulting charges depend on the orientation of the molecule relative to the ESP grid. PyPE_RESP addresses this by generating multiple orientations of each molecule and fitting all of them simultaneously. This is controlled by the `NUM_ORIENTATIONS` keyword:
 ```
NUM_ORIENTATIONS = 6    # Default: 6 (range: 1-6)

```
 With the default value of 6, PyPE_RESP generates the original orientation plus 5 additional rotations:
 OrientationRotation 1Original (no rotation) 290° around X-axis 390° around Y-axis 490° around Z-axis 590° around X-axis + 90° around Z-axis 690° around Y-axis + 90° around Z-axis Each orientation generates its own ESP grid, and all grids are combined into a single RESP fit. The resulting charges are more robust and less sensitive to the initial molecular orientation. Setting `NUM_ORIENTATIONS = 1` disables multi-orientation sampling.
 **Note:** Multi-orientation sampling is only available when using QUICK for ESP generation (i.e., `FORMAT = xyz` or `FORMAT = mol2`). When reading pre-computed Gaussian or GAMESS output, the ESP grid is fixed by the QM calculation.
6. ### QUICK Version Selection
 QUICK is available in multiple build variants, selected via the `QUICK_VERSION` keyword:
 ```
QUICK_VERSION = serial      # Default

```
 ValueDescription `serial`Single-threaded CPU execution (default) `mpi`MPI-parallel CPU execution `cuda`GPU-accelerated execution `cuda+mpi`Multi-GPU execution The choice depends on your available hardware. The serial version works everywhere and is sufficient for small to medium-sized molecules. For larger molecules or batch processing, the CUDA-enabled versions can significantly reduce computation time.
 Additional QUICK command-line options can be passed through with:
 ```
QUICK_OPTS = <additional options>

```
7. ### Understanding the Output with QUICK
 When using QUICK, the output directory includes additional files compared to the Gaussian/GAMESS workflow:
 ```
ethanol_quick/                                  # your working directory
  ├── ethanol.xyz                               # your coordinate file
  ├── _PyPE_RESP.log                            # full run log
  └── _PyPE_RESP/                               # PyPE_RESP working directory
      ├── QUICK_FILES/                          # QUICK inputs/outputs (prefix from the input filename)
      │   ├── ethanolxyz_orr0_quick_esp.in      # ESP input  (orientation 1)
      │   ├── ethanolxyz_orr0_quick_esp.out     # QUICK log  (orientation 1)
      │   ├── ethanolxyz_orr0_quick_esp.vdw     # ESP grid   (orientation 1)
      │   ├── ...                               # (orientations 2-6: _orr1 … _orr5)
      │   └── ethanolxyz_orr5_quick_esp.vdw
      │   # with GEOM_OPT = True, ethanolxyz_quick_opt.in/.out also appear here
      ├── ESP_FILES/                            # ESP data converted for PyRESP
      │   ├── ethanol_orr0.esp
      │   ├── ...
      │   └── ethanol_orr5.esp
      ├── RESP_RUN_FILES/                       # generated PyRESP inputs + raw outputs
      │   ├── First_stage.in / .out / .chg / .espout
      │   ├── Second_stage.in / .out / .chg / .espout
      │   └── Concatenated.esp                  # combined ESP data, all orientations
      ├── Mol3_files/
      │   └── ethanol_AFTER_RESP.mol2           # final charges, mol2 form
      └── Run_statistics/
          ├── RESP_charge_summary.log           # fit/charge summary (text)
          └── RESP_charge_summary.csv           # fit/charge summary (CSV)

```
 The **QUICK_FILES/** directory contains the QUICK input ( `.in`), QUICK log ( `.out`), and ESP grid ( `.vdw`) for each orientation; the filename prefix (here `ethanolxyz`) is derived from the coordinate file. If geometry optimization was enabled, the optimized coordinates are found in the additional `*_quick_opt.out` file. PyPE_RESP converts the `.vdw` grids into the per-orientation `.esp` files under **ESP_FILES/**, which are concatenated and fed to the two-stage RESP fit in **RESP_RUN_FILES/**. As in the Gaussian workflow, the final per-atom charges are written to **Mol3_files/ethanol_AFTER_RESP.mol2** (and `Second_stage.chg`), while **Run_statistics/** holds the fit summary (its filename follows `OUTPUT_NAME`, defaulting to `RESP_charge_summary`).
 **Downloads:** the complete worked run — coordinates, input, and the full `_PyPE_RESP/` tree (including the QUICK and ESP grid files) — is packaged as [ethanol_quick.tar.gz](https://ambermd.org/tutorials/basic/tutorial23/include/ethanol_quick.tar.gz). Individually: [ethanol.xyz](https://ambermd.org/tutorials/basic/tutorial23/include/ethanol.xyz), [ethanol_quick.in](https://ambermd.org/tutorials/basic/tutorial23/include/ethanol_quick.in), final charges [ethanol_AFTER_RESP.mol2](https://ambermd.org/tutorials/basic/tutorial23/include/ethanol_AFTER_RESP.mol2) and [Second_stage.chg](https://ambermd.org/tutorials/basic/tutorial23/include/Second_stage.chg), fit output [Second_stage.out](https://ambermd.org/tutorials/basic/tutorial23/include/Second_stage.out), and summary [RESP_charge_summary.log](https://ambermd.org/tutorials/basic/tutorial23/include/RESP_charge_summary.log) / [RESP_charge_summary.csv](https://ambermd.org/tutorials/basic/tutorial23/include/RESP_charge_summary.csv).
 **The fitted charges.** This run constrains the methyl group ( `MOL_1_CON_1 = C(H)(H)(H) : 0.0`), so the three methyl hydrogens and their carbon sum to zero. The final charges in `ethanol_AFTER_RESP.mol2` are:
 ```
@<TRIPOS>ATOM
      1 C        -0.75160    0.01980   -0.03460 c3     1 MOL     -0.047043705
      2 C         0.74050   -0.08040    0.12450 c3     1 MOL      0.379574960
      3 O         1.19360    1.20320   -0.24200 oh     1 MOL     -0.685657070
      4 H        -1.16630   -0.94970    0.24240 hc     1 MOL      0.015681235
      5 H        -1.05070    0.20910   -1.06840 hc     1 MOL      0.015681235
      6 H        -1.18150    0.79880    0.59600 hc     1 MOL      0.015681235
      7 H         1.00440   -0.28010    1.16770 h1     1 MOL     -0.051943611
      8 H         1.22190   -0.84560   -0.49220 h1     1 MOL     -0.051943611
      9 H         2.14990    1.18680   -0.14380 ho     1 MOL      0.409969330

```
 Atoms 1, 4, 5 and 6 (the methyl carbon and its three hydrogens) sum to zero as constrained, the hydroxyl oxygen (atom 3) carries the expected large negative charge (−0.686) and the hydroxyl hydrogen (atom 9) the large positive charge (+0.410). The corresponding `RESP_charge_summary.log` reports a relative RMSE of about 0.17 and confirms the constraint group sums to 0.0.
8. ### A Complete Example
 Below is a complete input file that uses all the QUICK-related features discussed in this tutorial:
 ```
# PyPE_RESP input: full automated workflow with QUICK
# Ethanol from XYZ coordinates, with optimization and multi-orientation RESP

MOL_1 = ethanol

FORMAT    = xyz

# Geometry optimization
GEOM_OPT   = True
GEOM_FUNC  = B3LYP
GEOM_BASIS = 6-31G*

# ESP calculation
ESP_FUNC        = HF
ESP_BASIS       = 6-31G*
ESP_GRID_SPACE  = 0.25

# Multi-orientation fitting (default is 6)
NUM_ORIENTATIONS = 6

# QUICK execution
QUICK_VERSION = serial

# Constraints (optional)
MOL_1_CON_1 = C(H)(H)(H) : 0.0

# General
OUTPUT_NAME = ethanol_charges
VERBOSE     = True

```
 ```
pype-resp.py -i ethanol_quick.in

```
 This single command will optimize the ethanol geometry with B3LYP/6-31G*, calculate the ESP at 6 orientations with HF/6-31G*, run the two-stage RESP fit with a methyl group constrained to neutral, and produce mol2 files with the final charges.
 The downloadable [ethanol_quick.in](https://ambermd.org/tutorials/basic/tutorial23/include/ethanol_quick.in) (whose output is shown in the [output section](#understanding-quick-output) above) is exactly this file with the optional blocks left commented out — so the provided run uses the defaults (no geometry optimization, HF/6-31G* ESP, 6 orientations, default `OUTPUT_NAME`) and keeps only the methyl constraint. Uncomment those blocks to reproduce the fully featured run shown here.
9. ### Parallel Execution
 When fitting charges for multiple molecules or conformers, PyPE_RESP can run QUICK calculations in parallel using the `NCPUS` keyword:
 ```
NCPUS = 4    # Run up to 4 QUICK jobs simultaneously

```
 This is particularly beneficial when you have many conformers or molecules, as the ESP calculations for each are independent and can proceed concurrently.
 ## Summary
 PyPE_RESP combined with QUICK provides a fully automated path from molecular coordinates to RESP charges. The key advantages over the traditional workflow are:
 - **No external QM software required** — QUICK is included with AmberTools.
- **Single input file** — coordinates, QM settings, constraints, and RESP options are all defined in one place.
- **Multi-orientation fitting** — produces charges that are less sensitive to molecular orientation.
- **Automated pipeline** — geometry optimization, ESP calculation, format conversion, and two-stage RESP fitting are handled in a single command.
 ## Next Steps
 - For using PyPE_RESP with pre-computed Gaussian or GAMESS output, see [RESP Charge Derivation with PyPE_RESP](https://ambermd.org/tutorials/basic/tutorial22/index.php).
- For background on the RESP algorithm and the underlying py_resp.py program, see the [PyRESP tutorial](https://ambermd.org/tutorials/basic/tutorial19/index.php).
- For details on using the fitted charges in a simulation, see the [Antechamber and GAFF tutorial](https://ambermd.org/tutorials/basic/tutorial4b/index.php).
 by Marco Lapsien

