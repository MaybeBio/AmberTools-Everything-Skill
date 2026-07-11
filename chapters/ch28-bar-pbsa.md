# Chapter 28: BAR/PBSA -- Post-processing Free Energy Analysis

## Core Commands & Syntax

### bar_pbsa.py Overview

`bar_pbsa.py` is a Python script that automates the preparation of trajectories from the decharging step of alchemical simulations for BAR/PBSA analysis of binding free energies. It is separated into four stages:

```bash
# Stage 1: Strip solvent and ions from explicit solvent trajectories
python bar_pbsa.py strip strip_input.yaml

# Stage 2: Prepare sander PBSA input files
python bar_pbsa.py prep prep_input.yaml

# Stage 3: Run sander in parallel (ligand and complex separately)
python bar_pbsa.py run lig_input.yaml -n 8
python bar_pbsa.py run com_input.yaml -n 8

# Stage 4: Calculate final decharging energies
python bar_pbsa.py calc lig_input.yaml
python bar_pbsa.py calc com_input.yaml
```

### Command Help

```bash
# Get help on bar_pbsa.py
python bar_pbsa.py -h
python bar_pbsa.py --help
```

## Key Namelists / Input Files

### Stage 1: Stripping Input YAML (strip_input.yaml)

Prepares explicit solvent trajectories for sander PBSA by stripping water and ions, concatenating replicate trajectories, and autoimaging/RMSD aligning snapshots to the first frame.

```yaml
dest_path: '1C5X'
complex_paths:
    - '/home/bar_pbsa_demo/1C5X/t1/complex'
    - '/home/bar_pbsa_demo/1C5X/t2/complex'
ligand_paths:
    - '/home/bar_pbsa_demo/1C5X/t1/ligands'
    - '/home/bar_pbsa_demo/1C5X/t2/ligands'
ligand_mask: ':DRG'
ion_decharge: True
last_half_frames: True
stride: 1
```

| Variable | Description |
|----------|-------------|
| dest_path | Location to output stripped trajectories |
| complex_paths | Paths to explicit solvent trajectories for complex decharging. Multiple replicates (t1/t2) should be used to improve convergence of free energy simulations |
| ligand_paths | Paths to explicit solvent trajectories for ligand decharging |
| ligand_mask | Amber mask to select ligand atoms |
| ion_decharge | Boolean. Whether counter-ion decharge should be performed to maintain charge neutrality. The decharged counter-ion is automatically identified by its lack of charge |
| last_half_frames | Boolean. Whether only frames from the last half of the trajectory should be kept. This allows removal of unequilibrated data |
| stride | Step to subsample frames and remove correlated data from sequential snapshots. Higher values minimize the amount of processing necessary |

### Stage 2: Preparation Input YAML (prep_input.yaml)

Writes sander PBSA input files with target radii scaling factors and solute dielectric for post-processing each trajectory.

```yaml
dest_path: '1C5X'
ligand_res: 'DRG'
istrng: 150
epsin: 1.0
radiscale: 1.0
protscale: 1.0
```

| Variable | Description |
|----------|-------------|
| dest_path | Path to stripped trajectories |
| ligand_res | Ligand residue name, used to apply radiscale |
| istrng | Salt concentration in mM |
| epsin | Protein dielectric constant. Higher values more strongly screen charged interactions |
| radiscale | Ligand atoms radii scaling factor |
| protscale | Protein atoms radii scaling factor |

### Stage 3: Run Input YAML (lig_input.yaml / com_input.yaml)

Runs BAR/PBSA sander calculations for neighboring lambdas with multiprocessing. This stage is carried out separately for the complex and ligand paths due to memory limitations.

```yaml
dest_path: '1C5X'
epsin: 1.0
radiscale: 1.0
protscale: 1.0
ligcom: 'complex'
del_traj: True
```

| Variable | Description |
|----------|-------------|
| dest_path | Path to stripped trajectories |
| epsin | Protein dielectric constant. Higher values more strongly screen charged interactions |
| radiscale | Ligand atoms radii scaling factor |
| protscale | Protein atoms radii scaling factor |
| ligcom | Option to indicate whether the setup is for 'ligand' or 'complex' trajectories |
| del_traj | Boolean to delete trajectories after post-processing to save memory space. The trajectories are copied for each new combination of radiscale, protscale, and epsin and are redundant |

### Stage 4: Calculate Input YAML

The input YAML to calculate decharging energies is identical to the input for the run stage (Stage 3).

### Required Folder Structure

The source folders for the complex and ligand trajectories must follow this structure:

```
1C5X/
  t1/
    complex/
      0.000/
        ti.parm7
        ti001.nc
      0.200/
        ti.parm7
        ti001.nc
      0.400/
        ti.parm7
        ti001.nc
      0.600/
        ti.parm7
        ti001.nc
      0.800/
        ti.parm7
        ti001.nc
      1.000/
        ti.parm7
        ti001.nc
    ligands/
      0.000/
        ti.parm7
        ti001.nc
      0.200/
        ti.parm7
        ti001.nc
      ... (same lambda windows)
      1.000/
        ti.parm7
        ti001.nc
  t2/
    complex/...
    ligands/...
```

## Common Workflows

### Full BAR/PBSA Post-processing Pipeline

The BAR/PBSA method addresses the problem that alchemical simulations in standard additive force fields cannot handle electronic polarization effects upon ligand transfer from water to the protein interior. This leads to inaccurate prediction of binding affinities for charged molecules/binding pockets.

**Theory**: The protein dielectric constant in Poisson-Boltzmann models correlates with the strength of electronic polarization in the protein environment. In typical PBSA calculations, the default `epsin` value of 1 is used for Amber force fields, as this is how the force fields are designed and calibrated for MD. With `epsin=1`, no electronic polarization is included, and Coulombic interactions are not screened by electronic polarization, resulting in exaggerated electrostatic interactions. The overestimation can be alleviated by increasing the solute `epsin` value to imitate the effect of electronic polarization that screens electrostatic interactions.

**Calibration Workflow**:

1. Ligand radii are first scaled by the `radiscale` PBSA keyword, optimized to minimize the absolute deviation between PBSA and explicit-solvent electrostatic free energies for the ligand alchemical simulations.

2. Given the optimized `radiscale` value, the protein radii are then scaled by the `protscale` PBSA keyword, optimized to minimize the absolute deviation between PBSA and explicit-solvent electrostatic free energies for the complex alchemical simulations.

3. Following calibration of the atomic radii (radiscale and protscale optimized), the BAR/PBSA method can then be utilized for the investigation of electronic polarization by varying the solute dielectric constant, `epsin`.

### Step-by-Step Command Sequence

```bash
# Follow the 4-stage pipeline:

# Step 1: Strip solvent and ions from trajectories
python bar_pbsa.py strip strip_input.yaml

# Step 2: Prepare sander input files with selected radiscale, protscale, and epsin
python bar_pbsa.py prep prep_input.yaml

# Step 3: Run sander in parallel for ligand and complex separately
# Use -n to select max number of cores for parallel jobs
python bar_pbsa.py run lig_input.yaml -n 8
python bar_pbsa.py run com_input.yaml -n 8

# Step 4: Calculate final decharging energies
# Uses same input files as run step. Log with final energies saved in run directory
python bar_pbsa.py calc lig_input.yaml
python bar_pbsa.py calc com_input.yaml
```

### Running Parameter Scans

To investigate the effect of electronic polarization, run multiple preparation and calculation cycles with different `epsin` values:

```yaml
# prep_epsin2.yaml
dest_path: '1C5X'
ligand_res: 'DRG'
istrng: 150
epsin: 2.0
radiscale: 1.0
protscale: 1.0
```

```yaml
# prep_epsin4.yaml
dest_path: '1C5X'
ligand_res: 'DRG'
istrng: 150
epsin: 4.0
radiscale: 1.0
protscale: 1.0
```

Then run each preparation and calculation cycle:
```bash
python bar_pbsa.py prep prep_epsin2.yaml
python bar_pbsa.py run lig_epsin2.yaml -n 8
python bar_pbsa.py run com_epsin2.yaml -n 8
python bar_pbsa.py calc lig_epsin2.yaml
python bar_pbsa.py calc com_epsin2.yaml
```

## Reference Tables

### BAR/PBSA Python Script Arguments

| Argument | Purpose |
|----------|---------|
| strip | Stage 1: Strip solvent and ions from explicit solvent trajectories |
| prep | Stage 2: Prepare sander PBSA input files with target parameters |
| run | Stage 3: Run parallel sander post-processing trajectories |
| calc | Stage 4: Calculate total decharging energies through BAR |
| -n N | Number of cores for parallel jobs (run stage) |
| -h, --help | Print usage summary |

### Stage 1: Strip Variables

| Variable | Type | Description |
|----------|------|-------------|
| dest_path | string | Location to output stripped trajectories |
| complex_paths | list | Paths to explicit solvent complex decharging trajectories |
| ligand_paths | list | Paths to explicit solvent ligand decharging trajectories |
| ligand_mask | string | Amber mask to select ligand atoms |
| ion_decharge | bool | Perform counter-ion decharge for charge neutrality |
| last_half_frames | bool | Keep only frames from last half of trajectory |
| stride | int | Subsample step to remove correlated data |

### Stage 2: Prep Variables

| Variable | Type | Description |
|----------|------|-------------|
| dest_path | string | Path to stripped trajectories |
| ligand_res | string | Ligand residue name for radiscale application |
| istrng | int | Salt concentration in mM |
| epsin | float | Protein dielectric constant |
| radiscale | float | Ligand atoms radii scaling factor |
| protscale | float | Protein atoms radii scaling factor |

### Stage 3/4: Run/Calc Variables

| Variable | Type | Description |
|----------|------|-------------|
| dest_path | string | Path to stripped trajectories |
| epsin | float | Protein dielectric constant |
| radiscale | float | Ligand atoms radii scaling factor |
| protscale | float | Protein atoms radii scaling factor |
| ligcom | string | 'ligand' or 'complex' |
| del_traj | bool | Delete trajectories after post-processing |

### Key PBSA Parameters for BAR/PBSA

| Parameter | Description | Typical Range |
|-----------|-------------|---------------|
| epsin | Solute dielectric constant | 1.0 - 4.0 (or higher) |
| radiscale | Ligand radii scaling factor | 0.8 - 1.2 |
| protscale | Protein radii scaling factor | 0.8 - 1.2 |
| istrng | Salt concentration (mM) | 0 - 150 |

### Theoretical Framework

The BAR/PBSA method is based on the theory that the protein dielectric constant in Poisson-Boltzmann models correlates with the strength of electronic polarization in the protein environment:

| Concept | Explanation |
|---------|-------------|
| epsin=1 | Default for Amber force fields; no electronic polarization included; Coulombic interactions unscreened |
| epsin=2 | Commonly used to approximately include electronic polarization screening |
| epsin=4 | Higher screening; tested for strongly charged binding pockets |
| radiscale | Optimized to minimize absolute deviation between PBSA and explicit-solvent electrostatic free energies for ligand |
| protscale | Optimized to minimize absolute deviation between PBSA and explicit-solvent electrostatic free energies for complex |

Note: The additive force fields are developed with effective partial charges that already include polarization responses to the environment (mostly in water) in an averaged, mean-field manner. They are not fully compatible with the theoretical dielectric constant of 2 because polarization is already partially accounted for in the effective partial charges. Thus, a scanning procedure to find the optimal solute dielectrics is necessary.

## Worked Example

### Example Run for a Single Ligand/Complex Pair

From the manual, an example run for processing a single ligand trajectory and a single complex trajectory can be found in `$AMBERHOME/AmberTools/test/bar_pbsa/`.

**Check the test directory:**
```bash
ls $AMBERHOME/AmberTools/test/bar_pbsa/
```

**Step 1: Strip source trajectories**
```bash
python bar_pbsa.py strip strip_input.yaml
```

**Step 2: Prepare sander input files**
```bash
python bar_pbsa.py prep prep_input.yaml
```

**Step 3: Run sander in parallel (ligand and complex separately)**
```bash
python bar_pbsa.py run lig_input.yaml -n 8
python bar_pbsa.py run com_input.yaml -n 8
```

**Step 4: Calculate final decharging energies**
```bash
python bar_pbsa.py calc lig_input.yaml
python bar_pbsa.py calc com_input.yaml
```

The log with the final energies will be saved in the directory of the run.

### Interpreting Results

The final decharging energies are calculated using the Bennett Acceptance Ratio (BAR) method at each lambda window, and the data is aggregated to obtain the decharging energy for the full process. The binding free energy can be obtained by subtracting the ligand decharging energy from the complex decharging energy:

```
DeltaG_binding = DeltaG_complex - DeltaG_ligand
```

The effect of electronic polarization is assessed by comparing binding free energies calculated at different `epsin` values. If the binding free energy changes significantly with `epsin`, it indicates that electronic polarization effects are important for the system.

## Key Takeaways

1. **BAR/PBSA post-processes alchemical simulation trajectories** to incorporate electronic polarization effects in a continuum manner. It addresses the limitation of standard additive force fields that cannot handle polarization effects upon ligand transfer.

2. **The four-stage pipeline** (strip, prep, run, calc) automates the entire workflow: stripping solvent/ions, preparing PBSA inputs, running parallel sander calculations, and calculating final decharging energies via BAR.

3. **Radii calibration is essential** before investigating electronic polarization. First optimize `radiscale` for the ligand, then `protscale` for the complex, minimizing absolute deviation between PBSA and explicit-solvent electrostatic free energies.

4. **Electronic polarization is probed by varying `epsin`**. The default `epsin=1` is consistent with additive force field design but includes no polarization screening. Higher values (2, 4, etc.) approximately incorporate electronic polarization, screening electrostatic interactions.

5. **The method requires pre-computed alchemical decharging trajectories** organized in a specific folder structure (lambda windows as subdirectories with `ti.parm7` and `ti001.nc` files). Multiple replicates (t1, t2, ...) improve convergence.

## Connects To

- [ch06](ch06-parmed-topology.md) -- ParmEd topology manipulation (radii assignment)
- [ch07](ch07-md-engines-input.md) -- MD engines and mdin reference (sander PBSA)
- [ch11](ch11-free-energy-ti.md) -- Free energy and TI calculations (alchemical simulations)
- [ch12](ch12-mmpbsa.md) -- MMPBSA.py (related PB-based free energy method)
- [ch20](ch20-implicit-solvent.md) -- Implicit solvent models (PBSA details)