# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Purpose

This repo is a staging area for generating an **AmberTools technical skill** — a code-centric, API-centric knowledge base distilled from the official Amber 2026 Reference Manual and the official tutorials. The target output is a skill (`.md` file consumable by Claude Code) that enables writing correct AmberTools commands, input files, and analysis scripts. Theory is included only when it clarifies API usage.

## Source Materials

- **`manuals/Amber26.pdf`** — Amber 2026 Reference Manual (1112 pages). The authoritative reference for all Amber programs, input file syntax, force fields, and API details.
- **`manuals/tutorials/`** — 84 files (markdown + 1 ipynb) scraped from https://ambermd.org/tutorials/. These are the official step-by-step walkthroughs covering the full MD workflow.

## Tutorials Organization

The tutorials follow the standard MD pipeline:

| Section | Directory | Content |
|---|---|---|
| 1 | `05_building_systems/` | System building: LEaP, PDB prep, membranes, explicit solvent, ionic liquids, materials |
| 2 | `06_Developing_Nonstandard_Parameters/` | Force field parameterization: antechamber, GAFF, RESP, mdgx, metal ions, custom residues |
| 3 | `07_Creating_Stable_Systems_and_Running_MD/` | Minimization, heating, relaxation (explicit/implicit), production MD with pmemd |
| 4 | `08_Trajectory_Analysis/` | CPPTRAJ: RMSD, PCA, clustering, Markov state models, T-REMD analysis |
| 5 | `09_Case_Studies/` | End-to-end examples: alanine dipeptide, GFP, DNA, 3D-RISM, GBION |
| 6 | `10_Sampling_Configuration_Space/` | Enhanced sampling: steered MD, middle thermostat, WESTPA, REMD |
| 7 | `11_Free_Energies/` | Free energy: TI, MM-PBSA, umbrella sampling, FEW, APR, NFE, GIST, EMIL |
| 8 | `12_Chemical_Reactions_and_Equilibria/` | Constant pH MD, redox potential, Marcus ET, quantum dynamics |

Index files: `01_tutorials.md`, `02_overview.md`, `03_tutorial0.md`, `04_case_study.md`, `13_tools.md`

## Key AmberTools Programs (for skill distillation)

The skill should cover these core programs and their input/output files:

- **tLEaP / xLEaP** — System building, force field loading, solvation, topology/coordinate generation
- **antechamber** — Small molecule parameterization for GAFF/GAFF2
- **parmed** — Topology manipulation, hydrogen mass repartitioning
- **sander / pmemd / pmemd.cuda** — MD engines (CPU and GPU)
- **cpptraj** — Trajectory analysis (the primary analysis tool, successor to ptraj)
- **mdgx** — Force field parameterization and conformational sampling
- **MCPB.py** — Metal center parameterization
- **pdb4amber** — PDB preparation for Amber
- **FEW** — Free Energy Workflow tool
- **3D-RISM** — Solvation thermodynamics

## File Format Conventions

Key Amber file types that must be accurately represented in the skill:
- `prmtop` / `parm7` — Topology/parameter file
- `inpcrd` / `rst7` — Coordinate/restart file
- `mdin` — MD input control file (namelist format)
- `mdcrd` / `nc` — Trajectory file (ASCII NetCDF)
- `frcmod` — Force field modification file
- `lib` / `prep` — Residue library/prep files
- `leaprc` — LEaP resource file

## Skill Generation Guidelines

1. **Code-first**: Every concept should be grounded in a working command, input file snippet, or script. If a feature is described, show the exact input flags and syntax.
2. **API precision**: Flag names, namelist variables, and command syntax must match the manual exactly. Do not paraphrase option names.
3. **Progressive disclosure**: Start with minimal working examples, then layer on advanced options.
4. **Cross-reference**: Tutorial examples should reference the relevant manual sections for exhaustive parameter documentation.
5. **Deprecation awareness**: Older tutorials may reference deprecated programs (ptraj, old force fields). The skill should prefer current equivalents (cpptraj, ff19SB/ff14SB).