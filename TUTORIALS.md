# Tutorial Cross-Reference

This file maps each skill chapter to the corresponding tutorial files in `manuals/tutorials/`. When you query a topic, use this to find both the distilled chapter AND the original worked example.

## Core Workflow

| Chapter | Tutorials |
|---------|-----------|
| ch01 - Installation & Quick Start | `01_tutorials.md`, `02_overview.md`, `04_case_study.md`, `13_tools.md` |
| ch02 - Force Fields | `01-2_Fundamentals of LEaP.md`, `01-8_Un-natural amino acids ff15ipq-m.md`, `01-10_Using the Pantetheine Force Field Library.md`, `02_Developing Nonstandard Parameters.md` |
| ch03 - LEaP System Building | `01_building_systems.md`, `01-1_Preparing Structure.md`, `01-2_Fundamentals of LEaP.md`, `01-3_Building Systems with CHARMM-GUI.md`, `01-5_Building a Peptide Sequence.md`, `01-6_Building Protein Systems in Explicit Water.md`, `01-6-1_Calculating Salt Molarity.md`, `01-7_Simulation of a protein crystal.md`, `01-11_Building and Simulating an Ionic Liquid.md`, `01-12_Material Systems.md`, `01-12-1_Protein-Metal.md`, `01-12-2_PET polymer.md`, `01-12-3_Hydroxyapatite-Water.md`, `03_tutorial0.md`, `05-1_Simple Simulation of Alanine Dipeptide.md` |
| ch04 - PDB Preparation | `01-1_Preparing Structure.md` |
| ch05 - Antechamber & GAFF | `02-1_Simulating a pharmaceutical compound with Antechamber and GAFF.md`, `02-1-1_Original antechamber tutorials.md`, `02-2_Setting up a DNA-Ligand System.md` |
| ch06 - parmed Topology | `01-4_Hydrogen Mass Repartitioning.md` |
| ch07 - sander Namelist Reference | `03_Creating Stable Systems and Running MD.md`, `03-3_Running MD with pmemd.md` |
| ch08 - Minimization & Relaxation | `03_Creating Stable Systems and Running MD.md`, `03-1_Relaxation of Explicit Water Systems.md`, `03-2_Relaxation of Implicit Solvent System GB.md`, `03_tutorial0.md`, `05-1_Simple Simulation of Alanine Dipeptide.md` |
| ch09 - pmemd & Production MD | `03-3_Running MD with pmemd.md`, `03-4_Running MD in Parallel.md`, `05-1_Simple Simulation of Alanine Dipeptide.md` |

## Trajectory Analysis

| Chapter | Tutorials |
|---------|-----------|
| ch10 - CPPTRAJ | `04-1_An Introduction to CPPTRAJ.md`, `04-2_RMSD Analysis in CPPTRAJ.md`, `04-3_PCA with CPPTRAJ.md`, `04-4_Combined Clustering Analysis with CPPTRAJ.md`, `04-5_AMBER-Hub CPPTRAJ website.md`, `04-6_Analysis of MAD2 PCA tICA Markov State Models.md`, `04-7_Analysis of T-REMD Simulations.md` |

## Free Energies

| Chapter | Tutorials |
|---------|-----------|
| ch11 - Thermodynamic Integration | `07_Free_Energies.md`, `07-1_TI using soft core potentials.md`, `07-2_Advanced TI Using ACES.md`, `07-3_pKa Calculations using TI.md` |
| ch12 - MM-PBSA / MMPBSA.py | `07-4_MM-PBSA.md` |
| ch13 - Umbrella Sampling & NFE | `07-5_Umbrella sampling alanine dipeptide.md`, `07-6_Umbrella sampling methanol membrane.md`, `07-12_The Nonequilibrium Free Energy (NFE) Toolkit for pmemd.md`, `06-1_Adaptive Steered Molecular Dynamics.md` |
| ch28 - BAR/PBSA Post-processing | `07_Free_Energies.md`, `07-4_MM-PBSA.md` |
| ch35 - Free Energies Complete Methodology | `07_Free_Energies.md`, `07-1_TI using soft core potentials.md`, `07-2_Advanced TI Using ACES.md`, `07-3_pKa Calculations using TI.md`, `07-5_Umbrella sampling alanine dipeptide.md`, `07-6_Umbrella sampling methanol membrane.md`, `07-7_Free energy estimation EMIL.md`, `07-8_Computing binding enthalpy values.md`, `07-9_FEW Free Energy Workflow Tool.md`, `07-10_GIST Factor Xa.md`, `07-10-1_GIST streptavidin-biotin.md`, `07-11_APR Method.md`, `07-12_The Nonequilibrium Free Energy (NFE) Toolkit for pmemd.md` |

## Enhanced Sampling & Equilibria

| Chapter | Tutorials |
|---------|-----------|
| ch14 - Enhanced Sampling | `06_Sampling_Configuration_Space.md`, `06-1_Adaptive Steered Molecular Dynamics.md`, `06-2_The unified middle thermostat scheme.md`, `06-3_Weighted Ensemble Methods using WESTPA.md`, `06-4_Temperature Replica Exchange MD REMD.md` |
| ch15 - Constant pH & Redox | `08_Chemical_Reactions_and_Equilibria.md`, `08-1_Constant pH MD Lysozyme.md`, `08-2_Constant pH and Redox Potential MD.md`, `08-2-1_Marcus electron transfer parameters.md` |

## Force Field Development

| Chapter | Tutorials |
|---------|-----------|
| ch16 - Force Field Development (mdgx, RESP) | `02_Developing Nonstandard Parameters.md`, `02-3_Simulating GFP and building a modified amino acid.md`, `02-3-1_Recent modified-residue tutorial Carlos Ramos.md`, `02-5_Electrostatic Parameterization with Pyresp.py.md`, `02-6_RESP Charge Derivation with PyPE_RESP.md`, `02-6-1_From prior ESP calculations Gaussian GAMESS.md`, `02-6-2_Fully automated with QUICK.md`, `02-7_Deriving Implicitly Polarized Charges with mdgx.md`, `02-8_Deriving Custom Force Field Parameters with mdgx.md`, `02-9_Adding Custom Extra Points to a Model.md`, `05-2_Simulating GFP and building a modified amino acid residue.md`, `05-3_Using mdgx to manipulate small molecules.md` |
| ch17 - Metal Ion Modeling | `02-4_Metal Ion Modeling Tutorial.md`, `01-12-1_Protein-Metal.md` |
| ch24 - paramfit | `02-8_Deriving Custom Force Field Parameters with mdgx.md` |

## System Types

| Chapter | Tutorials |
|---------|-----------|
| ch18 - Membrane Systems | `01-9_Building Membrane Systems-Overview.md`, `01-9-1_Using PACKMOL-Memgen.md`, `01-9-2_Lipi21 and PACKMOL-Memgen Setup.md`, `01-9-3_Lipid-building Desmond Maestro PyMOL.md` |
| ch19 - QM/MM Overview | `08_Chemical_Reactions_and_Equilibria.md` |
| ch20 - Implicit Solvent Overview | `03-2_Relaxation of Implicit Solvent System GB.md`, `05-5_Simulation of DNA with Implicit Solvent and Explicit Ions GBION.md`, `01-13_Using 3D-RISM and MOFT.md` |
| ch25 - ProPrep | `01-1_Preparing Structure.md` |
| ch26 - LES | (no matching tutorial) |
| ch27 - NAB Nucleic Acid Builder | `05-2-1_Simulating the MTR1 ribozyme.md`, `05_Case Studies.md` |
| ch34 - QM/MM Detailed Reference | `08-3_Quantum dynamical effects liquid water.md`, `05-2_Simulating GFP and building a modified amino acid residue.md` |

## Advanced Tools & Standalone Programs

| Chapter | Tutorials |
|---------|-----------|
| ch21 - Advanced Tools (FEW, GIST, APR, EMIL) | `07-7_Free energy estimation EMIL.md`, `07-8_Computing binding enthalpy values.md`, `07-9_FEW Free Energy Workflow Tool.md`, `07-10_GIST Factor Xa.md`, `07-10-1_GIST streptavidin-biotin.md`, `07-11_APR Method.md` |
| ch22 - NMR, CryoEM, SAXS | (no matching tutorial) |
| ch23 - sqm Semi-empirical QM | `02-1-1_Original antechamber tutorials.md` |
| ch29 - RISM Detailed | `01-13_Using 3D-RISM and MOFT.md`, `05-4_Calculating ion distributions around DNA using 3D-RISM.md` |
| ch30 - Torch PBSA | (no matching tutorial) |
| ch31 - GBNSR6 | `03-2_Relaxation of Implicit Solvent System GB.md`, `05-5_Simulation of DNA with Implicit Solvent and Explicit Ions GBION.md` |
| ch32 - External Library Interface | (no matching tutorial) |
| ch33 - PBSA Detailed | `07-4_MM-PBSA.md`, `07-8_Computing binding enthalpy values.md` |

## Tutorials by Section

### 05_building_systems/ → ch02, ch03, ch04, ch06, ch17, ch18, ch20, ch25, ch29
### 06_Developing_Nonstandard_Parameters/ → ch05, ch16, ch17, ch23, ch24
### 07_Creating_Stable_Systems_and_Running_MD/ → ch07, ch08, ch09
### 08_Trajectory_Analysis/ → ch10
### 09_Case_Studies/ → ch03, ch08, ch09, ch10, ch16, ch20, ch27, ch29, ch34
### 10_Sampling_Configuration_Space/ → ch13, ch14
### 11_Free_Energies/ → ch11, ch12, ch13, ch21, ch28, ch33, ch35
### 12_Chemical_Reactions_and_Equilibria/ → ch15, ch19, ch34
