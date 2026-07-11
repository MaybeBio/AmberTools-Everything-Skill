---
name: AmberTools-Everything-Skill
description: "Knowledge base from Amber 2026 Reference Manual and official tutorials. Code-first technical reference for AmberTools26/Amber26 — MD simulation, force fields, free energy, enhanced sampling, and trajectory analysis. Use when writing Amber commands, input files, or MD workflows."
--- 

<!-- argument-hint: [program name, topic, or chapter number] -->

# Amber 2026 Technical Reference
**Source**: Amber 2026 Reference Manual (1112 pages, full distillation) + 83 official tutorials | **Generated**: 2026-07-11 | **Chapters**: 35

## How to Use This Skill

- **Without arguments** — load this core reference
- **With a program** — ask about `tleap`, `cpptraj`, `pmemd`, `antechamber`, `parmed`, `mdgx`, `MMPBSA.py`, `FEW`, etc.
- **With a workflow** — ask about `building systems`, `free energy`, `enhanced sampling`, `constant pH`, `membrane setup`
- **With chapter** — ask for `ch03` (LEaP), `ch07` (MD engines), `ch10` (CPPTRAJ), etc.
- **Browse** — ask "what chapters do you have?" to see the full index

When you ask about a topic not covered in Core Frameworks below, I will read the relevant chapter file before answering.

---

## Core Frameworks

### The Amber MD Pipeline

```
PDB → pdb4amber → LEaP (prmtop+inpcrd) → sander/pmemd (mdin) → cpptraj (analysis)
```

**Key files**: `prmtop` (topology/parameters), `inpcrd`/`rst7` (coordinates), `mdin` (namelist input), `mdcrd`/`nc` (trajectory)

### Force Fields (Section 3)

| Force Field | Use Case | `leaprc` source |
|---|---|---|
| FF19SB | Proteins | `source leaprc.protein.ff19SB` |
| FF14SB | Proteins (older) | `source leaprc.protein.ff14SB` |
| OL3/OL15/OL21 | RNA | `source leaprc.RNA.OL3` |
| OL15/OL21 | DNA | `source leaprc.DNA.OL15` |
| GAFF2 | Small organic molecules | Used via antechamber |
| Lipid21 | Lipids | `source leaprc.lipid21` |
| GLYCAM_06j | Carbohydrates | `source leaprc.GLYCAM_06j` |
| ff15ipq-m | Unnatural amino acids | `source leaprc.ff15ipq` |

Water models: OPC (recommended for FF19SB), TIP3P, TIP4P-Ew, OPC3, SPC/E

### MD Engines

| Engine | Hardware | Key capability |
|---|---|---|
| `sander` | CPU | Full feature set, QM/MM, NMR |
| `sander.MPI` | CPU multi-node | Parallel sander |
| `pmemd` | CPU | Optimized, subset of sander |
| `pmemd.cuda` | GPU | High performance, recommended |
| `pmemd.cuda.MPI` | Multi-GPU | Multi-node GPU |

**Run command**: `pmemd.cuda -O -i mdin -o mdout -p prmtop -c inpcrd -r rst7 -x mdcrd -v mdvel -e mden -inf mdinfo -ref refcrd`

### mdin Namelist Structure

```
&cntrl        ! General MD control
  imin=0,     ! 0=MD, 1=minimization
  ntx=1,      ! 1=restart from inpcrd, 5=restart from rst7
  nstlim=100, ! Number of MD steps
  dt=0.002,   ! Timestep in ps (2 fs)
  ntt=3,      ! 1=Berendsen, 2=Andersen, 3=Langevin
  temp0=300.0, ! Target temperature (K)
  ntp=1,      ! 0=no pressure scaling, 1=isotropic
  ntb=2,      ! 1=no PBC, 2=constant P
  ntc=2,      ! 1=no SHAKE, 2=SHAKE on H-bonds
  ntf=2,      ! 1=all bonds, 2=omit H-bond forces
  cut=8.0,    ! Nonbonded cutoff (Å)
  ioutfm=1,   ! 0=ASCII, 1=NetCDF trajectory
  ntwx=5000,  ! Trajectory write frequency
  ntpr=500,   ! Print to mdout frequency
  ntwr=5000,  ! Restart write frequency
/
&ewald        ! Ewald/PME settings
  skinnb=2.0, ! Skin width for pairlist
/
```

### CPPTRAJ — Core Analysis Commands

```
cpptraj -p prmtop -i script.in
```

**Data loading**: `parm`, `trajin` (load topology/trajectory)
**Actions**: `rms`, `rmsd`, `rmsf`, `distance`, `angle`, `dihedral`, `hbond`, `clustering`, `pca`, `strip`, `autoimage`, `center`
**Output**: `writedata`, `run`, `go`

**Atom mask syntax**: `:1-10@CA` (residues 1-10, Cα atoms), `@/H` (strip hydrogens), `:WAT` (water), `:Na+` (sodium ions)

### Antechamber — Small Molecule Parameterization

```
antechamber -i lig.mol2 -fi mol2 -o lig.mol2 -fo mol2 -c bcc -s 2 -nc 0
parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod -s 2
```

Charge methods: `-c bcc` (AM1-BCC), `-c resp` (RESP), `-c gas` (Gasteiger)

### LEaP — System Building

```
tleap -f leaprc.protein.ff19SB
> source leaprc.water.opc
> mol = loadpdb protein.pdb
> solvateoct mol OPCBOX 10.0
> addions2 mol Na+ 0
> addions2 mol Cl- 0
> saveamberparm mol prmtop inpcrd
> quit
```

---

## Chapter Index

### Core MD Workflow (ch01–ch09)
| # | Title | Key Programs |
|---|-------|-------------|
| [ch01](chapters/ch01-installation-quickstart.md) | Installation & Quick Start | configure, cmake, conda, AMBERHOME |
| [ch02](chapters/ch02-force-fields.md) | Force Fields & Molecular Mechanics | LEaP, leaprc, frcmod, ff19SB, ff14SB, OL3/15/21, Lipid21, GAFF2, OPC |
| [ch03](chapters/ch03-leap-system-building.md) | LEaP: System Building | tleap, xleap, solvateoct, addions2, loadpdb, saveamberparm |
| [ch04](chapters/ch04-pdb-preparation.md) | PDB Preparation | pdb4amber, reduce |
| [ch05](chapters/ch05-antechamber-gaff.md) | Antechamber & GAFF | antechamber, parmchk2, sqm, bcc, resp |
| [ch06](chapters/ch06-parmed-topology.md) | parmed: Topology Manipulation | parmed, HMR, frcmod |
| [ch07](chapters/ch07-sander-reference.md) | sander: Complete Namelist Reference | sander, sander.MPI, full &cntrl + &ewald + &wt |
| [ch08](chapters/ch08-minimization-relaxation.md) | Minimization, Heating & Relaxation | sander, pmemd, restraints, ntr, belly |
| [ch09](chapters/ch09-production-md.md) | pmemd & Production MD | pmemd, pmemd.cuda, pmemd.cuda.MPI, GPU, HIP, NFE, NEB |

### Trajectory Analysis (ch10)
| [ch10](chapters/ch10-cpptraj-analysis.md) | CPPTRAJ: Trajectory Analysis | cpptraj, RMSD, PCA, clustering, tICA, MSM |

### Free Energies (ch11–ch13, ch28, ch35)
| # | Title | Key Programs |
|---|-------|-------------|
| [ch11](chapters/ch11-free-energy-ti.md) | TI: Thermodynamic Integration | sander, pmemd, softcore, ACES |
| [ch12](chapters/ch12-mmpbsa.md) | MM-PBSA & MMPBSA.py | MMPBSA.py, MM-GBSA, GB, PB |
| [ch13](chapters/ch13-umbrella-sampling-nfe.md) | Umbrella Sampling & NFE | umbrella, WHAM, NFE toolkit, SMD, ABMD |
| [ch28](chapters/ch28-bar-pbsa.md) | BAR/PBSA Post-processing | bar_pbsa.py, edgembar, decharging |
| [ch35](chapters/ch35-free-energy-detailed.md) | Free Energies: Complete Methodology | TI theory, softcore, GaMD-TI, quadrature |

### Enhanced Sampling & Equilibria (ch14–ch15)
| [ch14](chapters/ch14-enhanced-sampling.md) | Enhanced Sampling | REMD, GaMD, aMD, targeted MD, NEB, LMOD, WESTPA |
| [ch15](chapters/ch15-constant-ph-redox.md) | Constant pH & Redox Potential | cphmd, cein, C(pH)MD, Marcus ET |

### Force Field Development (ch16–ch17, ch24)
| [ch16](chapters/ch16-force-field-development.md) | Force Field Development | mdgx, py_resp.py, RESP, IPolQ |
| [ch17](chapters/ch17-metal-ion-modeling.md) | Metal Ion Modeling | MCPB.py, pyMSMT, 12-6-4 LJ |
| [ch24](chapters/ch24-paramfit.md) | paramfit: Parameter Fitting | paramfit, force constants |

### System Types (ch18–ch20, ch25–ch27, ch34)
| [ch18](chapters/ch18-membrane-systems.md) | Membrane Systems | PACKMOL-Memgen, Lipid21, LiPi21 |
| [ch19](chapters/ch19-qmmm.md) | QM/MM Overview | sander.MPI, QM/MM, QUICK |
| [ch20](chapters/ch20-implicit-solvent.md) | Implicit Solvent Overview | igb, pbsa, 3D-RISM, GBNSR6 |
| [ch25](chapters/ch25-proprep.md) | ProPrep: Protein Preparation | propred, interface modes |
| [ch26](chapters/ch26-les.md) | LES: Locally-Enhanced Sampling | les, multiple copy simulation |
| [ch27](chapters/ch27-nab.md) | NAB: Nucleic Acid Builder | nab, nabc, nab2c, sff, molecule building |
| [ch34](chapters/ch34-qmmm-detailed.md) | QM/MM: Detailed Namelist Reference | &qmmm, DFTB3, GFN2-xTB, QUICK, link atoms |

### Advanced Tools (ch21–ch22)
| [ch21](chapters/ch21-advanced-tools.md) | Advanced Tools | FEW, GIST, APR, EMIL, MoFT |
| [ch22](chapters/ch22-nmr-cryoem-nab.md) | NMR, CryoEM, SAXS Refinement | sander, cpptraj |

### Standalone Programs (ch23, ch29–ch33)
| [ch23](chapters/ch23-sqm.md) | sqm: Semi-empirical QM | sqm, PM3, AM1, MNDO, DFTB, DFTB3, PM6, PM7, GFN2-xTB |
| [ch29](chapters/ch29-rism-detailed.md) | RISM: Detailed Reference | rism1d, rism3d, 3D-RISM, rism3d.snglpnt, XRISM, DRISM |
| [ch30](chapters/ch30-torch-pbsa.md) | Torch PBSA | LibTorch, GPU-PB |
| [ch31](chapters/ch31-gbnsr6.md) | GBNSR6: GB with R6 | gbnsr6, GB equations, R6 integration |
| [ch32](chapters/ch32-external-library.md) | External Library Interface | sander API, external energy/forces |
| [ch33](chapters/ch33-pbsa-detailed.md) | PBSA: Detailed Reference | pbsa, PB solver, grid, nonpolar, membrane, GPU-PBSA |

## Topic Index

- **addions2** → ch03
- **antechamber** → ch05
- **APR (Attach-Pull-Release)** → ch21
- **atom masks** → ch10
- **BAR** → ch28
- **bar_pbsa.py** → ch28
- **BCC charges** → ch05
- **CHARMM-GUI** → ch03
- **constant pH** → ch15
- **constant redox potential** → ch15
- **CPPTRAJ** → ch10
- **crystal simulation** → ch03
- **cutoff** → ch07
- **DFTB/DFTB3** → ch23
- **DNA** → ch02, ch03, ch27
- **dvdl** → ch35
- **EMIL** → ch21
- **Enhanced sampling** → ch14
- **external library** → ch32
- **FEW (Free Energy Workflow)** → ch21
- **FF14SB** → ch02
- **FF19SB** → ch02
- **force field** → ch02
- **frcmod** → ch02, ch16
- **GAFF/GAFF2** → ch05
- **GaMD** → ch14
- **GB (Generalized Born)** → ch20
- **GBNSR6** → ch31
- **GFN2-xTB** → ch23, ch34
- **GIST** → ch21
- **GPU** → ch09
- **HMR (Hydrogen Mass Repartitioning)** → ch06
- **igb** → ch20
- **implicit solvent** → ch20
- **LEaP** → ch03
- **LES** → ch26
- **lipid** → ch18
- **LMOD** → ch14
- **MCPB.py** → ch17
- **mdgx** → ch16
- **mdin** → ch07
- **membrane** → ch18
- **minimization** → ch08
- **MM-PBSA** → ch12
- **MoFT** → ch21
- **NAB** → ch27
- **nabc** → ch27
- **namelist** → ch07
- **NEB** → ch14
- **NetCDF** → ch07, ch10
- **NFE** → ch13
- **NMR** → ch22
- **OPC water** → ch02
- **PACKMOL-Memgen** → ch18
- **paramfit** → ch24
- **parmchk2** → ch05
- **parmed** → ch06
- **PBSA** → ch33
- **PCA** → ch10
- **pdb4amber** → ch04
- **pmemd** → ch09
- **pmemd.cuda** → ch09
- **pmemd.cuda.MPI** → ch09
- **prmtop** → ch03, ch06
- **ProPrep** → ch25
- **py_resp.py** → ch16
- **pyMSMT** → ch17
- **QM/MM** → ch19, ch34
- **QUICK** → ch19
- **REMD** → ch14
- **RESP** → ch05, ch16
- **restraints** → ch08
- **RISM** → ch29
- **RMSD** → ch10
- **RNA** → ch02, ch03, ch27
- **sander** → ch07
- **SHAKE** → ch07
- **softcore** → ch11, ch35
- **solvateoct** → ch03
- **sqm** → ch23
- **steered MD** → ch14
- **Thermodynamic Integration** → ch11, ch35
- **tleap** → ch03
- **Torch PBSA** → ch30
- **umbrella sampling** → ch13
- **WESTPA** → ch14
- **WHAM** → ch13
- **XMIN** → ch07

## Tutorial Traceability

When answering a query:
1. Find the relevant chapter(s) in the index below
2. Cross-reference with [TUTORIALS.md](TUTORIALS.md) to find the original worked examples in `manuals/tutorials/`
3. Load the tutorial file(s) for concrete, step-by-step code examples to accompany the chapter reference

Example query flow: "metal ion simulation" → ch17 (Metal Ion Modeling) → TUTORIALS.md shows: `02-4_Metal Ion Modeling Tutorial.md` + `01-12-1_Protein-Metal.md` → load those tutorials for worked examples.

## Supporting Files

- [glossary.md](glossary.md) — all key terms with definitions
- [patterns.md](patterns.md) — common workflows and design patterns
- [cheatsheet.md](cheatsheet.md) — quick reference tables and decision guides
- [TUTORIALS.md](TUTORIALS.md) — chapter → tutorial cross-reference for source traceability

---

## Scope & Limits

This skill covers AmberTools26 and Amber26 content from the official Reference Manual and tutorials. It is code-first: every concept is grounded in working commands, input file snippets, and exact flag syntax. For theory and background, consult the Amber Reference Manual PDF directly. For project-specific implementation, combine with the user's actual system setup.