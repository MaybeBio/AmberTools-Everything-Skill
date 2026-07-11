# AmberTools Glossary

**3D-RISM** -- 3D Reference Interaction Site Model; implicit solvent model computing solvent density distribution on a 3D grid (Ch 20, 21)

**ACES** -- Alchemical Enhanced Sampling; advanced TI method for ligand binding free energy (Ch 11)

**AM1-BCC** -- Semi-empirical charge method; AM1 + bond charge correction (Ch 5)

**AMBERHOME** -- Environment variable pointing to Amber installation directory (Ch 1)

**antechamber** -- Small molecule parameterization tool; generates mol2/frcmod for GAFF (Ch 5)

**APR** -- Attach-Pull-Release method; absolute binding free energy via staged pulling (Ch 21)

**BAR** -- Bennett Acceptance Ratio; free energy estimator from energy differences (Ch 11, 21)

**BCC** -- Bond Charge Correction; AM1-BCC charge method (Ch 5)

**C(pH)MD** -- Continuous constant pH molecular dynamics (Ch 15)

**cphmd** -- Constant pH MD; titratable residue protonation states (Ch 15)

**cpptraj** -- CPPTRAJ; trajectory analysis program (Ch 10)

**cutoff** -- Nonbonded interaction cutoff; cut=8.0 for PME, cut=12.0+ for GB (Ch 7)

**CryoEM** -- Cryo-electron microscopy; EMAP restraints for map fitting (Ch 22)

**DFTB** -- Density Functional Tight Binding; SCC-DFTB semi-empirical QM (Ch 19)

**EMIL** -- Einstein Molecule; absolute free energy integration (Ch 21)

**Ewald** -- Ewald summation for long-range electrostatics (Ch 7)

**FEW** -- Free Energy Workflow tool; automated MM-PBSA/GBSA, LIE, TI setup (Ch 21)

**FF14SB** -- Amber protein force field (2014) (Ch 2)

**FF19SB** -- Amber protein force field (2019), current recommended (Ch 2)

**frcmod** -- Force field modification file; additional parameters (Ch 2, 5, 16)

**FRET** -- Forster Resonance Energy Transfer; FRETrest restraints (Ch 22)

**GAFF** -- General Amber Force Field; small organic molecules (Ch 5)

**GAFF2** -- Updated GAFF with improved torsion parameters (Ch 5)

**GB** -- Generalized Born implicit solvent model (Ch 20)

**GBION** -- GB with explicit ions; implicit solvent + explicit counterions (Ch 20)

**GBNSR6** -- GB model using R6 integral for effective radii (Ch 20)

**GIST** -- Grid Inhomogeneous Solvation Theory; water thermodynamics mapping (Ch 21)

**HMR** -- Hydrogen Mass Repartitioning; enables 4 fs timestep (Ch 6)

**inpcrd** -- Amber input coordinate file (Ch 3, 7)

**LEaP** -- System building program; tleap (text) and xleap (GUI) (Ch 3)

**leaprc** -- LEaP resource file; loads force field parameters (Ch 2, 3)

**lib** -- Amber library file; residue template (Ch 3)

**libsff** -- Simple Force Field library; used by NAB (Ch 22)

**Lipid21** -- Amber lipid force field (Ch 2, 18)

**MCPB** -- Metal Center Parameter Builder (Ch 17)

**mdcrd** -- Amber trajectory file (ASCII format) (Ch 7)

**mdgx** -- Force field parameter development tool (Ch 16)

**mdin** -- MD input file with namelists (Ch 7)

**MBAR** -- Multistate Bennett Acceptance Ratio (Ch 11, 21)

**middle thermostat** -- Unified scheme for configurational sampling (Ch 14)

**MM-PBSA** -- Molecular Mechanics Poisson-Boltzmann Surface Area; binding free energy (Ch 12)

**MMPBSA.py** -- Python script for MM-PBSA/GBSA analysis (Ch 12)

**MoFT** -- Analysis of volumetric data; density/Laplacian mapping (Ch 21)

**mol2** -- Tripos Mol2 file format; small molecule structure (Ch 5)

**NAB** -- Nucleic Acid Builder; molecular scripting language (Ch 22)

**nabc** -- NAB compiler; NAB to C translator (Ch 22)

**namelist** -- Fortran namelist input format (`&cntrl`, `&ewald`, etc.) (Ch 7)

**NetCDF** -- Binary trajectory format; ioutfm=1 (Ch 7, 10)

**NFE** -- Nonequilibrium Free Energy toolkit (Ch 13)

**NMR** -- Nuclear Magnetic Resonance; restraint-based refinement (Ch 22)

**NNP** -- Neural Network Potential; ML-based force field (Ch 19, 22)

**NPT** -- Constant Number, Pressure, Temperature ensemble (Ch 7)

**NVE** -- Constant Number, Volume, Energy ensemble (Ch 7)

**NVT** -- Constant Number, Volume, Temperature ensemble (Ch 7)

**off** -- Amber OFF library file format (Ch 3)

**PACKMOL-Memgen** -- Membrane system builder (Ch 18)

**parm7** -- Amber topology file; same as prmtop (Ch 3)

**parmed** -- ParmEd; topology manipulation tool (Ch 6)

**PBC** -- Periodic Boundary Conditions (Ch 7)

**PBSA** -- Poisson-Boltzmann Surface Area; continuum electrostatics (Ch 20)

**PCA** -- Principal Component Analysis; dimensionality reduction (Ch 10)

**pdb4amber** -- PDB cleanup and preparation tool (Ch 4)

**PME** -- Particle Mesh Ewald; efficient long-range electrostatics (Ch 7)

**pmemd** -- Particle Mesh Ewald MD engine (CPU) (Ch 7, 9)

**pmemd.cuda** -- GPU-accelerated pmemd; recommended for production (Ch 7, 9)

**prep** -- Amber PREP input file; residue definition (Ch 3)

**prepi** -- Amber PREP input file format variant (Ch 3)

**prmtop** -- Amber topology file; parameters, atom types, bonds (Ch 3)

**py_resp.py** -- Python RESP charge fitting tool (Ch 16)

**pyMSMT** -- Python Metal Site Modeling Toolbox (Ch 17)

**QM/MM** -- Quantum Mechanics / Molecular Mechanics hybrid (Ch 19)

**QUICK** -- Ab initio quantum chemistry program (Ch 19)

**REMD** -- Replica Exchange MD; enhanced sampling (Ch 14)

**RESP** -- Restrained Electrostatic Potential fit (Ch 5, 16)

**restraint** -- Positional or internal coordinate restraint (Ch 8)

**RISM** -- Reference Interaction Site Model (Ch 20, 21)

**RMSD** -- Root Mean Square Deviation; structural comparison (Ch 10)

**RMSF** -- Root Mean Square Fluctuation; per-residue flexibility (Ch 10)

**rst7** -- Amber restart file; coordinates + velocities (Ch 7)

**sander** -- Simulated Annealing with NMR-Derived Energy Restraints; core MD engine (Ch 7)

**SAXS** -- Small-Angle X-ray Scattering; solution structure (Ch 22)

**SCC-DFTB** -- Self-Consistent Charge DFTB (Ch 19)

**sff** -- Simple Force Field; NAB molecular mechanics library (Ch 22)

**SHAKE** -- Bond constraint algorithm; ntc=2, ntf=2 (Ch 7)

**softcore** -- Modified LJ/Coulomb potentials for alchemical TI (Ch 11)

**solvateoct** -- LEaP command; truncated octahedron water box (Ch 3)

**sqm** -- Semi-empirical quantum chemistry program (Ch 5, 19)

**steered MD** -- Non-equilibrium pulling; SMD and Jarzynski (Ch 14)

**TI** -- Thermodynamic Integration; alchemical free energy (Ch 11)

**tleap** -- Text-based LEaP; system building (Ch 3)

**topology** -- prmtop file; system connectivity, atom types, parameters (Ch 3)

**torch** -- LibTorch; ML runtime for PBSA (Ch 20)

**trajectory** -- Time series of coordinates; mdcrd or NetCDF (Ch 7, 10)

**umbrella sampling** -- Enhanced sampling along reaction coordinate (Ch 13)

**WESTPA** -- Weighted Ensemble Simulation Toolkit (Ch 14)

**WHAM** -- Weighted Histogram Analysis Method; umbrella sampling analysis (Ch 13)

**xleap** -- GUI version of LEaP (Ch 3)