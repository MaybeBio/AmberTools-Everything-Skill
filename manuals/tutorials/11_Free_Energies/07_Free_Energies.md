

# Free Energy Calculations in Amber
  [**Thermodynamic Integration using soft core potentials**](https://ambermd.org/tutorials/advanced/tutorial9/index.html)


 Thermodynamic Integration (TI) is a technique used to compute relative binding free enegy. This tutorial uses TI on two simple ligands, benzene and phenol, and their binding with the T4-lysozyme mutant L99A. 



  [**Advanced Thermodynamic Integration Using ACES**](https://ambermd.org/tutorials/advanced/tutorial39/index.php)


 This tutorial calculates the relative binding free energies between two inhibitors bound to the the human THr160-phospho CDK2/cyclin A protein. The methodology employed is a little more sophisticated that the standard TI. 



  [**pKa Calculations using Thermodynamic Integration**](https://ambermd.org/tutorials/advanced/tutorial6/index.php)


  This tutorial reproduces the calculation of the pKa value of the ASP residue in the protein thioredoxin as described in the following paper: Simonson, T., Carlsson, J., Case, D.A., *"Proton Binding to Proteins: pKa Calculations with Explicit and Implicit Solvent Models"*, **JACS** 2004, 126, pp 4167-4180. By Ross Walker and Mike Crowley. 



  [**Molecular Mechanics with a Poisson Boltzmann/Surface Area solvent:MM-PBSA**](https://ambermd.org/tutorials/advanced/tutorial3/index.php)


  This tutorial provides a step by step explanation of using the mm_pbsa script to calculate the binding energy of the RAS-RAF protein complex. It also includes instructions on using the mmpbsa_py script to perform these calculations as well. ( [Japanese translation](http://amber.tkanai-lab.org/)) By Ross Walker and Thomas Steinbrecher. 



  [**Umbrella Sampling: A Potential of Mean Force in alanine dipeptide Phi/Psi Rotation**](https://ambermd.org/tutorials/advanced/tutorial17/index.php)  In this tutorial we will learn how to use the AMBER software coupled with the Weighted Histogram Analysis Method (WHAM) of Alan Grossfield to generate potentials of mean force. Often one might want to know what the free energy profile is along a specific reaction coordinate. Such a profile is known as a potential of mean force and it can be very useful for identifying transition states, intermediates as well as the relative stabilities of the end points. At first thought one might think that you could generate a free energy along a specific reaction coordinate by just running an MD simulation and then looking at the probabilities of the states sampled. However, often the energy barrier of interest is many times the size of kbT and so the MD simulation will either remain in the local minimum it started in or cross to different minima but very very rarely sample the transition state. Umbrella sampling offers a way to effectively force the system to move through a transition state and reaction pathway that chemical knowledge of the system under study suggests is important. By Ross Walker and Thomas Steinbrecher. 



  [**Umbrella Sampling: Transferring methanol through a membrane**](https://github.com/callumjd/AMBER-Umbrella_COM_restraint_tutorial)  This tutorial illustrates how to use the *fxyz* capability to carry out umbrella sampling in specific directions. It uses the AMBER16 center-of-mass (COM) umbrella restraint code to determine the free energy of transfer profile for a methanol molecule through a DMPC membrane bilayer. The methanol molecule is first pulled from the center of the membrane out into the water phase. From the pulling step, we extract starting positions with methanol at 0, 2, 4, ..., 32 A from the membrane center. We run windows with methanol restrained at each of these positions. From the fluctuation in the z-position, we can construct the free energy profile using WHAM. Finally, we use the same information to derive the z-diffusion and z-resistance profiles and an estimate of the overall permeability coefficient. By Callum Dickson. 



  [**Free Energy Estimation for Conformers of Dialanine using EMIL**](https://ambermd.org/tutorials/advanced/tutorial19/index.php) [ *AMBER v14 and later*]


  This tutorial is a walk-through of absolute free energy calculations using EMIL. EMIL works by perturbing the "normal" atomistic representation of the system into a model for which the free energy is exactly known, thus it is a sort of thermodynamic integration tool. The example in this tutorial is an estimate of the four free energy basins of the alanine dipeptide.By Joshua T. Berryman. 



  [**Computing Binding Enthalpy Values**](https://ambermd.org/tutorials/advanced/tutorial21/index.php)


  In this tutorial, we will learn how to use AMBER to compute precise binding enthalpies using explicit water molecular dynamics simulations. Here we focus on the guest B2 binding to the host CB7. However, while the example in this tutorial focuses on a host-guest system, the technique is directly applicable to protein-ligand systems although we caution that the such systems will require considerably more sampling to converge than the host-guest system in this tutorial. By Andrew T. Fenley and Michael K. Gilson. 



  [**Free energy calculations with the Free Energy Workflow Tool (FEW)**](https://ambermd.org/tutorials/advanced/tutorial24/index.php)  This tutorial demonstrates the functionality of the Free Energy Workflow tool FEW. Using a sample data set of inhibitors of the protein Factor Xa it is shown how FEW can be used to easily prepare MD simulations and binding free energy calculations by the MM-PB(GB)SA, the linear interaction energy (LIE), and the thermodynamic integration (TI) method for multiple ligands binding to the same receptor. FEW provides an efficient way to setup and conduct binding free energy calculations with AMBER. By Nadine Homeyer and Holger Gohlke. 



  [**Grid Inhomogeneous Solvation Theory for water thermodynamics: in Factor Xa**](https://ambermd.org/tutorials/advanced/tutorial25/index.php)


  In this tutorial we will learn how to use the AMBER software coupled with the Grid Inhomogeneous Solvation Theory Method (GIST) to estimate thermodynamic values for the water molecules occupying the binding pocket of Factor Xa. Often one might want to estimate changes in hydration which are central to correctly describe biomolecular phenomena such as molecular recognition, drug binding, etc. However, it is difficult to estimate precisely the thermodynamic solvation contribution in these phenomena. By Romelia Salomon, Crystal Nguyen, Steven Ramsey, Jonathan D. Gough, and Ross Walker. 



  [**Computing Binding Free Energy using the Attach-Pull-Release (APR) Method**](https://ambermd.org/tutorials/advanced/tutorial29/index.php)


  This tutorial demonstrates how to use the attach-pull-release (APR) approach to compute the binding thermodynamics of a host-guest system with an explicit water model. Results from binding calculations using the APR approach showed moderate to strong correlations to the experimental measurements of cucurbit[7]uril (CB7), octa acid (OA), tetra-endo-methyl octa-acid (TEMOA), α- and β-cyclodextrin (α-CD and β-CD) with their guest molecules. By Jian Yin, Niel M. Henriksen, David R. Slochower, and Michael K. Gilson. 



  [**The Nonequilibrium Free Energy (NFE) Toolkit for PMEMD**](https://ambermd.org/tutorials/advanced/tutorial31/index.php)


  This tutorial introduces the Nonequilibrium Free Energy (NFE) toolkit, which is fully functional in SANDER and partially ported to PMEMD AMBER. The purpose of this tutorial is multifold. We review the current status of the porting of software to PMEMD AMBER v.16, and provide suitable patches for upgrading the toolkit to older versions of AMBER (v.14). Additionally, patches are provided for upgrading modules not released to PMEMD AMBER v.16, and we provide a step-by-step tutorial on how to write new collective variables for free energy calculations. Software templates for new collective variables are also given. Finally, a few examples of using the new PMEMD versions of the code is provided. By Feng Pan, Mahmoud Moradi, Christopher Roland, and Celeste Sagui.


