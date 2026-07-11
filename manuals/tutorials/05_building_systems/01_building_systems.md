
# Building Systems
  ## Building Basic Systems
 The first part of a MD simulation is to build the system. The types of forcefields you will use depends on the type of structures you are using as part of your system, and what interactions you are interested in studying. The Amber Manual describes the types of forcefields in great detail. Input files, such as PDB files, will often need to be modified to ensure compatibility with Amber. The following tutorials include useful methods for setting up simulations and show the thought process behind these choices. Usually the more care that is used in setting a system will result in better data.
 [**Preparing Structure**](https://ambermd.org/tutorials/basic/tutorial9/index.php)


  This page highlights the existence of *pdb4amber*, the recommended method for preparing a pdb file for use in LEaP. This program was written by Romain Wolf and modified later mainly by Hai Nguyen. 



  [**Fundamentals of LEaP**](https://ambermd.org/tutorials/pengfei/index.php)


  This tutorial demonstrates rudimentary system building in LEaP, the main program for preparing simulations in AMBER. The tutorial serves as a nice reference for different kinds of files used in the LEaP program, introduces the workflow, and provides a simple example of building a protein in water. 



  [**CHARMM-GUI**](https://ambermd.org/tutorials/CHARMM-GUI.php)


 Some people find it easier to set up a system for simulation in Amber using the CHARMM-GUI interface. These tutorials were written by Wonpil Im at Lehigh University. He also has YouTube videos on these tutorials. 



  [**Hydrogen Mass Repartitioning**](https://ambermd.org/tutorials/basic/tutorial12/index.php)


 This tutorial shows how parmed can be used to perfom Hydrogen Mass Repartitioning, a technique that allows you to use a larger time step in your simulation. It is adapted from Adrian Roitberg's work.
  ## Building Specific Types of Systems
 The following tutorials demonstrate how to set up systems for particular types of systems that are common or interesting. They should serve as food for thought when creating your own system. 

 

  [**Building a Peptide Sequence**](https://ambermd.org/tutorials/basic/tutorial10/index.php)


  This tutorial utilizes the sequence command within LEaP to build an extended structure of the 10 amino acid chignolin. It also sets up a topology file and starting coordinates for an implicit solvent calculation with the GBneck2 parameters. 



  [**Building Protein Systems in Explicit Water**](https://ambermd.org/tutorials/basic/tutorial7/index.php)


  This tutorial will go over how to build an explicitly solvated protein system in LEaP. First, users will use basic features in VMD to examine a protein structure. Because preparing experimentally-determined protein structures for simulations are not straightforward, this tutorial will also help users evaluate and change the pdb file as necessary. Finally, users will use LEaP to build a protein system in explicit solvent including counterions and a buffer of ions. 



  [**Calculating Salt Molarity in an Explicit Water System**](https://ambermd.org/tutorials/basic/tutorial8/index.php)


  This tutorial show users how to calculate the number of buffer ions to add in LEaP based on an overall molarity and the volume of a truncated octahedral box of explicit water. 



  [**Simulation of a Protein Crystal**](https://ambermd.org/tutorials/advanced/tutorial13/XtalTutor1.php)


 This tutorial shows how to set up protein crystal structure, which requires different programs than other periodic systems. By David Cerutti and Dave Case. 

 

  [**Un-natural amino acids: the Amber ff15ipq-m forcefield**](https://ambermd.org/tutorials/advanced/tutorial36/index.php)


 This tutorial builds a tetrapeptide using the Amber forcefield ff15ipq-m, which is used to simulate peptides with unnatural amino acids in explicit solvent. 

 

  [**Building Membrane Systems - Overview**](https://ambermd.org/tutorials/MembraneSystems.php)


  There are a number of ways to build a membrane to be compatible with Amber including PACKMOL-Memgen, CHARMM-GUI and Maestro. 

 

  [**Using the Pantetheine Force Field Library**](https://ambermd.org/tutorials/basic/tutorial20/index.php)


 The Pantetheine Force Field (PFF) library is an Amber-compatible force field library for various pantetheine-containing ligands (PCLs). The PFF library was parameterized using the Gasteiger, AM1-BCC, or RESP charging methods in combination with the gaff2 and ff14SB parameter sets. In this tutorial, we are going to demonstrate how to use the parameter files in the PFF library for MD simulations using the AMBER molecular dynamics engines on two protein-PCL systems: Phosphopantetheine adenylyltransferase (PPAT) bound with phosphopantetheine (Ppant); and 3-hydroxy-3-methylglutaryl synthase/acyl carrier protein complex (HGMS/ACP) with phosphopantetheinyl serine (Ppant-Ser) nonstandard amino acid residue. Note that in the PPAT-Ppant complex, Ppant is a standalone ligand that binds to the protein non-covalently. While in the HGMS/ACP-Ppant-Ser complex, Ppant-Ser is part of the polymer chain of ACP. This difference will lead to slightly different parameterizations for the two types of complexes. 



  [**Building and Simulating an Ionic Liquid**](https://ambermd.org/tutorials/advanced/tutorial15/Tutorial2.php)


 This tutorial illustrates the use of antechamber and sander to carry out some simple simulations of a room-temperature (non-biological!) ionic liquid. By Chris Lim. 



  [**Material Systems**](https://ambermd.org/tutorials/advanced/tutorial27/index.php)


  Together with the INTERFACE force field, this tutorial explores the functions of employing the AMBER software package to modeling material and interfacial systems. By Pengfei Li. 



    [**Using 3D-RISM and MOFT to place waters and ions**](http://dansindhikara.com/Tutorials/Entries/2012/1/1_Using_3D-RISM_and_PLACEVENT.html)


  This tutorial shows how the 3D-RISM method can be used to generate an initial configuration of water around a solute molecule. By Daniel Sindhikara.

