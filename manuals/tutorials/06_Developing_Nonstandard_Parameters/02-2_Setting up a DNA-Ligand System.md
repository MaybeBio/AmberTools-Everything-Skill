

# Setting up a DNA-Ligand System
 by Bryan Lelanda, David Paula, Brent Kruegera and Ross Walkerb 

aDept. of Chemistry, Hope College* bSan Diego Supercomputer Center, University of California, San Diego
## Table of Contents
 - [Derive charges for the dye and linker.](https://ambermd.org/tutorials/advanced/tutorial1/section1.php)
- [Create a single PDB of the dye and the linker for use in step 3.](https://ambermd.org/tutorials/advanced/tutorial1/section2.php)
- [Build dye and linker Leap library file and frcmod file containing custom parameters.](https://ambermd.org/tutorials/advanced/tutorial1/section3.php)
- [Prepare 5' end of polyAT for attachment to dye and linker, attach dye and linker to polyAT and save prmtop and inpcrd files.](https://ambermd.org/tutorials/advanced/tutorial1/section4.php)
 ## Introduction
 In this tutorial we will work through an advanced scenario for preparing a system, for simulation with sander, that contains several non-standard residues where the use of Antechamber and the GAFF force field may not be appropriate (see [http://www.ambermd.org/tutorials/basic/tutorial4b/index.php](http://www.ambermd.org/tutorials/basic/tutorial4/index.php) for an example of using Antechamber). To illustrate this process we will choose a real world research problem rather than an 'artificially' simple example. This of course comes with several caveats. While it provides a more realistic tutorial it does of course present issues inherent in any research project. That is that there are often multiple ways to accomplish something and one does not necessarily know in advance the best option to choose. This is true here in the case of deriving certain RESP charges as you will see later on. However, this approach should provide you with a solid foundation for preparing your own systems. You just need to be aware that when things get complex there is no such thing as a step by step guide and therefore you should use this tutorial merely as a basis for approaching a similar problem. However, substantial thought and literature searching will likely be required before you tackle your own system.
 This tutorial centers on creating a prmtop and inpcrd file for a fluorescein dye bound to a poly-AT DNA decamer. The 3D structure, along with 2D schematics of the topology are shown below.
 ![](https://ambermd.org/tutorials/advanced/tutorial1/images/dna_linker_pic1.jpg)


 Figure 1. 3D structure ![](https://ambermd.org/tutorials/advanced/tutorial1/images/flo-and-linker.gif)


 Figure 2. 2D schematic
 The DNA is shown in blue and the fluorescein dye and its attached linker are shown in red. This sort of setup is typically used in laser probe experiments. In this tutorial we will go through a number of steps in order to ultimately create a prmtop and inpcrd file.
 *Note: This tutorial assumes that you have access to certain programs that are not part of the Amber suite. This includes Gaussian ( [http://www.gaussian.com](http://www.gaussian.com)) and Sirius ( [http://www.ngbw.org/sirius/](http://www.ngbw.org/sirius/)). Sirius is free of charge and you can download and install it but Gaussian is commercial. If you don't have Gaussian you can still complete this tutorial since all input and output files are provided. It is also possible to use the freely available GAMESS or Orca packages to calculate the ESPs for the RESP fit instead of Gaussian. However, this is not covered in this tutorial.*
 Last updated 2008.

