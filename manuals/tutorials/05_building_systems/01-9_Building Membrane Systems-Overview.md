

# Building Membrane Systems - Overview
 Building Membrane Systems is tricky in Amber. A [force field](https://ambermd.org/AmberModels.php) for the each molecule type including each lipid component must be chosen. Various members of the Amber community have developed tutorials, which are listed here. Common methods to build a membrane system include using: PACKMOL-Memgen, CHARMM, and Maestro.
## Lipid Force Fields
 Please see the force fields page to learn about [Lipid Force Fields.](https://ambermd.org/AmberModels_lipids.php)
 ## Tools and Tutorials to Build Membrane Systems
   [**Building Systems with CHARMM-GUI**](https://ambermd.org/tutorials/CHARMM-GUI.php)


 Some people find it easier to set up a system for simulation in Amber using the CHARMM-GUI interface. These tutorials were written by Wonpil Im at Lehigh University. He also has YouTube videos on these tutorials.
 [**Using PACKMOL-Memgen to Setup and Run a Membrane Simulation**](https://ambermd.org/tutorials/advanced/tutorial38/index.php)


 This tutorial will show you alternative ways of building a membrane system to perform molecular dynamics simulations within Amber, with or without an embedded protein, and with the optional inclusion of ligands. This is not a trivial task, as membrane environments are anisotropic (i.e. they are not the same in every dimension of the system), and, compared to typical explicit water simulations, the initial condition is important to obtain a representative structure in a reasonable simulation time. By Stephan Schott Verdugo
 [**Lipi21 and PACKMOL-Memgen Setup**](https://github.com/callumjd/AMBER-Membrane_protein_tutorial)


 PACKMOL-Memgen is built into Amber. This tutorial uses the CCG MOE GUI for protein preparation (add hydrogens, cap termini, set protonation states). You may have access to MOE or a similar GUI/tool to achieve the same task. PACKMOL-Memgen also has prep options. By JD Callum
 [**Lipid-building tutorial using Desmond/Maestro, PyMOL**](https://github.com/ParkerdeWaal/AMBER-Maestro-lipid-tutorial)


 This tutorial is a ptocol to build protein-lipid systems using the Maestro bundled with Desmond. The Maestro system has a high level of control of system size and orientation, allowing for simulation optimization. By Parker de Waal
  [**An Amber Membrane Simulation Tutorial: Lipid14**](https://ambermd.org/tutorials/advanced/tutorial16/index.php)  

 This tutorial explains how to set up and simulate lipid bilayers with the Lipid14 force field and CHARMM-GUI. This is a more general tutorial that walks you through how to build a system, equilibrate it, run MD and monitor properties. By Benjamin D. Madej and Ross C. Walker


