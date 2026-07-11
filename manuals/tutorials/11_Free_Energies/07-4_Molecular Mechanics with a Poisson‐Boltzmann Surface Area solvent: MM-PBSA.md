

# MM-PBSA
 ![](https://ambermd.org/tutorials/advanced/tutorial3/images/ras-raf.jpg)
 ## Table of Contents
 - [Section 1](https://ambermd.org/tutorials/advanced/tutorial3/section1.php) : Build the starting structure and run a simulation to obtain an equilibrated system.
- [Section 2](https://ambermd.org/tutorials/advanced/tutorial3/section2.php) : Run the production simulation and obtain an ensemble of snapshots.
- [Section 3.1](https://ambermd.org/tutorials/advanced/tutorial3/py_script/section1.php) : Calculate the binding free energy of a protein-protein complex (Ras-Raf).
- [Section 3.2](https://ambermd.org/tutorials/advanced/tutorial3/py_script/section2.php) : Calculate the binding free energy of a protein-ligand complex (Estrogen Receptor and Raloxifene).
- [Section 3.3](https://ambermd.org/tutorials/advanced/tutorial3/py_script/section3.php) : Calculate the binding free energy of Ras-Raf and use Alanine Scanning to compare to the binding energy of a mutant Ras-Raf complex that has had a residue mutated to alanine and analyze the results.
- [Section 3.4](https://ambermd.org/tutorials/advanced/tutorial3/section4.php) : Calculate the binding free energy of Ras-Raf in parallel using three processors.
- [Section 3.5](https://ambermd.org/tutorials/advanced/tutorial3/section5.php) : Calculate the entropy of the Estrogen Receptor and Raloxifene complex using Normal Mode Analysis (Nmode).
- [Section 3.6](https://ambermd.org/tutorials/advanced/tutorial3/section6.php) : Decomposing the free energy contributions to the binding free energy of Ras-Raf in a per-residue or pairwise per-residue basis.
 ## Introduction
 In this tutorial we will use the MM-PBSA method to calculate the binding free energy for the association of two proteins.
 The overall objective of the MM-PBSA method and it's complementary MM-GBSA method is to calculate the free energy difference between two states which most often represent the bound and unbound state of two solvated molecules or alternatively to compare the free energy of two different solvated conformations of the same molecule.
 ![](https://ambermd.org/tutorials/advanced/tutorial3/images/equation1.gif)
 Ideally we would like to calculate this free energy of binding directly as shown in the figure below:
 ![](https://ambermd.org/tutorials/advanced/tutorial3/images/figure1.gif)
 However, in such a simulation of these solvated states the majority of the energy contributions would come from solvent-solvent interactions and the fluctuations in total energy would be an order of magnitude larger than binding energy. Thus the calculation would take an inordinate amount of time to converge. Thus a more effective method is to divide up the calculation according to the following thermodynamic cycle:
 ![](https://ambermd.org/tutorials/advanced/tutorial3/images/figure2.gif)
 Evidently from this diagram the binding free energy delta-Gbind,solv can be calculated by:
 ![](https://ambermd.org/tutorials/advanced/tutorial3/images/equation2.gif)
 In the MM-PBSA approach the different contributions to the binding free energy above are calculated in various ways:
  In this tutorial we will demonstrate the use of the MM/PB(GB)SA scripts included with Amber and AmberTools to automatically perform all the necessary steps to estimate the binding free energy of a protein-protein complex (RAS and RAF) and a protein-ligand complex (Estrogen Receptor and Raloxifene) using both MM-GBSA and MM-PBSA methods in serial and parallel. Furthermore, we will be demonstrating the use of Alanine Scanning and Normal Mode entropy calculations using the script. In principle, the calculation of the binding free energy described above would require three independent MD simulations of the complex and both individual proteins. However, typically one makes the approximation that no significant conformational changes occur upon binding so that the snapshots for all three species can be obtained from a single trajectory. This is the 'single trajectory approach' and is what we will use in this tutorial.
 Perl Version by Ross Walker and Thomas Steinbrecher 

Python Version by Dwight McGee, Bill Miller III, and Jason Swails


