

# Developing Nonstandard Parameters
   ### Each entity in the simulation or calculation must have parameters.
 **Standard Biomolecules,Common Ions and Water:** See the [Force Field page](https://ambermd.org/AmberModels.php). 

 

 **Metal Ions (especially transition metals) and Nonstandard Residues:** The following tutorials give examples on how to parameterize systems that are complicated.  [**Simulating a pharmaceutical compound using antechamber and the Generalized Amber Force Field**](https://ambermd.org/tutorials/basic/tutorial4b/index.php)


  Antechamber is a set of tools in Amber that can be used to prepare input files for organic molecules, which can then be read into LEaP and used to create prmtop and inpcrd files. The Antechamber suite is designed for use with the [**G**eneral **A**MBER **F**orce **F**ield (GAFF)](http://ambermd.org/antechamber/gaff.html) and is ideal for setting up simulations involving organic pharmaceutical compounds or other organic molecules. In this tutorial we will use antechamber to create a leap input file for BMS's HIV reverse transcriptase inhibitor sustiva (efavirenz). Then we set up a simulation of sustiva bound to HIV-RT. [Japanese translation](http://amber.tkanai-lab.org/)) 



  [**Original antechamber tutorials**](https://ambermd.org/antechamber/example.html)


 These are a set of original antechamber tutorials. They contain examples of how to generate topology files for normal molecules, non-standard amino acid resues, and nucleotides in gas simulations. 



  [**Setting up a DNA-Ligand System**](https://ambermd.org/tutorials/advanced/tutorial1/index.php)


  This tutorial covers setting up an advanced system. In this case it shows you how to set up a dye system that is covalently bound to DNA. It also includes manually running multiconformational RESP fits, building custom units and assigning parameters manually. 



  [**Simulating the Green Fluorescent Protein and building a modified amino acid residue**](https://ambermd.org/tutorials/basic/tutorial5/index.php)


 This tutorial prepares a system with a modified amino acid. It differs from other examples that parametrize small organic compounds because the custom residue is part of a polymer chain. 

 

  [**Metal Ion Modeling Tutorial**](https://ambermd.org/tutorials/advanced/tutorial20/index.php)


 In this tutorial we will delineate several modeling strategies of metal ions in mixed systems (proteins and nucleic acids) using the AmberTools package. Both the bonded model and nonbonded model are illustrated. For the bonded model, MCPB and MCPB.py are used to facilitate the modeling. While for the nonbonded model modeling strategies for the 12-6 Lennard-Jones (LJ) and 12-6-4 LJ-type nonbonded models are presented. ( [Japanese translation](http://amber.tkanai-lab.org/)) 



  [**Electrostatic Parameterization with PyRESP**](https://ambermd.org/tutorials/basic/tutorial19/index.php)


 In this tutorial we will understand the limitations of additive force fields and the need for polarizable models like pGM-ind and pGM-perm. Then this tutorial will show you the steps involved in the parameterization process of the PyRESP algorithm. 



  [**RESP Charge Derivation with PyPE_RESP**](https://ambermd.org/tutorials/basic/tutorial22/index.php)


 PyPE_RESP is a command-line pipeline that automates the entire two-stage RESP fitting procedure through a single input file. It supports charge constraints via atom indices and SMARTS patterns, multi-molecule RESP fitting with intermolecular constraints, and cross-molecule atom equivalencies. 

   [**Deriving Implicitly Polarized Charges in `mdgx`**](https://ambermd.org/tutorials/advanced/tutorial28/index.php)


  This example will guide users through the process of making implicitly polarized charges for glycerol, appropriate for simulations in liquid water. This functionality in the `mdgx` program offers a self-contained and highly adapatable way for users to create charges tailored for specific environments and understand the level of accuracy. The same procedures that make IPolQ charges can be leveraged to perform traditional ESP fitting. 



  [**Deriving custom force field parameters with `mdgx`**](https://ambermd.org/tutorials/advanced/tutorial32/index.php)


  This example showcases the `mdgx` valence parameter fitting capabilities, taking the glycerol from the previous tutorial and also including a more complicated diol. Parameters are derived in a streamlined, highly automated procedure that puts users firmly in control of the molecular model building. Generational learning improves the outcome and ensures that the model can guide simulations while maintaining agreement with its quantum benchmark. 



  [**Adding custom extra points to a model**](https://ambermd.org/tutorials/advanced/tutorial35/index.php)


  This example showcases an expanded extra point palette debuted in Amber22. Models enhanced with customized extra points will run in either `pmemd` or `pmemd.cuda` based on topologies modified by `mdgx`. More accurate electrostatic and hydrogen bonding models are within reach of Amber simulations. 



   ### Deprecated
 [**Building your own Custom Residues**](https://ambermd.org/tutorials/advanced/tutorial1_orig/index.php) (old version)


  This tutorial is somewhat replaced by the tutorials directly above and below. However, it is kept here since it shows our original method for creating a system containing nonstandard residues, in this case a metal atom. A new unit is generated in `xleap` in order to simulate the plastocyanin system with a bound copper ion. There are two versions of this tutorial. The first is a [simple version](https://ambermd.org/tutorials/advanced/tutorial1_orig/index.php) which creates just a new copper residue and approximates it as a +1 ion. There is also a more [advanced version](https://ambermd.org/tutorials/advanced/tutorial1_adv/index.php) where new special histidine and methionine residues are created so that different charges and bond / angle and dihedral parameters can be used. 



  [**[Deprecated] Generating Force Field Parameters with Paramfit**](https://ambermd.org/tutorials/advanced/tutorial23/index.php)


  In this tutorial we will generate force field parameters for two small molecules from *ab-initio* quantum calculations using the AmberTools program Paramfit. This tutorial generates the phi and psi dihedral potentials over two different small peptide chains, and details each step of the parameter generation process from preparation of a conformational sampling of each structure to generation of quantum data to evaluating the quality of the resulting parameters.


