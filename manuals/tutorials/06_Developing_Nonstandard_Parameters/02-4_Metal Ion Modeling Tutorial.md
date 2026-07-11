

# Metal Ion Modeling Overview
 ## Learning Outcomes
 In this tutorial we will delineate several modeling strategies of metal ions in mixed systems (proteins and nucleic acids) using the AmberTools package. As the user, you should choose the best model for your simulation needs.
- [Bonded Model](#bonded)
- [Nonbonded Model](#nonbonded)
- [Cationic Dummy Atom Model](#cation)
 ### Resources
 - A comprehensive review about metal containing system modeling has been published, [[1]](#ref1) which covers various strategies in quantum modeling, classical force field methods, and polarizable force field models.
- Readers can consult the paper for more details: [Metal Ion Modeling Using Classical Mechanics](http://pubs.acs.org/doi/abs/10.1021/acs.chemrev.6b00440), which is open accessed (free to download). Meanwhile, this review is also a good resource to learn about the basics and history of force field.
- If you are new to Amber, strongly suggest you to read [the LEaP tutorial](https://ambermd.org/tutorials/pengfei/index.php) first, which covers some basic concepts of modeling in Amber.
 ## Introduction
 Approximately, 40% of the protein structures contained in the Protein Data Bank have metal ions. They play significant roles in the structural, catalytic and electron transfer processes involved in the metabolism of an organism.
 In Molecular Dynamics (MD) simulation, there are multiple ways to model metal ions in protein systems. For example, the bonded model, [[2]](#ref2) nonbonded model, [[3]](#ref3) cationic dummy atom model, [[4]](#ref4) etc. have all been described. The bonded model treats the interaction between metal ions and its ligating residues via bond, angle, dihedral, electrostatic and van der Waals (VDW) terms. These parameters can be readily obtained from QM calculation or where available experimental results. The nonbonded model treats the metal ion-protein interaction via electrostatic and VDW terms only. The cationic dummy atom model places the charges between the metal ions and surrounding residues to mimic covalent bonds and offer a more sophisticated electrostatic model.
 Each of these models have their associated pluses and minuses. For example, the bonded model cannot simulate coordination number (CN) changes and ligand exchange processes. The nonbonded model qualitatively simulates the covalent bond and the electron transfer effect due to its lack of bonded terms and the integer charge on the metal sites. The cationic dummy atom model requires an intricate parameterization process due to the many empirical parameters needed. For example, not only the charge and VDW parameters on the metal ion are needed, but also the charge and VDW on each site as well.
 ![](https://ambermd.org/tutorials/advanced/tutorial20/images/Metal_site.png)
 Here we present several tutorials about different modeling strategies (both bonded and nonbonded models) for user's need:
 ## Bonded Model
 While ligand exchange of CN switching process certainly do happen in protein metal sites, these occur on much slower timescales than what is observed for metal ions in aqueous solution. For example, for solvated metal ion complexes, the residence time of the first solvation shell water molecules are in the range of 5×10-4 microsecond (μs) to 1.3×104 μs for +2 metal ions, and in the range of 0.05 μs to 3.2×1013 μs for +3 metal ions. The residence time of ligands in the first coordination shell of a metal site in many proteins is far longer due site pre-organization among other factors. Since typical MD simulation timescales are at the nanosecond or microsecond levels (circa 2014), the bonded model is an effective approach to study processes where ligand exchange does not occur over the timescales employed.
 ### 1. Metal Center Parameter Builder(MCPB)
 In proteins (which contain 20 standard amino acids) and nucleic acids (contains 5 standard nucleosides), there are standard libraries and parameter files in the Amber force field. The metal sites in proteins vary quite a bit, so it is very hard to parameterize a standard force field for them. Taking these considerations into account, Dr. Martin Peters in Merz research group developed the Metal Center Parameter Builder (MCPB) [[2]](#ref2) in the Modeling ToolKit++ (MTK++) software package in AmberTools. It a Semi-automated workflow, which facilitates the bonded model parameterization process for metal sites in protein systems.
 ### 2. MCPB.py
 MCPB.py is a python version of MCPB. It uses an optimized workflow with a modest number of steps to perform the modeling. [[5]](#ref5) It is available only in AmberTools 15 or higher. The current version supports more than 80 metal ions and various Amber force fields. It was developed by Pengfei Li in the Merz research group in Michigan State University.
 [Click here to go to: Ion Modeling using MCPB/MCPB.py](https://ambermd.org/tutorials/advanced/tutorial20/bonded_model.php)
 ## Nonbonded Model
 ### 1. 12-6 Lennard-Jones (LJ) Nonbonded Model
 The 12-6 LJ nonbonded model is widely used due to its simple form and excellent transferability. Li et al. have parameterized the 12-6 model for more than 60 ions spanning from monovalent to tetravalent across the periodic table. [[3]](#ref3) [[6]](#ref6) [[7]](#ref7) These parameters are available in current the Amber force field (Please check the **Amber force fields-->Molecular mechanics force fields-->Ions** section in the manual for more details).
 ### 2. 12-6-4 LJ-Type Nonbonded Model
 The 12-6-4 LJ-type nonboned model was proposed and parameterized for divalent metal ions by Li and Merz. [[8]](#ref8) A C4 term was added to the 12-6 LJ nonbonded model to represent the ion-induced dipole interaction, which is proportional to r-4. Later on Li et al. parameterized the 12-6-4 model for monovalent, [[6]](#ref6) trivalent [[7]](#ref7) and tetravalent [[7]](#ref7) ions. It was shown that the 12-6-4 model simultaneously reproduced several experimental values (HFE, IOD and CN) coupled with excellent transferability properties for mixed systems. [[6]](#ref6) [[7]](#ref7) [[8]](#ref8) In a recent paper of Panteva et al. they found that the 12-6-4 model combined with the SPC/E water model preformed the best for Mg(II) among 17 different nonbonded models investigated. [[9]](#ref9)
 ### 3. Modified 12-6-4 LJ-Type Nonbonded Model
 The modified 12-6-4 LJ-type nonboned model was proposed and parameterized for ion-ligand interactions. Polarizability of the ligand was modified so the final C4 values will be different for ligands like imidazole [[10]](#ref10) or acetate [[11]](#ref11).
 [Click here to go to: Ion Modeling using Nonbonded Model](https://ambermd.org/tutorials/advanced/tutorial20/nonbonded_model.php)
 [To help the users select correct C4 values, we provided a C4 lookup table so the users do not need to calculate from scratch using the polarizability values](https://ambermd.org/tutorials/advanced/tutorial20/c4table.php)
 ## Cationic Dummy Atom Model
 Users can review the tutorial here:
 [Click here to go to: Tutorial for the Cationic Dummy Atom Model for the Zinc ion](http://mayoresearch.mayo.edu/mayo/research/camdl/zinc_protein.cfm)
  References:
 **[1]** Pengfei Li and Kenneth M. Merz, Jr. "Metal Ion Modeling Using Classical Mechanics", *Chem. Rev.*, **2017**, 117, 1564-1686
 **[2]** Martin B. Peters, Yue Yang, Bing Wang, László Füsti-Molnár, Michael N. Weaver, and Kenneth M. Merz, Jr. "Structural Survey of Zinc-Containing Proteins and Development of the Zinc Amber Force Field (ZAFF)", *J. Chem. Theory Comput.*, **2010**, 6, 2935-2947
 **[3]** Pengfei Li, Benjamin P. Roberts, Dhruva K. Chakravorty, and Kenneth M. Merz, Jr. "Rational Design of Particle Mesh Ewald Compatible Lennard-Jones Parameters for +2 Metal Cations in Explicit Solvent", *J. Chem. Theory Comput.*, **2013**, 9, 2733-2748
 **[4]** Yuan-Ping Pang "Successful molecular dynamics simulation of two zinc complexes bridged by a hydroxide in phosphotriesterase using the cationic dummy atom method", *Proteins: Struct., Funct., Bioinf.*, **2001**, 45, 183-189
 **[5]** Pengfei Li and Kenneth M. Merz, Jr. "MCPB.py: A Python Based Metal Center Parameter Builder." *J. Chem. Inf. Model.*, **2016**, 56, 599-604
 **[6]** Pengfei Li, Lin F. Song and Kenneth M. Merz, Jr. "Systematic Parameterization of Monovalent Ions Employing the Nonbonded Model", *J. Chem. Theory Comput.*, **2015**, 11, 1645-1657
 **[7]** Pengfei Li, Lin F. Song and Kenneth M. Merz, Jr. "Parameterization of Highly Charged Metal Ions Using the 12-6-4 LJ-Type Nonbonded Model in Explicit Water", *J. Phys. Chem. B*, **2015**, 119, 883-895
 **[8]** Pengfei Li and Kenneth M. Merz, Jr. "Taking into Account the Ion-Induced Dipole Interaction in the Nonbonded Model of Ions", *J. Chem. Theory Comput.*, **2014**, 10, 289-297
 **[9]** Maria T. Panteva, George M. Giambasu, Darrin M. York "Comparison of Structural, Thermodynamic, Kinetic and Mass Transport Properties of Mg2+ Ion Models Commonly Used in Biomolecular Simulations", *J. Comput. Chem.*, **2015**, 13, 970-982
 **[10]** Zhen Li, Lin Frank Song, Gaurav Sharma, Basak Koca Fındık, and Kenneth M. Merz Jr. "Accurate Metal–Imidazole Interactions" *J. Chem. Theory Comput.*, **2023**, 19(2), 619-625, DOI: 10.1021/acs.jctc.2c01081
 **[11]** Majid Jafari, Zhen Li, Lin Frank Song, Luca Sagresti, Giuseppe Brancato, and Kenneth M. Merz Jr., "Thermodynamics of Metal–Acetate Interactions" *J. Phys. Chem. B*, **2024**, 128(3), 684-697, DOI: 10.1021/acs.jpcb.3c06567
 by Pengfei Li, Majid Jafari, Zhen Li & Kenneth M. Merz Jr. 2024


