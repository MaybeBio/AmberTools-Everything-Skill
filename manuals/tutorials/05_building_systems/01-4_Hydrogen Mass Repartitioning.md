

# Hydrogen Mass Repartitioning
 Implementation of methods found in: 

 Hopkins, C.W., Le Grand, S.L., Walker, R.C. and A.E. Roitberg, Long-Time-Step Molecular Dynamics through Hydrogen Mass Repartitioning, [*J. Chem. Theory Comput.*, **2015**, *11(4)*, 1864-1874.](https://pubs.acs.org/doi/abs/10.1021/ct5010406) 

 

 ![Hydrogen Mass Repartitioning](https://ambermd.org/tutorials/basic/tutorial12/include/hmass.gif)


 Figure 1: *Hydrogen Mass Repartitioning example with Alanine Dipeptide* (a) Dialanine peptide atomic masses. (b) Dialanine peptide masses after repartitioning. Notice that the overall mass remains the same.
 

 ## Learning Outcomes
 Understand the reasoning behind repartitioning the mass of a molecule to make hydrogen heavier
 Be able to use **parmed** to do hydrogen mass repartitioning on an existing topology file (parm7)
 ## Introduction
 Hydrogen mass repartitioning (HMR) is very useful because it allows a larger time step to be used for your simulation by redistributing some of the mass from heavy atoms connected to hydrogen into the bonded hydrogens. This means that simulations can accurately represent a longer time frame without encountering instability-related errors caused by high-frequency hydrogen motion. Some biological studies are interested in large motions of biomolecules. These require longer simulation times, which can be accomplished by using hydrogen mass repartitioning. For instance, if one wanted to run a simulation with a 4 fs time step instead of a more customary 2 fs time step, then one could use parmed (below) to perform hydrogen mass repartitioning on an existing parm7 topology file and then run equilibration and MD.




 For more information about hydrogen mass repartitioning, please read [Roitberg and coworkers](https://pubs.acs.org/doi/abs/10.1021/ct5010406). ## Process
 1. Use parmed to load a topology file built in LEaP


 

 For this tutorial, we will use the topology file made during the [5.1 Simple Simulation of Alanine Dipeptide](https://ambermd.org/tutorials/basic/tutorial0/index.php#Save_the_Amber_parm7_and_rst7_input_) tutorial.




 The alanine dipeptide **parm7** file is available here:


 [diala.parm7](https://ambermd.org/tutorials/basic/tutorial12/include/diala.parm7)


 ```
 parmed diala.parm7ParmEd: a Parameter file EditorLoaded Amber topology file parm7Reading input from STDIN...>
```
 2. Perform hydrogen mass repartitioning with the *hmassrepartition* command 

 

 This command both triples the mass of all hydrogens on alanine dipeptide and scales down the mass of all other atoms on alanine dipeptide, preserving the total mass (See Figure 1). To learn more about the *hmassrepartition* command, please look through the **HMassRepartition** subsection on page **271** of the [Amber 2021 Manual](https://ambermd.org/doc12/Amber21.pdf#page=285). ```
 hmassrepartitionRepartitioning hydrogen masses to 3.024 daltons. Not changing water hydrogen masses.>
```
 3. Save the new topology file




 Give this file a new name that says " *hmass*" so you know it was modified. Close parmed with the *quit* command. ```
 outparm diala_hmass.parm7Outputting Amber topology file diala_hmass.parm7> quitDone!
```
 

 The new file hmass.parm7 should look like this [diala_hmass.parm7](https://ambermd.org/tutorials/basic/tutorial12/include/diala_hmass.parm7)


 By Jan Ziembicki, Carlos Simmerling and Maria Nagan


