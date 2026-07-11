

## Creating Stable Systems Through Relaxation
 After building the system, you will then relax the system through a process of minimization, heating, and reducing restraints gradually. This is used to ensure that the system is stable and representative of a real environment. This was previously called equilibration. 

 

  [**Relaxation of Explicit Water Systems**](https://ambermd.org/tutorials/basic/tutorial13/index.php)


 This tutorial demonstrates how to relax the system in explicit water solvents such as TIP3P and OPC. 

 

  [**Relaxation of Implicit Solvent System (GB)**](https://ambermd.org/tutorials/basic/tutorial15/index.php)


 Implicit solvent increases the speed of MD simulations by reducing the number of atoms in the simulation. This tutorial shows how to relax systems in implicit solvent such as the GB model.
 ## Running MD
 The production run is where data is collected. As part of your production run, you will create an input file that sets the parameters for the time, temperature, and pressure of the environment you would like to simulate your structure in. During the production run, Amber will save snapshots of your structure (called frames) at a user-given iteration. These frames will be linked together to create a movie of the changes in your structure (called a trajectory). This trajectory can then be analyzed. 

 

  [**Running MD with pmemd**](https://ambermd.org/tutorials/basic/tutorial14/index.php)


 This tutorial describe the elements needed for a production run with a protein system in explicit solvent. It also offers suggestions on how to analyze the resulting data. 

 

  [**Running MD in parallel**](https://ambermd.org/tutorials/basic/tutorial21/index.php)


 This tutorial describes how to run multiple independent molecular dynamics (MD) simulations in parallel using pmemd.MPI or pmemd.cuda.MPI.

