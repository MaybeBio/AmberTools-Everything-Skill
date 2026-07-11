

# Calculating Salt Molarity in an Explicit Water System
 By Abigail Held1 and Maria Nagan2


 1from Bill Miller III's lab, Truman State University, 2Stony Brook University
 # Learning Outcomes
 - Obtain the volume from LEaP
- Calculate the number of buffer ions to add in LEaP based on an overall molarity
 # Introduction
 Ion concentration and ionic strength need to be properly set for MD simulations. The concentration of NaCl in cells is approximately 150 mM, so we will be calculating how may sodium and chloride ions are required for this concentration in the protein system from [Tutorial 1.4](https://ambermd.org/tutorials/basic/tutorial7/index.php). By convention, these ions are added beyond what is required to neutralize the solute.


 

 1. Find the Volume of the box from leap.log file. 



Open the leap.log file generated from the *tleap* script without a buffer present for the volume of the box in Å3. The relevant portion of the leap.log output for the 2YX8 system from tutorial 1.4 is shown below. You can technically calculate the volume of the truncated octahedral box but LEaP also does this for you. Find the line with "Volume".




 ```
> solvateOct ramp OPCBOX 10.0
Scaling up box by a factor of 1.413369 to meet diagonal cut criterion
  Solute vdw bounding box:              45.400 31.108 30.346
  Total bounding box for atom centers:  73.667 73.667 73.667
      (box expansion for 'iso' is  55.9%)
  Solvent unit box:                     18.865 18.478 19.006
The number of boxes:  x= 4  y= 4  z= 4
  Volume: 208141.839 A^3 (oct)
  Total mass 110188.212 amu,  Density 0.879 g/cc
  Added 5569 residues.

```
 

 The volume of the truncated octahedron is 208141.839 Å3. 



2. Convert the volume of the system in Å3 to liters using the following conversion factors.




 Be sure to distribute the exponents when performing the calculations.




 ![A^3 to liters](https://ambermd.org/tutorials/basic/tutorial8/include/Molarity_Liters_calculation.PNG)
 

 3. Determine how many chloride ions are present in one liter of solution at a concentration of 150 mM.




 ![Cl atoms per liter](https://ambermd.org/tutorials/basic/tutorial8/include/Molarity_Clatoms_calculation.PNG)
 

 4. Determine how many chloride ions are needed in the system. 



Multiply the volume of the box by the concentration of Cl- ions and round to the nearest whole number.




 ![Cl ions to add](https://ambermd.org/tutorials/basic/tutorial8/include/Molarity_Clions.PNG)
 

 5. Determine the number of sodium ions needed for 150 mM NaCl. 



The amount of chloride ions and sodium ions needed are the same. If 19 Cl- ions are added, then 19 Na+ ions should also be added to the system. 

 

 Note 1: Upon equilibration of the system, the volume of the box will decrease as the solute-solvent interactions come into effect. Therefore, the equilibrated salt concentration is larger than the starting cocentration (150 mM NaCl) for NPT simulations.


 

 Note 2: [Matias Machado](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953) has proposed a possibly better way to calculate charges and buffers, the instructions for which can be found on [this page](http://archive.ambermd.org/202002/0194.html) of the AMBER mailing list archive. This tutorial was last updated July 16, 2021.


