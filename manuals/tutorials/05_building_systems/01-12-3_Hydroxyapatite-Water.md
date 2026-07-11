

Copyright (C) Pengfei Li 2016
 ## Setting Up A Hydroxyapatite Slab in Water Box
  ![](https://ambermd.org/tutorials/advanced/tutorial27/images/HAP_water.png)
 Some references about the apatite force field in the INTERFACE force field could be found here. [[1]](#ref1) [[2]](#ref2) [[3]](#ref3)
 1. Generate a monomer structure:
 > You can find a hydroxyapatite unit cell car file (hap_unit_cell.car) under the "INTERFACE_FF_1_5/MODEL_DATABASE/HYDROXYAPATITE" directory, which is from the INTERFACE force field version 1.5 you can download from: [INTERFACE-MD](https://bionanostructures.com/interface-md/). It has six PO43- groups, two OH- groups, and ten Ca2+ ions in the unit. We can check the car file and get that each PO43- group has a partial charge of -2.2e, each OH- group has a partial charge of -0.9e, while each Ca2+ has a partial charge of +1.5e, yielding a total charge of zero for the unit cell.
>  In order to build related systems in LEaP, we need to have related PDB file and mol2/prepi/lib files available for the system.
>  We need to generate a PDB and a mol2 file for the HAP unit cell. And the generated mol2 file should have the atom types of the AMBER compatible INTERFACE force field. Herein we will use the car_to_files.py program in pyMSMT to generate the related PDB and mol2 files. The program is available in AmberTools17 but not in AmberTools16.
>  You can use the car_to_files.py program to convert the car file into a mol2 file and a PDB file by using the following command:
>  car_to_files.py -i hap_unit_cell.car -m hap.mol2 -p hap.pdb -r HAP
>  It will generate a [hap.mol2](https://ambermd.org/tutorials/advanced/tutorial27/files/hap.mol2) file and a [hap.pdb](https://ambermd.org/tutorials/advanced/tutorial27/files/hap.pdb) with residue name as HAP. Due to car_to_files.py will re-assign the bond information itself, sometimes the bond information may not be correct, it is important to check and fix it by using VMD. Add/Removing bonds can be done through "Mouse"->"Add/Remove Bonds" option in the "VMD Main" control panel, and afterwards users can save a new mol2 file by right click the molecule in the "VMD Main" window, and click "Save Coordinates", in "Selected atoms" choose "all", and in the "File type" choose "mol2", and then click "Save..." to save the mol2 file. It is fine for the bond linkages for the mol2 file herein through visualizing the mol2 file using VMD, so we don't need to do any changes.
 1. Build the water solvated hydroxyapatite cluster:
 > Herein we use the PropPDB to get a bigger hydroapatite cluster by duplicating the unit cell.
>  PropPDB -p hap.pdb -o hap16.pdb -ix 4 -iy 4 -iz 1
>  Here is the hap16.pdb file generated: [hap16.pdb](https://ambermd.org/tutorials/advanced/tutorial27/files/hap16.pdb). Afterwards we can use the tleap program to build the system, here is the input file: [hap_tleap.in](https://ambermd.org/tutorials/advanced/tutorial27/files/hap_tleap.in), in which we use the SPC/E water model to model the solvents. Then we run the command:
>  tleap -s -f hap_tleap.in > hap_tleap.out
>  Then we will get the [hap16_water.prmtop](https://ambermd.org/tutorials/advanced/tutorial27/files/hap16_water.prmtop) and [hap16_water.inpcrd](https://ambermd.org/tutorials/advanced/tutorial27/files/hap16_water.inpcrd) files, which can be used to do the minimization and MD simulations.
 1. Perform minimization and MD simulations:
 > Herein we will perform the minimization and MD simulations under liquid phase, therefore we can take advantage of the steep of pmemd.cuda to perform these modelings. To just show an example, in general herein we only performed four steps: minimization of water, full minimization, heating, and production. However, practical research may need much more steps to get meaningful data.
>  Here are the three input files for different stages: [hap16_water_min1.in](https://ambermd.org/tutorials/advanced/tutorial27/files/hap16_water_min1.in), [hap16_water_min2.in](https://ambermd.org/tutorials/advanced/tutorial27/files/hap16_water_min2.in), [hap16_water_heat.in](https://ambermd.org/tutorials/advanced/tutorial27/files/hap16_water_heat.in), and [hap16_water_prod.in](https://ambermd.org/tutorials/advanced/tutorial27/files/hap16_water_prod.in)
>  And then we can perform the following commands sequentially to perform these simulations:
>  pmemd.cuda -O -i hap16_water_min1.in -o hap16_water_min1.out -p hap16_water.prmtop -c hap16_water.inpcrd -r hap16_water_min1.rst -ref hap16_water.inpcrd
>  pmemd.cuda -O -i hap16_water_min2.in -o hap16_water_min2.out -p hap16_water.prmtop -c hap16_water_min1.rst -r hap16_water_min2.rst
>  pmemd.cuda -O -i hap16_water_heat.in -o hap16_water_heat.out -p hap16_water.prmtop -c hap16_water_min2.rst -r hap16_water_heat.rst -x heat.netcdf -ref hap16_water_min.rst
>  pmemd.cuda -O -i hap16_water_prod.in -o hap16_water_prod.out -p hap16_water.prmtop -c hap16_water_heat.rst -r hap16_water_prod.rst -x prod.netcdf -ref hap16_water_heat.rst
 [Click here to go back to the main page](https://ambermd.org/tutorials/advanced/tutorial27/index.php)
 **[1]** Hendrik Heinz, Tzu-Jen Lin, Ratan K. Mishra, and Fateme S. Emami, "Thermodynamically consistent force fields for the assembly of inorganic, organic, and biological nanostructures: The INTERFACE force field", *Langmuir*, **2013**, 29, pp. 1754-1765
 **[2]** Tzu-Jen Lin, "Force field parameters and atomistic surface models for hydroxyapatite and analysis of biomolecular adsorption at aqueous interfaces", *PhD Thesis, The University of Akron*, **2013**
**[3]** Hendrik Heinz, "The role of chemistry and pH of solid surfaces for specific adsorption of biomolecules in solution-accurate computational models and experiment", *Journal of Physics: Condensed Matter*, **2014**, 26, pp. 244105

