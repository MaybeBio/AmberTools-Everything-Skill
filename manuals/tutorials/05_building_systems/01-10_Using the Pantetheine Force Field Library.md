

# Using the Pantetheine Force Field library
 By Shiji Zhao. Adapted for Amber website by Elizabeth Sebastian.
 ## 1. Introduction
 The Pantetheine Force Field (PFF) library is an Amber-compatible force field library for various pantetheine-containing ligands (PCLs). The PFF library was parameterized using the Gasteiger, AM1-BCC, or RESP charging methods in combination with the gaff2 and ff14SB parameter sets. In this tutorial, we are going to demonstrate how to use the parameter files in the PFF library for MD simulations using the `AMBER` molecular dynamics engines on two protein-PCL systems: Phosphopantetheine adenylyltransferase (PPAT) bound with phosphopantetheine (Ppant); and 3-hydroxy-3-methylglutaryl synthase/acyl carrier protein complex (HGMS/ACP) with phosphopantetheinyl serine (Ppant-Ser) nonstandard amino acid residue. Note that in the PPAT-Ppant complex, Ppant is a standalone ligand that binds to the protein non-covalently. While in the HGMS/ACP-Ppant-Ser complex, Ppant-Ser is part of the polymer chain of ACP. This difference will lead to slightly different parameterizations for the two types of complexes.
 Note: This tutorial also relies on the interactive visualization, preparation and analysis software [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/), but it is not necessary. All operations involving Chimera could be done with text editors.
 ## 2. Preparing the PPAT-Ppant System Using Ppant Parameter Files
 In this section, we will show the preparation process for MD simulations of phosphopantetheine adenylyltransferase bound with phosphopantetheine (PPAT-Ppant system) using the PFF parameter files of Ppant. Since Ppant is standalone ligand of the system, this tutorial also applies to the usage of PFF parameter files for CoA PCLs.




 We start with the PDB file [1od6.pdb](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/PNS/1od6.pdb), which contains the PPAT with a Ppant ligand. The PDB file could be downloaded from the [RCSB website](https://www.rcsb.org/structure/1OD6). Besides protein and PCL, the original PDB file contains solvent molecules such as water (HOH) and sulfate ion (SO4). You can use either UCSF Chimera or text editors to remove those unwanted molecules. [Here](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/PNS/1od6_clean.pdb) is the cleaned PDB file. 



 Another problem with the PDB file is that amino acid residues 38-42 are missed. This is a common problem of the flexible loop regions of X-Ray crystal protein structures. Fortunately, Chimera provides a missing loop modeling feature under **Tools -> Structure Editing -> Model/Refine Loops** by interfacing with `Modeller`, a powerful tool for homology or comparative modeling of protein structures. [Here](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/PNS/1od6_clean_fix.pdb) is the generated PDB file with the missing loop fixed.
 ### 2.2 Creating the topology and coordinate files of the PPAT-Ppant system
 `AMBER` built-in force field parameter sets include those for proteins, DNAs, RNS, carbohydrates, lipids, etc. In order to perform MD simulations on the PPAT-Ppant system, we need extra parameter files for the Ppant ligand. That is why you are using the PFF library! In this [page](https://github.com/ShijiZ/PFF-Library/tree/main/Ppant-lib/PNS), you can find the parameter files for Ppant charged with three different charging methods. For example, if you want to use the parameters charged with the AM1-BCC method, you need to download [PNS-BCC.lib](https://github.com/ShijiZ/PFF-Library/blob/main/Ppant-lib/PNS/AM1BCC/PNS-BCC.lib) and [PNS-BCC.frcmod](https://github.com/ShijiZ/PFF-Library/blob/main/Ppant-lib/PNS/AM1BCC/PNS-BCC.frcmod). 



 Now we have all the files needed to create the topology and coordinate files for sander or pmemd! We just need to load these files into LEaP to create these files. We choose the *ff19SB* field for the protein, the *gaff2* force field for the non-electrostatic terms of Ppant, and the OPC explicit water model. Then, we load the parameter files from the PFF library, and load the PDB file prepared in **section 2.1**. The molecular system is solvated in an octahedral box of OPC water molecules with thickness extending 8 Å from the protein surface, and is neutralized by adding or Cl- counter ions. Finally, we save the topology and coordinate files, as well as a PDB file for visual checking. This is shown in the script below, which we will save to [tleap_1od6.in](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/PNS/tleap/tleap_1od6.in).
 ```
source leaprc.protein.ff19SB 
source leaprc.gaff2
source leaprc.water.opc 

loadamberparams PNS-BCC.frcmod
loadoff PNS-BCC.lib

mol = loadpdb 1od6_clean_fix.pdb

solvateoct mol OPCBOX 8.0
addionsrand mol Cl- 0

saveamberparm mol 1od6-bcc.prmtop 1od6-bcc.inpcrd
savepdb mol 1od6-bcc.pdb

quit

```
 We can then run `tleap` with the following command:
 ```
tleap -f tleap_1od6.in
```
 After this step, you should have your finished topology and coordinate files, and you are ready to run simulations! You can download the [1od6-bcc.prmtop](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/PNS/tleap/1od6-bcc.prmtop), [1od6-bcc.inpcrd](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/PNS/tleap/1od6-bcc.inpcrd), and [1od6-bcc.pdb](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/PNS/tleap/1od6-bcc.pdb) files we created to compare with your own.
 ## 3. Preparing the HGMS/ACP-Ppant-Ser System Using Ppant-Ser Parameter Files
 In this section, we focus on the preparation of PDB files of HMG synthase (HGMS) in complex with acyl carrier protein (ACP) using the PFF parameter files of Ppant-Ser. Different from the PPAT-Ppant system where the Ppant ligand is standalone in the complex, the Ppant ligand in this case is covalently linked to a conserved serine of ACP, so there are additional complications involved with preparing the files. 



 We start with the PDB file [5kp7.pdb](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/PPS/5kp7.pdb), which contains the structure of HGMS/ACP complex. Chain B is the ACP, where Ppant 100 is covalently linked to serine 39. The PDB file could be downloaded from the [RCSB website](https://www.rcsb.org/structure/5KP7). Similar to [section 2](#2), we need to remove solvent molecules including water (HOH) and sulfate ion (SO4), and fix the missing loop using `Modeller`. [Here](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/PPS/5kp7_clean_fix.pdb) s the cleaned and fixed PDB file, in which serine 39 and Ppant 100 in the original PDB file has been renumbered as serine 457 (SER 457) and Ppant (PNS 498) during the loop modeling process. 



 Note that in the generated PDB file, SER 457 and PNS 498 are two independent residues. However, the Ppant-Ser parameter files in the PFF library treat the covalently linked Ppant and serine amino acid as a single residue, named as PNS. Therefore, an important step is to “stick” SER 457 and PNS 498 together into a new residue. This step has to be done manually. First, we cut and copy the coordinate information of PNS 498 right below those of SER 457. Second, we rename both SER 457 and PNS 498 to PNS 457. Now, these two separate residues have become a single residue named as PNS. [Here](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/PPS/5kp7_clean_fix_modify.pdb) is the modified PDB file. You may find the atoms numbers are out of order know. Don’t worry, `tleap` will fix this problem automatically in the next section.
 ### 3.2 Creating the topology and coordinate files of the HGMS/ACP-Ppant-Ser system
 Now we are ready to create the topology and coordinate files for subsequent MD simulations. Like in [section 2.2](#2.2), the parameter files for Ppant-Ser charged with three different charging methods could be found in this [page](https://github.com/ShijiZ/PFF-Library/tree/main/Ppant-Ser-lib/PNS). Note that there are two .frcmod files associated with each charging method, with names end with .frcmod1 and .frcmod2 respectively. This is due to the fact that non-electrostatic parameters of Ppant-Ser were first derived from the *ff14SB* force field (.frcmod1 file) then from the *gaff2* force field (.frcmod2 file), following this [tutorial](https://ambermd.org/tutorials/basic/tutorial5/index.php) of preparing parameter files for modified amino acid residues. For example, if you want to use the parameters charged with the AM1-BCC method, you need to download [S-PNS-BCC.lib](https://github.com/ShijiZ/PFF-Library/tree/main/Ppant-Ser-lib/PNS/AM1BCC/S-PNS-BCC.lib), [S-PNS-BCC.frcmod1](https://github.com/ShijiZ/PFF-Library/tree/main/Ppant-Ser-lib/PNS/AM1BCC/S-PNS-BCC.frcmod1), and [S-PNS-BCC.frcmod2](https://github.com/ShijiZ/PFF-Library/tree/main/Ppant-Ser-lib/PNS/AM1BCC/S-PNS-BCC.frcmod2). 



 Now we have all the files needed to create the topology and coordinate files for sander or pmemd! We just need to load these files into LEaP to create these files. We choose the *ff19SB* force field for the protein, the *gaff2* force field for the non-electrostatic terms of Ppant-Ser, and the OPC explicit water model. (Although Ppant-Ser in the PFF library was parameterized using *ff14SB* force field, the parameteres for Ppant-Ser are compatible with the *ff19SB* force field.) While loading the parameter files from the PFF library, it is very important that we load *S-PNS-BCC.frcmod2* file first, followed by *S-PNS-BCC.frcmod1* to make sure that all the *gaff2* parameters are overwritten by available *ff14SB* parameters. he molecular system is solvated in an octahedral box of OPC water molecules with thickness extending 8 Å from the protein surface, and is neutralized by adding or Na+ counter ions. Finally, we save the topology and coordinate files, as well as a PDB file for visual checking. This is shown in the script below, which we will save to [tleap_5kp7.in](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/PPS/tleap/tleap_5kp7.in).
 ```
source leaprc.protein.ff19SB 
source leaprc.gaff2
source leaprc.water.opc

loadamberparams PNS.frcmod2
loadamberparams PNS.frcmod1
loadoff PNS.lib 

mol = loadpdb 5kp7_clean_fix_modify.pdb

solvateoct mol OPCBOX 8.0
addionsrand mol Na+ 0

saveamberparm mol 5kp7-bcc.prmtop 5kp7-bcc.inpcrd
savepdb mol 5kp7-bcc.pdb

quit

```
 We can then run tleap with the following command:
 ```
tleap -f tleap_5kp7.in
```
 After this step, you should have your finished topology and coordinate files, and you are ready to run simulations! You can download the [5kp7-bcc.prmtop](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/PPS/tleap/5kp7-bcc.prmtop), [5kp7-bcc.inpcrd](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/PPS/tleap/5kp7-bcc.inpcrd), and [5kp7-bcc.pdb](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/PPS/tleap/5kp7-bcc.pdb), we created to compare with your own.
 ## 4. MD Simulations
 In this section we provide a sample workflow for MD simulations on either the PPAT-Ppant system prepared in [section 2](#2), or the HGMS/ACP-Ppant-Ser system prepared in [section 3](#3). In your own project, you are of course free to choose different minimization, heating, and equilibration routine than what we use here.
 ### 4.1 Minimization
 Following this [tutorial](https://ambermd.org/tutorials/basic/tutorial1/section5.php), we use a two-stage minimization. In the first stage we will keep the protein-PCL system fixed and just optimize the positions of the water and ions. Then in the second stage we will minimize the entire system. Below are the input files of the two stages: 



 [min1.in](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/MD/min1.in):
 ```
initial minimization solvent + ions
 &cntrl
    imin = 1,
    maxcyc = 1000,
    ncyc = 500,
    ntb = 1,
    ntr = 1,
    cut = 10.0
/
Hold PPAT-Ppant system fixed
500.0
RES 1 161
END
END

```
 Note: Change the residue numbers to 1-497 for HGMS/ACP-Ppant-Ser.
 [min2.in](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/MD/min2.in):
 ```
minimization of the whole system
&cntrl
    imin = 1,
    maxcyc = 2500,
    ncyc = 1000,
    ntb = 1,
    ntr = 0,
    cut = 10.0
/
```
 Run the following two commands in sequence to perform the two-stage minimizations:
 ```
pmemd -O -i min1.in -p *.prmtop -c *.inpcrd -o min1.out -r min1.rst -ref *.inpcrd
pmemd -O -i min2.in -p *.prmtop -c min1.rst -o min2.out -r min2.rst
```
 Note that the CPU program *pmemd* instead of the GPU program is used for minimization, for the sake of [accuracy consideration](https://ambermd.org/GPULogistics.php#Accuracy). You should visualize the resulting structure ( *min2.rst*) to make sure that nothing obviously bad happened, as well as the output file to make sure that everything looks OK (e.g., that the structure remained intact and that the total energy and maximum gradient steadily decreased during the minimization).
 ### 4.2 Heating, Equilibration and Production
 Next, we first run 200ps of heating that increase the temperature from 0 K to 300 K. Then the systems are equilibrated for 100 ns under constant pressure and temperature (NPT) to adjust the system density; followed by 100 ns production simulations under constant volume and temperature (NVT) conditions. Below are the input files of heating, equilibration, and production stages:




 [heat.in](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/MD/heat.in):
 ```
200ps heating with restraint on protein-PCL
&cntrl
    imin = 0,
    irest = 0,
    ntx = 1,
    ntb = 1,
    cut = 10.0,
    ntr = 1,
    ntc = 2,
    ntf = 2,
    tempi = 0.0,
    temp0 = 300.0,
    ntt = 3,
    gamma_ln = 1.0,
    nstlim = 100000, dt = 0.002,
    ntpr = 1000, ntwx = 1000, ntwr = 10000
/
Keep protein-PCL fixed with weak restraints
10.0
RES 1 161
END
END
```
 Note: Change the residue numbers to 1-497 for HGMS/ACP-Ppant-Ser.
 [equil.in](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/MD/equil.in):
 ```
100ns equilibration MD (NPT)
&cntrl
    imin = 0, irest = 1, ntx = 7,
    ntb = 2, pres0 = 1.0, ntp = 1,
    taup = 2.0,
    cut = 10.0, ntr = 0,
    ntc = 2, ntf = 2,
    tempi = 300.0, temp0 = 300.0,
    ntt = 3, gamma_ln = 1.0,
    nstlim = 50000000, dt = 0.002,
    ntpr = 100000, ntwx = 100000, ntwr = 1000000
/
```
 [prod.in](https://github.com/ShijiZ/PFF-Library/blob/main/tutorial/include/MD/prod.in):
 ```
100ns Production MD (NVT)
&cntrl
    imin = 0, irest = 1, ntx = 7,
    ntb = 1, ntp = 0,
    cut = 10.0, ntr = 0,
    ntc = 2, ntf = 2,
    tempi = 300.0, temp0 = 300.0,
    ntt = 3, gamma_ln = 1.0,
    nstlim = 50000000, dt = 0.002,
    ntpr = 100000, ntwx = 100000, ntwr = 100000
/
```
 Run the following three commands in sequence to perform heating, equilibration, and production:
 ```
pmemd.cuda_SPFP -O -i heat.in -p *.prmtop -c min2.rst -ref min2.rst -o heat.out -r heat.rst -x heat.nc
pmemd.cuda_SPFP -O -i equil.in -p *.prmtop -c heat.rst -o equil.out -r equil.rst -x equil.nc
pmemd.cuda_SPFP -O -i prod.in -p *.prmtop -c equil.rst -o prod.out -r prod.rst -x prod.nc
```
 Note we use *pmemd.cuda_SPFP* to run our simulation, which finishes much faster than the CPU program *pmemd*. And as always, check the resulting structure and trajectory with your favorite visualization tool to make sure nothing obviously bad happened.


