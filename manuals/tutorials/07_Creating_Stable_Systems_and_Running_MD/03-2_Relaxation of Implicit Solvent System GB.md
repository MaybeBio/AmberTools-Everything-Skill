

# Relaxation of Implicit Solvent Systems (GB)
 This tutorial has inputs for relaxing an implicit solvent system very similar to what was outlined in: 

 Nguyen, H., Maier, J., Huang, H., Perrone, V. and C. Simmerling, Folding Simulations for Proteins with Diverse Topologies are Accessible in Days with Physics-Based Force Field and Implicit Solvent, [*J. Am. Chem. Soc.*, **2014,** *136(40),* 13959–13962.](https://pubs.acs.org/doi/full/10.1021/ja5032776)
 ## Learning Outcomes
 Be able to relax a system under implicit solvent conditions with a GB model.
 ## Introduction
 An implicit solvent model mimics the presence a solvent in an average manner. In contrast, an explicit solvent model literally represents the solvent with lots (usually thousands) of individual atomistic solvent molecules surrounding the solute. Common explicit solvent models employed in Amber are found in the [Amber Reference Manual](https://ambermd.org/doc12/Amber23.pdf#page=52). Videos by Chris Cramer outline the difference between [explicit solvent](https://www.youtube.com/watch?v=caxpxBT8JsM&list=PLkNVwyLvX_TFBLHCvApmvafqqQUHb6JwF&index=34) and [implicit solvent](https://www.youtube.com/watch?v=VN6D5lgJvUc&list=PLkNVwyLvX_TFBLHCvApmvafqqQUHb6JwF&index=37) environments in general.
 In Amber, the most computationally efficient method for modeling solvent is the Generalized Born (GB) method. For a review of Generalized Born (GB) methods, please see:
- Onufriev, A.V. and D.A. Case "Generalized Born Implicit Solvent Models for Biomolecules" [*Annu. Rev. Biophys.* **2019,** *58,* 275-296](https://www.annualreviews.org/doi/10.1146/annurev-biophys-052118-115325).
  Implicit solvent drastically accelerates molecular dynamics simulations by reducing the number of particles in the system. In Amber, this GB model is the most commonly used implicit solvent model. It has been widely tested on the force fields ff99SB and ff14Bonlysc. In this tutorial, minimization, heating, and equilibration steps are shown.
 ## Process
 In this tutorial, you will: - [Prepare a pdb file](#pdb)
- [Minimize bad contacts](#min)
- [Heat the system](#heat)
- [Relax the system (2 different atom restraint)](#relax)
 ## Prepare a pdb file
  The tc5b variant of the trp-cage mini-protein will be used in this tutorial. The pdb file for the protein can be found [here at the PDB databank (PDBID:1L2Y)](https://www.rcsb.org/structure/1L2Y).
This pdb file should be converted into amber-happy format and prepared with hydrogen mass repartitioning by LEaP.
To set up your pdb file:
 ### 1. Set AMBERHOME. For more information see [this](http://archive.ambermd.org/201006/0689.html).
 ```
setenv AMBERHOME=/path/amber/

```
 or ```
export AMBERHOME="/path/amber/"

```
 ### 2. Use [pdb4amber](https://ambermd.org/tutorials/basic/tutorial9/index.php)
 ```
pdb4amber 1l2y.pdb>1l2y.amber.pdb
```
 ### 3. Use leap to set radii and make a topology file and starting coordinates.
 This uses the [ff14SBonlysc](https://ambermd.org/AmberModels_proteins.php) force field. ```
tleap -f tleap.in > tleap.out
```
 **tleap.in file**
 ```
source leaprc.protein.ff14SBonlysc
tc5b = loadpdb 1l2y.amber.pdb
set default pbradii mbondi3
saveamberparm tc5b tc5b.1l2y.parm7 tc5b.1l2y.rst7

```
 ### 4. [Perform hmass repartitioning with parmed.](https://ambermd.org/tutorials/basic/tutorial9/index.php)
 ```
parmed tc5b.1l2y.parm7

```
 parmed will write some info. ```
hmassrepartition

```
 parmed will write out "Repartitioning hydrogen masses to 3.024 daltons. Not changing water hydrogen masses." ```
outparm tc5b.1l2y.hmass.parm7
quit

```
 After those steps, you should now have: - Topology file after hydrogen repartitioning: [tc5b.1l2y.hmass.parm7](https://ambermd.org/tutorials/basic/tutorial15/include/tc5b.1l2y.hmass.parm7)
- Coordinate file: [tc5b.1l2y.rst7](https://ambermd.org/tutorials/basic/tutorial15/include/tc5b.1l2y.rst7)
 ## Run Minimization - Eliminate Bad Contacts
 Energy minimization step is necessary because bad contacts created by addition of water, hydrogens, and other atoms need to be removed. Without this minimization step, the energy in that bad contact region can be too high enough to crash the simulation. Minimization moves the atoms based on the forces and finds the closest structure at its energy minimum.
 ### 1. Create an input for minimization "min.in":
  ```
energy minimization
 &cntrl
  imin = 1, 
  maxcyc=100,  
  ntx = 1, 
  ntwr = 100, 
  ntpr = 10, 
  ioutfm=0, 
  ntxo=1, 
  cut = 1000.0, 
  ntb=0, 
  igb = 8, 
  gbsa=3, 
  surften=0.007, 
  saltcon = 0.0, 
 &end

```
 **Settings**
 ```
imin = 1      :Choose a minimization run
maxcyc=100    :Maximum minimization cycles
ntx = 1,      :Read coordinates with no velocities from the restart file
ntwr = 100    :Write "restrt" file in every 100 steps
ntpr = 10     :Print to the Amber "mdout" output file every 10 cycles
ioutfm=0      :Specify format of coordinate file as binary NetCDF file
ntxo=1        :Specify the format of the final coordinates, 
                 velocities in the restart file to be ASCII format
cut = 1000.0  :Nonbounded cutoff distance in angstroms. Higher the 
                 number, better accuracy with increased computational cost
ntb=0         :No periodic boundaries imposed on the system during the 
                 non-bounded interactions calculation 
igb = 8       :Specify GBneck2, the implicit solvent model
gbsa=3        :Compute surface area using a fast pairwise approximation 
                 suitable for GPU machine 
surften=0.007 :Set the surface tension used to calculate the nonpolar 
                 contribution to the free energy of solvation 
saltcon = 0.0 :Set the concentration (M) of 1-1 mobile counterions in 
                 the solution

```
 ### 2. Run the minimization.
  Here is an example script for running minimization with pmemd.cuda
 ```
#!/bin/bash
export AMBERHOME=/home/packages/amber20
/home/packages/amber20/bin/pmemd.cuda -O -i min.in\
 -p tc5b.1l2y.hmass.parm7 -c tc5b.1l2y.rst7 -o min.out\
 -x min.crd -inf min.info -r min.rst7\


```
 ### 3. Check the output file for minimization.
 To know if this process ran correctly, it is important to check the output information. The output (min.out) for this minimization process should look like this:
 ```
          -------------------------------------------------------
          Amber 20 PMEMD                              2020
          -------------------------------------------------------

| PMEMD implementation of SANDER, Release 18

|  Compiled date/time: Mon May  4 04:10:34 2020
| Run on 07/30/2021 at 23:04:25

|   Executable path: /home/packages/amber20/bin/pmemd.cuda
| Working directory: /home/jchong
|          Hostname: Unknown
  [-O]verwriting output

File Assignments:
|   MDIN: min.in                                                                
|  MDOUT: min.out                                                               
| INPCRD: tc5b.1l2y.rst7                                                        
|   PARM: tc5b.1l2y.hmass.parm7                                                 
| RESTRT: min.rst7                                                              
|   REFC: refc                                                                  
|  MDVEL: mdvel                                                                 
|   MDEN: mden                                                                  
|  MDCRD: min.crd                                                               
| MDINFO: min.info                                                              
|  MDFRC: mdfrc                                                                 


 Here is the input file:

energy minimization                                                            
 &cntrl                                                                        
  imin = 1, maxcyc=100,                                                        
                                                                               
  ntx = 1,                                                                     
                                                                               
  ntwr = 100, ntpr = 10, ioutfm=0, ntxo=1,                                     
                                                                               
  cut = 1000.0,                                                                
                                                                               
  ntb=0, igb = 8, gbsa=3, surften=0.007, saltcon = 0.0,                        
                                                                               
 &end                             

```
 (..and more output until you see this info about your simulated system...)
 ```
----------------------------------------------------------------------------
   1.  RESOURCE   USE: 
----------------------------------------------------------------------------

 NATOM  =     304 NTYPES =      12 NBONH =     150 MBONA  =     160
 NTHETH =     346 MTHETA =     219 NPHIH =     700 MPHIA  =     656
 NHPARM =       0 NPARM  =       0 NNB   =    1701 NRES   =      20
 NBONA  =     160 NTHETA =     219 NPHIA =     656 NUMBND =      53
 NUMANG =     124 NPTRA  =     136 NATYP =      26 NPHB   =       0
 IFBOX  =       0 NMXRS  =      24 IFCAP =       0 NEXTRA =       0
 NCOPY  =       0

 Implicit solvent radii are ArgH and AspGluO modified Bondi2 radii (mbondi3)                                
 Replacing prmtop screening parameters with GBn2 (igb=8) values

```
 (...then you will see a list of all of the settings, including ones you did not specify in your input--the defaults are used...)
 Result section should look like this:
 ```
----------------------------------------------------------------------------
   4.  RESULTS
----------------------------------------------------------------------------



   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
      1      -5.7351E+02     4.3317E+00     1.6622E+01     CG        165

 BOND    =       12.3844  ANGLE   =       39.0680  DIHED      =     237.2271
 VDWAALS =     -119.8407  EEL     =    -1381.8316  EGB        =    -322.8611
 1-4 VDW =       62.7666  1-4 EEL =      862.3303  RESTRAINT  =       0.0000
 ESURF   =       37.2492


   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
     10      -5.8438E+02     8.0802E-01     5.0348E+00     C         302

 BOND    =       10.3356  ANGLE   =       39.0395  DIHED      =     235.8701
 VDWAALS =     -122.4424  EEL     =    -1377.1130  EGB        =    -328.7565
 1-4 VDW =       61.6567  1-4 EEL =      859.7844  RESTRAINT  =       0.0000
 ESURF   =       37.2492

```
 This goes on upto NSTEP of 100.
 ## Heating
 In this step, it sets the temperature to run in the simulation. Notice that it is not heated directly to the desired temperature, but rather heating slowly from the lower temperature. Along with the heating, this job script also inputs the thermostat options and water restraint. In this tutorial, the temperature was controlled with a Langevin thermostat with collision frequency gamma=1.0ps-1, and the desired temperature is 300K.
 ### 1. Create an input for heating (heat.in):
 ```
heating with backbone restraints
 &cntrl
   ntx = 1,                             
   ntwx = 5000, 
   ntwe = 0, 
   ntwr = 500,     
   ntpr = 5000,                 
   ioutfm=0, 
   ntxo=1,                      
   imin = 0,
   nstlim = 500000, 
   dt = 0.002, 
   ntt = 3,
   gamma_ln=1.,
   temp0 = 100.0,   
   nscm = 1000,
   ig = -1,                  
   ntc= 2, 
   ntf = 2, 
   cut = 1000,           
   igb=8,
   gbsa=3,
   surften=0.007,
   ntb = 0, 
   saltcon=0., 
   ntr=1,
   restraintmask="@CA,N,C,O",     
   restraint_wt=10.0,  
   nmropt=1,                              
 &end
 &wt
   TYPE="TEMP0", istep1=0, istep2=500000, 
   value1=100., value2=300.,              
 &end
 &wt
   TYPE="END",
 &end

```
 **Settings**
 ```
ntx = 1         :Read coordinates with no velocities from the restart file
ntwx = 5000     :Write Amber trajectory file mdcrd every 5000 steps
ntwe = 0        :Write energies and temperatures in"mden" file in every 
                    "ntwe" steps. 0 means no mden file will be written.
ntwr = 500      :Write "restrt" file every 500 steps. 
ntpr = 5000     :Print to the Amber "mdout" output file every 5000 cycles.           
ioutfm=0        :Specify format of coordinate file to be binary NetCDF file
ntxo=1          :Specify the format of the final coordinates, velocities 
                       in the restart file to be ASCII format
imin = 0        :Choose a MD run (no minimization)
nstlim = 500000 :Number of MD steps in run
dt = 0.002      :Time length of each step in picoseconds
ntt = 3         :Temperature control with Langevin thermostat
gamma_ln=1.     :Langevin thermostat collision frequency. 1.0ps-1
temp0 = 100.0   :Set the first thermostat temperature
nscm = 1000     :Remove of transitional and rotational center of mass motion 
                       at regular intervals.
ig = -1         :Randomize the seed for the pseudo-random number generator
ntc= 2          :Turn on the SHAKE algorithm for hydrogen atoms to perform 
                       bond length constraints. 
ntf = 2         :Force evaluation. Usually ntc=ntf, and by ntf=2, bond 
                       interactions involving hydrogen atoms are omitted; 
                       “cut” determines the bond cutoff in angstrom
cut = 1000      :Nonbounded cutoff distance in angstroms
igb=8           :Specify GBneck2, the implicit solvent model
gbsa=3          :Compute surface area using a fast pairwise approximation 
                        suitable for GPU machine. 
surften=0.007   :Set the surface tension used to calculate the nonpolar 
                        contribution to the free energy of solvation. 
ntb = 0         :No periodic boundaries imposed on the system during the 
                        non-bounded interaction calculation. 
saltcon=0.      :Set the concentration (M) of 1-1 mobile counterions in 
                        the solution. 
ntr=1           :Restrain specified atoms
restraintmask="@CA,N,C,O"    :A list of atoms for restraining.     
restraint_wt=10.0   :The weight of restraint, the spring force constant, 
                             is set to be 10 kcal/mol/angstrom2. 
nmropt=1        :Read inputs between "&wt"
istep1=0        :The first step at which the temperature raise starts
istep2=500000   :The last step at which the temperature raise ends
value1=100      :Initial temperature
value2=300      :Target, final temperature 

```
 Here is a sample submission script.
 ```
#!/bin/bash
export CUDA_VISIBLE_DEVICES=3
/home/packages/amber21/bin/pmemd.cuda -O -i heat.in\
 -p tc5b.1l2y.hmass.parm7 -c min.rst7 -ref min.rst7\
 -o heat.out -x heat.crd -inf heat.info -r heat.rst7\


```
 ### 2. Check the output file for the heating process.
 By checking your output file, you are able to see if the temperature raised from 100K to 300K overtime. Also, check if the temperature has fluctuations over time while raising to 300K, not going straight up to 300K. The output (heat.out) for this heating process should look like this:
```
      -------------------------------------------------------
          Amber 20 PMEMD                              2020
          -------------------------------------------------------

| PMEMD implementation of SANDER, Release 18

|  Compiled date/time: Tue May  4 00:51:49 2021
| Run on 08/01/2021 at 21:45:12

|   Executable path: /home/packages/amber21/bin/pmemd.cuda
| Working directory: /home/jchong/Amber_2021summer
|          Hostname: Unknown
  [-O]verwriting output


```
 (and so on) ## Relaxation steps with different atomic restraints
 In the MD simulation, atoms in the system have to undergo a relaxation process to reach to a stationary state in the silico environment. The initial nonstationary segment of the simulated trajectory are typically discarded in the calculation of equilibrium properties in this process.
 Keeping the temperature from the heating step constant, MD is now run with different atom restraints. In this case, we want to create 2 stages with different restraints.
 ### 1. Create an input for the first relaxation process "eq1.in"
 ```
 &cntrl
   ntx = 5,                              
 
   ntwx = 5000, ntwe = 0, ntwr = 500,     
   ntpr = 5000,                 
 
   ioutfm=0, ntxo=1, 
                     
   imin = 0, nstlim = 125000, dt = 0.002, 

   ntt = 3, gamma_ln=1.0, temp0 = 300.0,   
   nscm = 1000, ig = -1, irest = 1,                  
 
   ntc= 2, ntf = 2, cut = 1000,           
 
   igb=8, gbsa=3, surften=0.007, ntb = 0, 
   saltcon=0., 
 
   ntr=1, restraintmask="@CA,N,C,O",      
   restraint_wt=1.0,  
 /

```
 All flags are almost identical to the heating run. However, since this run is continuing from the heating as an equilibration run, the velocity information from the restart file is used. "ntx" and "irest" are the flags to tell Amber to start equilibration from the restart file.
 **Additional Settings**
 ```
ntx = 5   :Reading coordinates and velocities from the previous run.
irest = 1 :Restarting the simulation, reading coordinates and 
                velocities from a previously saved restart file.

```
 The weight for the positional restraints of the specified atoms (CA, N, C, and O’’) are now 1.0 kcal/mol/angstrom2 instead of 10. This stage is also 250 ps by taking 62500 steps with 0.004 ps for each.
 ### 2. Create a second input for the relaxation process "eq2.in"
 Now we want to decrease the restraints for the atoms from the first equilibration. This stage is a 250 ps stage, and the temperature is kept constant after the heating. “Restraint_wt” is now set to be 0.1 kcal/mol/angstrom2.
```
 &cntrl
   ntx = 5,                              
   ntwx = 5000, ntwe = 0, ntwr = 500,     
   ntpr = 5000,                 

   ioutfm=0, ntxo=1,                     
 
   imin = 0, nstlim = 125000, dt = 0.002, 
 
   ntt = 3, gamma_ln=1., temp0 = 300.0,   
   nscm = 1000, ig = -1, irest = 1,                 
 
   ntc= 2, ntf = 2, cut = 1000,           
 
   igb=8, gbsa=3, surften=0.007, ntb = 0, 
   saltcon=0., 
 
   ntr=1, restraintmask="@CA,N,C,O",      
   restraint_wt=0.1,  
 /

```
 Here is a sample job script for the relaxation steps.
 ```
#!/bin/bash
export AMBERHOME=/home/packages/amber20/
export CUDA_VISIBLE_DEVICES=3
/home/packages/amber20/bin/pmemd.cuda -O -i eq1.in\
 -p tc5b.1l2y.hmass.parm7 -c heat.rst7 -ref heat.rst7\
 -o eq1.out -x eq1.crd -inf eq1.info -r eq1.rst7\

/home/packages/amber20/bin/pmemd.cuda -O -i eq2.in\
 -p tc5b.1l2y.hmass.parm7 -c eq1.rst7 -ref eq1.rst7\
 -o eq2.out -x eq2.crd -inf eq2.info -r eq2.rst7\


```
 After the system is equilibrated and the outputs are checked, you should visualize the coordinates of your equilibrated structure to make sure it looks OK, the system didn't explode and your peptide looks reasonable. Then you are ready to go onto a production run with implicit solvent.


