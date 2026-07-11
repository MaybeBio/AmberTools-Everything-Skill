

# Temperature Replica Exchange MD
 ## Learning Outcomes
 - Be able to set up and run temperature replica exchange in parallel on GPUs
 ## Introduction
 This tutorial is a basic guideline for running Temperature Replica Exchange Molecular Dynamics (T-REMD) simulations with the current implementation of Amber.
 **Assumptions:** This tutorial assumes that you have experience building systems using tleap, running MD simulations in parallel with pmemd.MPI or pmemd.cuda.MPI, and know about the terminology and typical variables used for standard MD simulations. If you are not familiar with setting up a system, equilibration and running MD simulations please refer to the more basic tutorials before attempting to run T-REMD simulations.
 You should also be familiar with the REMD theory. Make sure to read the paper by [Sugita and Okamoto ("Replica-exchange molecular dynamics method for protein folding." Chemical Physics Letters 314(1-2): 141-151, 1999)](https://ambermd.org/tutorials/advanced/tutorial7/files/remd_chem_phys_lett.pdf) and the relevant sections in the [AMBER manual](https://ambermd.org/doc12/Amber21.pdf#page=478) to learn about REMD.
 Note: T-REMD with explicit solvent can be run with constant pressure or volume (V) conditions. However, if you choose to run at higher temperatures (≥400 K) you probably should run with constant V so your sample doesn't boil. In addition, any type of REMD in Amber requires access to multiple processors (one per replica) simultaneously available.
 All the files needed for this tutorial can be found [here](https://ambermd.org/tutorials/advanced/tutorial42/include/ReplicaExchange).
 ## Process
 ### 1. Setting up the system: Building and Equilibration
 In this section the user should build the system using tLeap, and run minimization using sander. For this tutorial we will use:
 - the peptide Ala10 ( `m = sequence { ACE ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA NHE }`)
- the ff19SB force field ( `source leaprc.protein.ff19SB`)
- the mbondi2 GB radii ( `set default PBRadii mbondi2`)
- the `igb=5` option in sander.
 Make sure to check the tleap output and log file for any errors. For consistency in this tutorial, please build:
- a parmtop named `ala10.prmtop`
- input coordinates named `ala10.inpcrd`
- a pdb file named `ala10.pdb`.
  Next, run minimization to relax the input structure.
 ### 2. Determining the number of replicas and temperatures for each replica
 Before starting the replica exchange part of the simulation, we need to determine the number of replicas needed and the temperatures to which to set each replica.
To obtain reasonable exchange probabilities (20 % to 30 % acceptance rates) we need to have overlap between the potential energy distributions of neighboring replicas. Closer temperature spacing would increase this overlap and result in higher probabilities of exchange acceptance, but require more replicas to cover the desired temperature range.
 Finding an optimum temperature distribution can be tricky as there are many variables. Some things to think about include:
- the size of your system and if you can reduce it,
- if your system behaves the way you need it to in implicit solvent vs. requiring explicit solvent
- the temperature range you will need to effectively study the region of phase space you are interested in. - For example, is sampling at 277 K necessary, or will 300 K suffice? How about 325 K?
- Are we above the Tm of the protein?
- How high do we need to increase T to sample important conformations?
 You are advised to read the literature for examples.
 **Calculating Temperature Spacing** Temperature spacing is directly calculable from and proportional to the number of degrees of freedom in one’s system. The webserver here: [https://virtualchemistry.org/remd-temperature-generator/](https://virtualchemistry.org/remd-temperature-generator/) and python code on github here: [https://github.com/dspoel/remd-temperature-generator](https://github.com/dspoel/remd-temperature-generator) are available to help you select initial spacing [“Alexandra Patriksson and David van der Spoel, A temperature predictor for parallel tempering simulations Phys. Chem. Chem. Phys., 10 pp. 2073-2077 (2008) [https://dx.doi.org/10.1039/b716554d”](https://dx.doi.org/10.1039/b716554d)]. We will use a set of temperatures calculated using the above server.
 To determine the number of replicas and temperature distribution, one needs to determine the number of atoms in the system. This information can be found in both the topology and coordinate files. The first number in an ASCII formatted coordinate (or restart) file (ala10.inpcrd) shows the number of atoms.
 ```
ACE
   109
   2.0000010   1.0000000  -0.0000013   2.0000010   2.0900000   0.0000001
   1.4862640   2.4538490   0.8898240   1.4862590   2.4538520  -0.8898200
   3.4274200   2.6407950  -0.0000030   4.3905800   1.8774060  -0.0000066

```
 Alternatively, the number of atoms can be found in the topology file (ala10.prmtop):
 ```
%VERSION  VERSION_STAMP = V0001.000  DATE = 06/10/25  13:24:30
%FLAG TITLE
%FORMAT(20a4)
ACE
%FLAG POINTERS
%FORMAT(10I8)
     109       7      55      53     119      73     229     211       0       0
     561      12      53      73     211       9      18      15       8       0

```
 For both cases the number of atoms, 109, is highlighted. The number of replicas and temperature distribution is related to the square root of the number of atoms in the system.
 Usually REMD simulations are run for a temperature range between 270 K and 600 K. Depending on your system, a different temperature range may be required. An even number of replicas is required since exchanges are always attempted pair-wise. For the purpose of this tutorial, we will use 8 replicas and we will use the following temperature distribution (in K): 270.00, 303.23, 339.62, 379.48, 423.14, 470.95, 523.41, 580.91.
 You should copy the list of temperatures into a new file called `temperatures.dat` which will be used by the setup scripts later on.
 temperatures.dat
```
270.00
303.23
339.62
379.48
423.14
470.95
523.41
580.91

```
  ### 3. Generating chirality restraints
 To enhance conformational sampling, T-REMD runs typically go to very high temperatures for proteins (500 K to 600 K). High temperatures can cause unwanted rotations around the peptide bond, leading to non-physical chiralities. To prevent this, we need to use chirality restraints on the backbone. We can generate the restraints using the "makeCHIR_RST" executable that is provided in the AMBERHOME bin directory and should be in your PATH. To do this with the ala10.pdb we have already built using tleap:
```
> makeCHIR_RST ala10.pdb ala10_chir.dat

```
  To understand what this file is doing, you can find more information about distance, angle, and dihedral restraints in the AMBER manual under the NMR refinement section.
 ### 4. Equilibration
 Before starting the T-REMD simulations one needs to run short simulations with each replica set at its individual temperatures. This step uses Multisander to perform multiple independent (non-exchanging) MD runs. This type of simulation will run efficiently across a parallel cluster if you have access to one. See the tutorial, Running MD with PMEMD – Multisander.
 Generating the input files for multisander Since we have 8 replicas with temperatures from 270.00 K to 580.91 K we need 8 input files. The script "setup_equilibrate_input.sh" will generate the input files (mdin). This script uses a template called “equilibrate.mdin” for the input parameters, and reads the temperatures from the “temperatures.dat” file and adjusts temp0 of each mdin file to match. The script also ensures that the seed for the random number generator is different for each structure. This is very important, especially when running Langevin simulations, because without “ig=-1” set, unintended correlations can occur during the MD (see: https://pubs.acs.org/doi/10.1021/ct800573m). This script also sets up the groupfile.
 (Note: when running multisander it is important that every file that sander writes to be defined with a unique name. If this is not done, files are overwritten AND performance is poor. The user must specify unique names for mdinfo, mdout, rst and mdcrd.)
 For this section we will run 200 ps (nstlim=100000, dt=0.02) of simulation for each replica, heating them to their respective temperatures using a Langevin thermostat (ntt=3, gamma_ln=1.0). To read the chirality restraints we use "nmropt = 1" to turn the restraints on, and "DISANG=ala10_chir.dat" to define the file to read the restraints from. Since all simulations will use the same restraints we do not need to duplicate this file. We start by making the template input file called “equilibrate.mdin,” shown below. Note the " temp0=XXXXX" and "ig=RANDOM_NUMBER", which will be automatically adjusted to the selected Temperatures and a random number by the setup_equilibrate_input.sh script.
 equilibrate.mdin
```
Equilibration
 &cntrl
   irest=0, ntx=1, 
   nstlim=100000, dt=0.002,
   ntt=3, gamma_ln=1.0,
   temp0=XXXXX, ig=RANDOM_NUMBER,
   ntc=2, ntf=2, nscm=1000,
   ntb=0, igb=5,
   cut=999.0, rgbmax=999.0,
   ntpr=500, ntwx=500, ntwr=100000,
   nmropt=1,
 /
 &wt TYPE='END'
 /
 DISANG=ala10_chir.dat

```
  Since we will run 8 independent simulations in one multisander run we need a "groupfile" to define the individual input and output files for each simulation. An explanation of the groupfile format can be found in the AMBER manual and in Tutorial Running MD with pmemd - Multisander. Essentially it is a set of n replica lines specifying what we would normally have included on the command line for a regular MD run. The "setup_equsetup_equilibrate_input.sh" script will produce this for us. You can copy the text below to a new file, and make it executable by using:
```
> chmod +x setup_equilibrate_input.sh

```
  setup_equilibrate_input.sh
```
#!/bin/bash -f

if [ -f equilibrate.groupfile ];
then
  rm equilibrate.groupfile
fi

nrep=`wc temperatures.dat | awk '{print $1}'`
echo $nrep
count=0
for TEMP in `cat temperatures.dat`
do
  let COUNT+=1
  REP=`printf "%03d" $COUNT`
  echo "TEMPERATURE: $TEMP K ==> FILE: equilibrate.mdin.$REP"
  sed "s/XXXXX/$TEMP/g" equilibrate.mdin > temp
  sed "s/RANDOM_NUMBER/$RANDOM/g" temp > equilibrate.mdin.$REP
  echo "-O -rem 0 -i equilibrate.mdin.$REP -o equilibrate.mdout.$REP -c min.rst\
 -r equilibrate.rst.$REP -x equilibrate.mdcrd.$REP -inf equilibrate.mdinfo.$REP\
 -p ala10.prmtop" >> equilibrate.groupfile

  rm -f temp
done
echo "#" >> equilibrate.groupfile

echo "N REPLICAS  = $nrep"
echo " Done."

```
  Once we have generated the template, equilibrate.mdin file and created the setup_equilibrate_input.sh script and made it executable, we can then run the setup_equilibrate_input.sh to produce the multiple input files and the groupfile:
```
> ./setup_equilibrate_input.sh

```
  This should produce the following output to the terminal window:
```
8
TEMPERATURE: 270.00 K ==> FILE: equilibrate.mdin.001
TEMPERATURE: 303.23 K ==> FILE: equilibrate.mdin.002
TEMPERATURE: 339.62 K ==> FILE: equilibrate.mdin.003
TEMPERATURE: 379.48 K ==> FILE: equilibrate.mdin.004
TEMPERATURE: 423.14 K ==> FILE: equilibrate.mdin.005
TEMPERATURE: 470.95 K ==> FILE: equilibrate.mdin.006
TEMPERATURE: 523.41 K ==> FILE: equilibrate.mdin.007
TEMPERATURE: 580.91 K ==> FILE: equilibrate.mdin.008
N REPLICAS  = 8
Done.

```
  Along with a groupfile name equilibrate.groupfile that should be as follows: equilibrate.groupfile
```
-O -rem 0 -i equilibrate.mdin.001 -o equilibrate.mdout.001 -c min.rst -r equilibrate.rst.001 -x equilibrate.mdcrd.001 -inf equilibrate.mdinfo.001 -p ala10.prmtop
-O -rem 0 -i equilibrate.mdin.002 -o equilibrate.mdout.002 -c min.rst -r equilibrate.rst.002 -x equilibrate.mdcrd.002 -inf equilibrate.mdinfo.002 -p ala10.prmtop
-O -rem 0 -i equilibrate.mdin.003 -o equilibrate.mdout.003 -c min.rst -r equilibrate.rst.003 -x equilibrate.mdcrd.003 -inf equilibrate.mdinfo.003 -p ala10.prmtop
-O -rem 0 -i equilibrate.mdin.004 -o equilibrate.mdout.004 -c min.rst -r equilibrate.rst.004 -x equilibrate.mdcrd.004 -inf equilibrate.mdinfo.004 -p ala10.prmtop
-O -rem 0 -i equilibrate.mdin.005 -o equilibrate.mdout.005 -c min.rst -r equilibrate.rst.005 -x equilibrate.mdcrd.005 -inf equilibrate.mdinfo.005 -p ala10.prmtop
-O -rem 0 -i equilibrate.mdin.006 -o equilibrate.mdout.006 -c min.rst -r equilibrate.rst.006 -x equilibrate.mdcrd.006 -inf equilibrate.mdinfo.006 -p ala10.prmtop
-O -rem 0 -i equilibrate.mdin.007 -o equilibrate.mdout.007 -c min.rst -r equilibrate.rst.007 -x equilibrate.mdcrd.007 -inf equilibrate.mdinfo.007 -p ala10.prmtop
-O -rem 0 -i equilibrate.mdin.008 -o equilibrate.mdout.008 -c min.rst -r equilibrate.rst.008 -x equilibrate.mdcrd.008 -inf equilibrate.mdinfo.008 -p ala10.prmtop

```
  The “-rem 0” in the equilibrate.groupfile specifies that we do not wish to do actual replica exchange. The script should also produce eight mdin files named “equilibrate.mdin.001 to equilibrate.mdin.008.” You should check, by opening the file, or by using the command line “grep temp equilibrate.mdin.*”, that each mdin file produced has the desired temperature and that each contains a different value for the random number seed (ig=$RANDOM).
 To run the equilibration, we will use either sander.MPI, pmemd.MPI, or pmemd.cuda.MPI. When specifying the number of CPUs for the mpirun command, it is essential to ensure that it is a multiple of the number of replicas (in this case, 8). While a single processor can be used for testing purposes, using 8 cores will significantly improve run efficiency. If additional CPUs are available, you can specify more than 8 (e.g., 16 CPUs will allocate 2 CPUs per replica). As a best practice, it is recommended to benchmark your job during setup to optimize the use of your HPC resources. The example command line used is provided below; however, you may need to adjust it according to your MPI implementation, queuing system, and cluster architecture.
```
mpirun -np 8 $AMBERHOME/bin/sander.MPI -ng 8 -groupfile equilibrate.groupfile
```
  Here the -ng 8 option specifies that there are 8 groups and the -groupfile specifies the name of the groupfile to be read and indicates to sander.MPI that we want to run in multisander mode. Depending on the speed of your system this could take anywhere from a couple of minutes to several hours. You can check the progress by opening a separate terminal and tailing one of the mdout files. When it is complete you should have 8 rst files, 8 mdcrd files and 8 mdout files. These files are available [here](https://ambermd.org/tutorials/advanced/tutorial42/include/ReplicaExchange/00EQUIL/).
 You should check each output file to make sure that it read the restraints correctly (i.e., the “RESTRAINT” energy entry is greater than zero and generally increases with replica temperature), ran without problems (exited cleanly), and had the correct target temperature. You should also visualize the mdcrd or restart files to make sure there are no problems. This can be done by combining the restart files into a trajectory, then looking at them in VMD or Chimera.
 Now we have equilibrated our initial structures and we can move on to running the actual REMD simulations.
 ### 5. REMD Simulations
 Running REMD is essentially the same process as running the multisander equilibration described above. During REMD each multisander job (replica) will periodically communicate and attempt an exchange. We need to use additional variables in the input files for each replica, and the groupfile needs to be changed, to run REMD.
 Input files for REMD 

 As with multisander we need an input file for each replica. The input files are the same except for the target temperature (temp0) and the random number seed (ig). There are also some other minor changes compared with non-REMD simulations. Here are the necessary changes and explanations:
- nstlim = 500: During REMD this variable determines the number of MD steps between each exchange attempt. In this case we set it to 500 steps which, with a 2 fs timestep, gives us 1 ps between exchange attempts.
- numexchg = 1000: This specifies the number of exchange attempts during the simulation.
  Note that unlike normal MD simulations, the REMD nstlim parameter governs the steps between exchanges; the total simulation length will be numexch * nstlim * dt, which yields 1 ns of REMD simulation for the above parameter settings.
 We need to generate the input files again. First we create a generic mdin file:


 remd.mdin
```
Equilibration
 &cntrl
   irest=0, ntx=1, 
   nstlim=500, dt=0.002,
   ntt=3, gamma_ln=1.0,
   temp0=XXXXX, ig=RANDOM_NUMBER,
   ntc=2, ntf=2, nscm=1000,
   ntb=0, igb=5,
   cut=999.0, rgbmax=999.0,
   ntpr=100, ntwx=1000, ntwr=100000,
   nmropt=1,
   numexchg=1000,
 /
 &wt TYPE='END'
 /
 DISANG=ala10_chir.dat

```
  We then have a slightly modified version of the setup script above to build the multiple mdin files and the groupfile.


 setup_remd_input.sh
```
#!/bin/bash -f

if [ -f remd.groupfile ];
then
  rm remd.groupfile
fi

nrep=`wc temperatures.dat | awk '{print $1}'`
echo $nrep
count=0
for TEMP in `cat temperatures.dat`
do
  let COUNT+=1
  REP=`printf "%03d" $COUNT`
  echo "TEMPERATURE: $TEMP K ==> FILE: remd.mdin.$REP"
  sed "s/XXXXX/$TEMP/g" remd.mdin > temp
  sed "s/RANDOM_NUMBER/$RANDOM/g" temp > remd.mdin.$REP
  echo "-O -rem 1 -remlog rem.log -i remd.mdin.$REP -o remd.mdout.$REP -c equilibrate.rst.$REP -r remd.rst.$REP -x remd.mdcrd.$REP -inf remd.mdinfo.$REP -p ala10.prmtop" >> remd.groupfile" 

  rm -f temp
done
echo "#" >> remd.groupfile

echo "N REPLICAS  = $nrep"
echo " Done."

```
  Note that we have made some changes to the groupfile generation part of this script:
- `-rem 1`: Run Temperature REMD
- `-remlog rem.log`: Save the log file into "rem.log" to keep track of each replica
- `-c equilibrate.rst.$REP`: We are using the restart files from the equilibration as the input coordinates for each replica.
  Running this script creates the necessary input files, including a groupfile named remd.groupfile and the input files for production replica exchange MD.
```
> mpirun -np 8 $AMBERHOME/bin/sander.MPI -ng 8 -groupfile remd.groupfile

```
  Depending on the speed of your system this could take anywhere from a couple of minutes to several hours using sander.MPI. You can speed this up by using pmemd.MPI if you have > 1 CPU per replica, or by using pmemd.cuda.MPI, where you can run all 8 replicas on 1 GPU or (best case scenario) use one GPU per replica (8 GPUs total). When complete you should have 8 rst files, 8 mdcrd files, 8 mdout files and a rem.log file. These files are available [here](https://ambermd.org/include/ReplicaExchange/00PROD).
 The "rem.log" has summary information about each replica and exchanges. The columns provide (in this order) the replica number (corresponding to the order of files listed in the groupfile), the velocity scaling factor if an exchange has occurred (-1 otherwise), the instantaneous system temperature (but note that the thermostat temperature is used for the exchange calculation), the potential energy, the thermostat bath temperature during the previous MD stage, the new bath temperature for the next MD stage, the overall running average success rate between the Temp0 temperature and the next higher temperature (0 to 10, but it can be greater than 1 at the beginning since we divide the # of successes by 0.5* the total attempts since each pair is attempted only every other exchange), and finally a last integer that is used only for reservoir REMD, and is -1 otherwise.
 The log file looks like:
```
# Replica Exchange log file
# numexchg is       1000
# REMD filenames:
#   remlog= rem.log
#   remtype= rem.type
# Rep#, Velocity Scaling, T, Eptot, Temp0, NewTemp0, Success rate (i,i+1), ResStruct#
# exchange        1
 1     -1.00      0.00     25.08    269.50    269.50      0.00      -1
 2     -1.00      0.00     16.21    300.00    300.00      0.00      -1
 3     -1.00      0.00     49.19    334.00    334.00      0.00      -1
 4      1.06      0.00     57.15    371.80    413.90      2.00      -1
 5      0.95      0.00     66.86    413.90    371.80      0.00      -1
 6     -1.00      0.00     60.58    460.70    460.70      0.00      -1
 7     -1.00      0.00     67.98    512.90    512.90      0.00      -1
 8     -1.00      0.00     93.67    570.90    570.90      0.00      -1
# exchange        2
 1     -1.00    244.55      0.95    269.50    269.50      0.00      -1
 2     -1.00    251.98     13.77    300.00    300.00      0.00      -1
 3      1.06    241.54     31.16    334.00    371.80      1.00      -1
 4      0.95    311.46     29.45    371.80    334.00      1.00      -1
 5      1.06    297.63     39.26    413.90    460.70      1.00      -1
 6      0.95    381.98     34.99    460.70    413.90      0.00      -1
 7     -1.00    430.62     54.54    512.90    512.90      0.00      -1
 8     -1.00    424.99     80.56    570.90    570.90      0.00      -1
# exchange        3
 1     -1.00    238.63     11.17    269.50    269.50      0.00      -1
 2     -1.00    291.79     18.81    300.00    300.00      0.00      -1
 3     -1.00    325.13     23.99    334.00    334.00      0.67      -1
 4     -1.00    371.31     41.80    371.80    371.80      0.67      -1
 5     -1.00    375.01     54.78    413.90    413.90      0.67      -1
 6     -1.00    364.77     64.48    460.70    460.70      0.00      -1
 7     -1.00    449.13     81.14    512.90    512.90      0.00      -1
 8     -1.00    599.45     92.80    570.90    570.90      0.00      -1

```
  At this point, the REMD simulation has run. If you wish to restart a REMD simulation, you will need to change the irest flag in the input files to 1, and the ntx flag to 5, to indicate you wish to restart a run and should look in the coordinate file for velocities. You will also need to update the groupfile to use the first simulation’s restart files as the input coordinates for the starting point of the next run. After running the number of exchanges you desire, you can proceed to analysis.
 This tutorial is written by Christina Bergonzo and Akanksha Manghrani. 

 Other versions by Dan Roe, Asim Okur and Carlos Simmerling.
 ## Disclaimer
 This article was funded by the National Institute of Standards and Technology and is not subject to U.S. Copyright. Certain commercial equipment, instruments, materials or software are identified in this paper to foster understanding. Such identification does not imply recommendation or endorsement by the National Institute of Standards and Technology, nor does it imply that the materials or equipment identified are necessarily the best available for the purpose.


