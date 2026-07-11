

# Running MD in parallel – multisander and multipmemd
 ## Learning Outcomes
 - Be able to run multiple independent simulations in parallel using pmemd.MPI on central processing units (CPUs) or pmemd.cuda.MPI on graphical processing units (GPUs)
 ## Prerequisites
 - Understanding of the basics of molecular dynamics simulations and empirical forcefields
- Access to AMBER executables in a Linux environment built with MPI support - with optional SLURM queueing system - [the Amber Reference Manual](https://ambermd.org/Manuals.php) for installation instructions
- Familiarity with using a text editor - such *vim, ecmacs, nano,* or similar - to create and modify plain-text files - see [these tutorials](https://ambermd.org/tutorials/Tools.php)
- System prepared from [Building Protein Systems in Explicit Solvent](https://ambermd.org/tutorials/basic/tutorial7/index.php) to run the simulation - this tutorial only covers the production run, and does not go over minimization, heating, or equilibration.
 ## Introduction
 This tutorial will introduce you to running multiple independent molecular dynamics (MD) simulations in parallel using a **Groupfile** with the **`pmemd.MPI`** or **`pmemd.cuda.MPI`** executables. Groupfile input is required for many advanced simulation approaches in Amber, such as multiple (non-interacting) copies of a system to improve statistics for a simulation, Replica Exchange molecular dynamics, and nudged elastic band. See the “multisander (and multipmemd)” section of the [Amber Reference Manual](https://ambermd.org/Manuals.php) (21.12. multisander (and multipmemd)) for additional information. You can use the protein system prepared in [Building Protein Systems in Explicit Solvent](https://ambermd.org/tutorials/basic/tutorial7/index.php) to run the simulation. Note that this tutorial only covers the production run, and does not go over minimization, heating, or equilibration.
 ## Running on GPUs with pmemd (previous tutorial)
 Particle Mesh Ewald Molecular Dynamics (pmemd) is the primary engine for running MD simulations with Amber. The `sander` executable only runs on cpus but it will be much slower than `pmemd`. The `pmemd.CUDA` executable requires [GPUs](https://ambermd.org/GPUSupport.php) that greatly increase the simulation performance (typically reported in terms of ns/day). The functionality of `pmemd` supports Particle Mesh Ewald simulations, Generalized Born simulations, Isotropic Periodic Sums, ALPB solvent, and gas-phase simulations. Note that `pmemd` is not a complete implementation of `sander`. See the pmemd section of the [Amber Reference Manual](https://ambermd.org/Manuals.php) for more information.
 ## Initial configurations for independent runs
 Starting from the [relaxed system](https://ambermd.org/tutorials/BuildingSystems.php), we can set up multiple independent runs to improve the statistics of production-run observables. Each of these runs should be both distinct, with respect to each other, and representative of the thermodynamic state variables as achieved during relaxation. There are various ways to setup the initial coordinates for each independent run. For example, you could repeat the build and relaxation for the desired number of runs using the random placement of solvent and ions in [tleap](https://ambermd.org/tutorials/BuildingSystems.php). Note that for some advanced techniques, a single shared topology file is required. For implicit solvent simulations, the user can choose four different coordinate frames from relaxation of the initial structure ( [Relaxation of Implicit Solvent System (GB)](https://ambermd.org/tutorials/basic/tutorial15/index.php)). Before proceeding to next step, please download the files required for this tutorial ( [here](https://ambermd.org/tutorials/basic/tutorial21/include/pmemd/)).
 In this tutorial, we will use the RAMP1.prmtop topology file and four different starting coordinates of taken from the system relaxation step of [Building Protein Systems in Explicit Solvent](https://ambermd.org/tutorials/basic/tutorial7/index.php):
 - RAMP1_equil.01.rst7
- RAMP1_equil.02.rst7
- RAMP1_equil.03.rst7
- RAMP1_equil.04.rst7
 Different starting coordinates can include different conformations of the same molecule (that is, the same topology will be used to define the system). In the most typical case, differentiating systems is achieved after building in tleap, where ions are randomized N number of times, with N=the number of copies the user ultimately wants to create. Then, each copy is individually minimized and equilibrated.
 ## Process
 ### 1. Preparing the input file
 

 md.in ```
Explicit solvent MD, constant pressure, 10 ns MD, print every 10 ps
 &cntrl
   imin=0, irest=1, ntx=5, 
   ntpr=5000, ntwx=5000, ntwr=5000, nstlim=5000000, 
   dt=0.002, ntt=3, tempi=300, 
   temp0=300, gamma_ln=1.0, ig=-1, 
   ntp=1, ntc=2, ntf=2, cut=9, 
   ntb=2, iwrap=1, ioutfm=1, 
/ 

```
 These settings are appropriate for a short simulation (10 ns) of a protein system in a truncated octahedral box solvated in explicit water, but you may need to change some aspects of it for your own system.
 ### 2. Preparing the Groupfile
 The ‘groupfile’ combines the command line instructions for running pmemd for each copy into a single file.
 Previously, in the ‘Running MD with pmemd’ tutorial, you used this line to submit pmemd jobs to a GPU:
 ```
$AMBERHOME/bin/pmemd.cuda -O -i md.mdin -p RAMP1.prmtop -c RAMP1_eq7.rst7\
 -ref RAMP1_eq7.rst7 -o RAMP1_md.mdout -r RAMP1_md.rst7 -x RAMP1_md.nc

```
 Now, we will need to bundle everything after the executable into a file, and create a newline entry containing the instructions for pmemd for each of the four independent runs. We will have four copies, so four simulations, and our groupfile will look like the following:
 ```
-O -i md.in -p RAMP1.prmtop -c RAMP1_equil.01.rst7 -ref RAMP1_equil.01.rst7 -x md.nc.01\
 -r 01.rst7 -o md.out.01 -inf md.info.01 -l logfile.01
-O -i md.in -p RAMP1.prmtop -c RAMP1_equil.02.rst7 -ref RAMP1_equil.02.rst7 -x md.nc.02\
 -r 02.rst7 -o md.out.02 -inf md.info.02 -l logfile.02
-O -i md.in -p RAMP1.prmtop -c RAMP1_equil.03.rst7 -ref RAMP1_equil.03.rst7 -x md.nc.03\
 -r 03.rst7 -o md.out.03 -inf md.info.03 -l logfile.03
-O -i md.in -p RAMP1.prmtop -c RAMP1_equil.04.rst7 -ref RAMP1_equil.04.rst7 -x md.nc.04\
 -r 04.rst7 -o md.out.04 -inf md.info.04 -l logfile.04

```
 *Note that all the above commands are on a single line per copy, so four lines total for the four copies you will run.
 ### 3. Running the multipmemd MD
 Now we have all of the necessary files to run the multipmemd MD simulation. On a designated machine or interactive session (on HPC resource with queueing system), you can do a test run by editing the md.in to have a smaller number of steps (nstlim=5000) BEFORE running pmemd.cuda.MPI with mprirun using the following command:
 ```
mpirun -np 4 $AMBERHOME/bin/pmemd.cuda.MPI -O -ng 4 -groupfile groupfile

```
 Running 10 ns of this system takes four or more hours on cpus. However, it is more likely that you will be using an HPC resource with a queueing system, such as SLURM or PBS. While proper usage of HPC and queueing systems is outside the scope of this tutorial, we included an example slurm script for running multipmemd on GPUs that you can adapt to your HPC resources:
 ```
#!/bin/bash
#SBATCH -t 46:00:00         #wallclock time requested
#SBATCH --gres=gpu:4        #request for 4 GPUs
#SBATCH --cpus-per-gpu=1    #request 1 CPU per GPU
#SBATCH --mem=4G		    #set compute node memory to 4G
#SBATCH --nodelist=md-gpu-5 #request this to run on my md node
#SBATCH -J testjob	    #jobname

#Set location of Amber
source ~/amber22.env.sh     #source the environment script to set $AMBERHOME
SANDER=pmemd.cuda.MPI	    #set pmemd.cuda.MPI to a variable named SANDER

cd $SLURM_SUBMIT_DIR

# Run executable specified above ($SANDER) using mpi (mpirun), 4 copies (-ng) 
# each on its own processor (-np 4), using the groupfile named ‘groupfile’ 
mpirun -np 4 $SANDER -O -ng 4 -groupfile groupfile

exit 0			    #exit

```
 Submit job to SLURM:
 ```
sbatch qsub.sh

```
 Wait for job to finish. It should take around 4 h to 7 h (see mdinfo), depending on the GPU you are using.
 All of the files used and generated in this tutorial are available [here](https://ambermd.org/include/pmemd/):
 - RAMP1.prmtop
- RAMP1_equil.01.rst7, RAMP1_equil.02.rst7, RAMP1_equil.03.rst7, RAMP1_equil.04.rst7
- md.in
- groupfile
- qsub.sh
- md.out.01, md.out.02, md.out.03, md.out.04
- 01.rst7, 02.rst7, 03.rst7, 04.rst7
- md.nc.01, md.nc.02, md.nc.03, md.nc.04
 ### 4. Analysis
 The next steps would be to analyze the trajectory. Tutorials on using cpptraj to perform analysis can be found at in [Section 4 Trajectory Analysis](https://ambermd.org/tutorials/TrajectoryAnalysis.php). It is also sensible to visualize the trajectory. [VMD](https://ambermd.org/tutorials/VMD.php) and [Chimera](https://ambermd.org/tutorials/Chimera.php) are common visualization software.
 This tutorial was written by Christina Bergonzo, Akanksha Manghrani, Abigail Held and Maria Nagan. 

 Last updated on June 17, 2025.
 ## Disclaimer
 This article was funded by the National Institute of Standards and Technology and is not subject to U.S. Copyright. Certain commercial equipment, instruments, materials or software are identified in this paper to foster understanding. Such identification does not imply recommendation or endorsement by the National Institute of Standards and Technology, nor does it imply that the materials or equipment identified are necessarily the best available for the purpose.


