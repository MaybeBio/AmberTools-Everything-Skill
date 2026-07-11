

# Replica Exchange MD - Analyzing Results
 ## Learning Outcomes
 - Be able to analyze results from T-REMD simulations
 ## Introduction
  This tutorial will introduce you to analyzing a Temperature replica exchange molecular dynamics (T-REMD) simulation (that you have completed in the previous [Replica Exchange MD Tutorial](https://ambermd.org/tutorials/advanced/tutorial42/index.php)). The scripts needed for this tutorial can be found [here](https://ambermd.org/include/scripts).  The output for T-REMD follows the REPLICAS, not the temperatures; each replica writes to its own output files as they exchange through temperatures. This means that we need to carry out extra post-processing on the trajectory files in order to obtain data for ensembles sampled at a specific temperature. This processing is described below.
 mdout files


 In a REMD run, the mdout file will have additional data that describes the information being used for the exchange calculation. An example is shown below, where the normal sander energy output (specified by ntpr) is written, along with information about the replica number (REPNUM) and exchange attempt counter (EXCHANGE#), as well as the current “TEMP0” since, remember, replica trajectories follow structures and sample different temperatures.
```
NSTEP =      200   TIME(PS) =       0.400  TEMP(K) =   199.60  PRESS =     0.0
 Etot   =       -29.6551  EKtot   =        53.9446  EPtot      =       -83.5997
 BOND   =        10.7369  ANGLE   =        32.4393  DIHED      =        17.5192
 UB     =         0.0000  IMP     =         0.0000  CMAP       =         4.8406
 1-4 NB =        21.2545  1-4 EEL =       642.1998  VDWAALS    =       -29.8651
 EELEC  =      -706.5033  EGB     =       -76.2257  RESTRAINT  =         0.0041
 EAMBER (non-restraint)  =       -83.6038
 TEMP0  =       270.0000  REPNUM  =              1  EXCHANGE#  =              1
 ------------------------------------------------------------------------------

 NMR restraints: Bond =    0.000   Angle =     0.000   Torsion =     0.004
===============================================================================

```
  remlog


 During a replica exchange MD simulation, the remlog will be written (the “-remlog” flag in the groupfile specifies this file name). In this file one will find details about what each replica is doing during an exchange attempt. Details about the potential energy calculation governing exchange acceptance are printed, as well as the former and current temperatures and success rates of the swaps between replicas. The remlog file can also be used by Cpptraj to sort trajectories by coordinate index (see relevant section in the Amber manual).
 One part of the remlog to check is the last exchange to make sure the exchange rate is successful across all replicas, as shown below for exchange 1000: all replicas show a 0.26 to 0.37 success rate (except for one replica, which is 0, which is the current highest temperature replica and cannot exchange with anything above itself).
```
# Replica Exchange log file
# numexchg is       1000
# REMD filenames:
#   remlog= rem.log
#   remtype= rem.type
# Rep#, Velocity Scaling, T, Eptot, Temp0, NewTemp0, Success rate (i,i+1), ResStruct#
# exchange        1
 1     -1.0000      0.00         -56.51    270.00    270.00      0.00      -1
 2     -1.0000      0.00         -44.95    303.23    303.23      0.00      -1
 3     -1.0000      0.00         -25.39    339.62    339.62      0.00      -1
 4     -1.0000      0.00         -36.38    379.48    379.48      0.00      -1
 5     -1.0000      0.00         -16.57    423.14    423.14      0.00      -1
 6      1.0542      0.00          39.23    470.95    523.41      2.00      -1
 7      0.9486      0.00          11.82    523.41    470.95      0.00      -1
 8     -1.0000      0.00          50.25    580.91    580.91      0.00      -1
# exchange        2
 1     -1.0000    219.66         -73.18    270.00    270.00      0.00      -1
 2     -1.0000    285.79         -65.88    303.23    303.23      0.00      -1
 3      1.0571    291.78         -48.27    339.62    379.48      1.00      -1
 4      0.9460    288.13         -49.10    379.48    339.62      0.00      -1
 5      1.0550    309.40         -22.02    423.14    470.95      1.00      -1
 6     -1.0000    415.91         -10.16    523.41    523.41      0.00      -1
 7      0.9479    359.98         -13.93    470.95    423.14      1.00      -1
 8     -1.0000    507.72           7.14    580.91    580.91      0.00      -1
.
.
.
# exchange     1000
 1     -1.0000    533.87           0.13    523.41    523.41      0.35      -1
 2      1.0571    311.19         -43.34    339.62    379.48      0.33      -1
 3     -1.0000    359.74         -34.70    423.14    423.14      0.34      -1
 4     -1.0000    472.18          17.86    470.95    470.95      0.37      -1
 5     -1.0000    271.58         -49.41    303.23    303.23      0.27      -1
 6     -1.0000    683.46          41.76    580.91    580.91      0.00      -1
 7     -1.0000    270.45         -78.15    270.00    270.00      0.26      -1
 8      0.9460    418.55         -35.50    379.48    339.62      0.31      -1

```
  ## Process
 REMD will produce trajectories for each replica. Each trajectory represents the snapshots sampled in the continuous MD run, which samples different thermostat temperatures during the run. Furthermore, this trajectory has additional information about the temperature at each snapshot so that Cpptraj can reconstruct the data at individual temperatures. In the next section, we provide scripts and explanations for each of these processes: extracting temperature trajectories from the entire set of REMD trajectories, calculating end-to-end distances for replica trajectories and temperature trajectories, and evaluating convergence of the REMD simulation.
 1. Extracting temperature trajectories


 The following script can be used with Cpptraj to read in the REMD trajectories (remd.mdcrd.001) and then write out trajectories at a single temperature. By using the “ensemble” command, Cpptraj will automatically recognize that it should look for trajectories from a REMD run. Denoting the lowest numbered replica (here, remd.mdcrd.001) is important, and Cpptraj will automatically look for the rest of the trajectory files with the same file name format, without having to specify how many replicas there are. In the example below, a distance is also calculated, and the output file will have the value of this distance for each corresponding temperature replica trajectory in each of 8 columns. This script assumes that Cpptraj is in your PATH.
 extract.sh
 ```
#!/bin/bash

cpptraj ala10.prmtop <<EOF
ensemble remd.mdcrd.001
trajout temp.crd nobox

distance d1 out dist.ensemble.dat :2@CA :8@CA

go
EOF

```
 Running this script will create the following files: 

 dist.ensemble.dat, temp.crd.0, temp.crd.1, temp.crd.2, temp.crd.3, temp.crd.4, temp.crd.5, temp.crd.6, temp.crd.7 

 The output from running the above script should show that Cpptraj successfully recognized the eight REMD replicas and is processing the ensemble.
 ```
CPPTRAJ: Trajectory Analysis. V6.29.2 (GitHub)
    ___  ___  ___  ___
     | \/ | \/ | \/ |
    _|_/\_|_/\_|_/\_|_

| Date/time: 07/02/25 12:11:43
| Available memory: 489.193 MB

        Reading 'ala10.prmtop' as Amber Topology
        Radius Set: H(N)-modified Bondi radii (mbondi2)
INPUT: Reading input from 'STDIN'
  [ensemble remd.mdcrd.001]
        Found 8 replicas.
        Reading 'remd.mdcrd.001' as Amber NetCDF
        Replica dimensions:
                1: Temperature
  [trajout temp.crd nobox]
        Writing ensemble member 'temp.crd.0' as Amber Trajectory
  [distance d1 out dist.ensemble.dat :2@CA :8@CA]
    DISTANCE: :2@CA to :8@CA, center of mass.
  [go]
---------- RUN BEGIN -------------------------------------------------

ENSEMBLE INFO:
  Ensemble size is 8
  Ensemble Indices Map:
        { 1 } -> 0
        { 2 } -> 1
        { 3 } -> 2
        { 4 } -> 3
        { 5 } -> 4
        { 6 } -> 5
        { 7 } -> 6
        { 8 } -> 7

```
 2. Post-processing replica trajectories


 To evaluate sampling across the replicas, a user may want to post-process all of the replica trajectories. In the distance.sh script below, we use the ‘ensemble’ command in Cpptraj to calculate our distance for each of the replica trajectories, and print the output to a file named “distance.replicas.dat.” We add the “nosort” keyword to the end of the ensemble command to tell Cpptraj we do not want to sort the trajectories by temperature. Running this script will create the following files:


 dist.ensemble.dat.
 distance.sh
 ```
#!/bin/bash

cpptraj ala10.prmtop <<EOF
ensemble remd.mdcrd.001 nosort

distance d1 out dist.replicas.dat :2@CA :8@CA

go
EOF

```
 Here is the histogram plot of the dist.ensemble.dat and dist.replicas.dat:
![Histogram plot of the dist.ensemble.dat and dist.replicas.dat with 1000 exchanges](https://ambermd.org/tutorials/advanced/tutorial43/include/Fig1.jpg)
 **Figure 1.** *Histogram plot of the dist.ensemble.dat and dist.replicas.dat with 1000 exchanges*
 Remember, we only ran 1000 exchanges (numexchg=1000), and exchanged every 1 ps (nstlim=500, dt=0.002), so this is 1 ns of dynamics per replica. Although this is a short simulation, we can still see that the higher temperatures sample longer distances, consistent with protein unfolding. However, a look at the replica trajectories shows that each replica samples different distances. In a converged simulation, all replicas should sample the same distribution to call the simulations converged in the measured distribution. An example of a converged simulation histogram is provided below.
 ![Histogram plot of the dist.ensemble.dat and dist.replicas.dat with 1,000,000 exchanges](https://ambermd.org/tutorials/advanced/tutorial43/include/Fig2.jpg)
 **Figure 2.** *Histogram plot of the dist.ensemble.dat and dist.replicas.dat with 1,000,000 exchanges*
 This shows the same system, where exchanges were attempted every 1 ps, but instead of 1000 exchanges 1,000,000 exchanges were performed (1 μs of sampling per replica).
 A user may want to ONLY extract a specific temperature from the ensemble. To do so, one can use the following script, extract_300K.sh.
 Extract_300k.sh
 ```
#!/bin/bash

cpptraj ala10.prmtop <<EOF
trajin remd.mdcrd.001 remdtraj remdtrajtemp 303.23
trajout temp.303.23.crd nobox

go
EOF

```
 If the user diffs the temp.303.23.crd trajectory with temp.crd.1 output from using the ensemble command to sort the temperature replicas, there should be NO DIFFERENCES. This is a convenient way to extract single temperatures in case systems are large, explicitly solvated, and/or the user doesn’t want to use the additional storage space to store all temperature trajectories.
 This tutorial was written by Christina Bergonzo and Akanksha Manghrani. 

 Previous versions were written by Dan Roe, Asim Okur and Carlos Simmerling
 ## Disclaimer
 This article was funded by the National Institute of Standards and Technology and is not subject to U.S. Copyright. Certain commercial equipment, instruments, materials or software are identified in this paper to foster understanding. Such identification does not imply recommendation or endorsement by the National Institute of Standards and Technology, nor does it imply that the materials or equipment identified are necessarily the best available for the purpose.


