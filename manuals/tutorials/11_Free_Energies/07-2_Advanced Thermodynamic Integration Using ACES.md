

# Advanced Thermodynamic Integration (TI) Methods to Mutate a Ligand Bound to a Protein
 **by Shi Zhang**a 

 *aDarrin York group, Rutgers University, New Brunswick, NJ* 

 ## Learning Outcomes
 - [Understand recent Amber choices available for TI calculations.](#TICalc)
- [Install the most up-to-date version of `FE-Workflow`.](https://ambermd.org/tutorials/advanced/tutorial39/install_feworkflow.php)
- [Setup all initial files for a ligand mutation.](#Setup_System)
- [Set up an alchemical enhanced sampling (ACES) TI calculation using softcore potentials.](#FEWorkflow_Setup)
- [Run an ACES TI calculation with a single edge.](#Running)
- [Use FE-Toolkit to analyze the data from the TI calculation.](#FETools)
- For additional discussion and case studies, please visit the [Free Energy Workshop tutorial](https://ambertutorials-rutgerslbsr-c744272d5a9c1169e0dc9e19b8d800019105.gitlab.io/workshop/05_afeSims/04_rbfe.html).
 Implementation Notes:  # Introduction
 This tutorial will guide you through calculating the relative binding free energies (ΔΔGbinding) between two ligands bound to the same protein. The system will be the human Thr160-phospho CDK2/cyclin A protein bound to two different inhibitors (NU6094 - PDBID: 1H1Q and NU6086 - PDBID: 1H1R).**[ [1](#References)]** This is a relatively conservative chemical change from H to Cl with only one “edge”. An edge is one alchemical transformation between two molecules. If you have more ligands or multiple function group changes between two ligands, you may need to perform multiple alchemical changes, in which case you should use MCS to define the softcore atoms for each edge.**[ [2](#References)]**
 ![](https://ambermd.org/tutorials/advanced/tutorial39/include/1h1q_1h1r.png) Figure 1. Ligands NU6094 (1H1Q) and NU6086 (1H1R) bound to the the human Thr160-phospho CDK2/cyclin A protein.
 This tutorial assumes you know general Amber usage such as protein system setup, choices when relaxing a system and running production run molecular dynamics simulations. As well, it also assumes you know the fundamental issues surrounding free energy calculations. If you need more guidance, please see these tutorials:
 - [Building Protein Systems in Explicit Solvent](https://ambermd.org/tutorials/basic/tutorial7/index.php)
- [Relaxation of Explicit Water Systems](https://ambermd.org/tutorials/basic/tutorial13/index.php)
- [Running MD with pmemd](https://ambermd.org/tutorials/basic/tutorial14/index.php)
- [Thermodynamic Integration using softcore potentials](https://ambermd.org/tutorials/advanced/tutorial9/index.html)
 ## Advanced TI Calculations
 Free energies will be computed with thermodynamic integration (TI) with softcore potentials, a method used to ensure smooth free energy curves, as implemented in `pmemd.cuda` which runs on a GPU.**[ [3](#References)]** This implementation is a little more sophisticated than the standard TI. It will use `FE-Workflow` to create a number of inputs that will run the GPU-driven Hamiltonian Replica Exchange Molecular Dynamics (HREMD) TI calculation method known as alchemical enhanced sampling (ACES)[ **[4](#References)**] with recent flags available in Amber, including smoothstep softcore potentials and enhanced sampling methods.
 In this tutorial, we’re going to use `FE-Workflow` to set up the TI calculations, which will employ:
 - a second-order smoothstep softcore potential.
- enhanced sampling `git_add_sc` functionalities.
- Hamiltonian REMD with lambda TI as implemented in the ACES code.
 ## Noteworthy TI Flags
 ### Smoothstep Softcore Potentials
 Many issues with linear interpolation between end states result in spurious energies and poor statistics (reviewed in ref [ [5](#References)] introduction). Concerted softcore potentials “soften” the Lennard-Jones and electrostatic interactions with increasing values of α and β parameters, respectively. The default values[6](#References) in Amber are **`scalpha=0.5`** and **`scbeta=12`**. 

 

 A second-order smoothstep softcore potential, SSC(2), was developed to alleviate some of the remaining energy issues, replacing the original softcore potentials and setting variables **`scalpha=0.2`** and **`scbeta=50`**. The SSC(2) function with these values makes the free energy derivatives vanish at the transformation endpoints and allows smooth adjustment of the λ-dependent terms in the potential.[5](#References)
 ### Enhanced Sampling for Softcore Ligand Energies
 In theory, sampling will not affect the converged value of the ΔΔGbinding. However, in practice, for example when rotating a ligand from cis to trans, there may be steric hindrances in the protein binding pocket during calculation of the alchemical change. Enhanced sampling methods for softcore energy calculations have been developed to improve convergence. [4](#References)
![](https://ambermd.org/tutorials/advanced/tutorial39/include/gti_add_sc.png) Figure 2. Softcore enhanced sampling flags in Amber. For more explanation, see text in [ref 4](#References).
 If you think there is an issue with convergence, for instance, in your ΔΔGbinding calculation, you can scale interactions between the common core (e.g. pharmacophore) and your softcore region (the atoms that changed). Figure 2 lists `gti_add_sc` options, where S means that the interactions between the common core and softcore region will be scaled in lambda windows, P means energies will stay constant for all lambda windows, P* means only rotatable bonds are scaled with lambda. The default flag in Amber TI calculations is `git_add_sc=1`. Theoretically, you should obtain the same ΔΔGbinding with any of the 1-6 `gti_add_sc` values because G is a state function. However, in practice, the York group currently recommends `gti_add_sc=25`.
 ### ACES Method on the Ligand
 The alchemical enhanced sampling (ACES)[ **[4](#References)**] method is a HREMD Thermodynamic Integration method that performs exchanges between lambda windows. This increases the sampling between neighboring lambda windows and should increase convergence of the calculated ΔΔGbinding between trials. This will require running on multiple GPUs on a single node ( [see below](#GPU) for an estimation of how many).
 ### Reviews
 More information is also available through these articles. [ **[7, 8](#References)**]
# System Setup
 ## Prepare the Protein•Ligand PDB Files – 1H1Q and 1H1R
 1. Download the 1h1q structure Download the [1h1q](https://www.rcsb.org/structure/1H1Q). structure from the PDB.  2. Edit the pdb file to remove extraneous info. Open the 1h1q.pdb with any text editor. Remove everything but the protein-complex structure, i.e., remove the structural information at the beginning, remove all the solvent molecules and ions. Just use one monomer and one ligand. Remove any CONNECT info as well.  3. Add hydrogens. Use tleap to add H's to the protein.  Repeat *Steps 1-3* for the [1h1r](https://www.rcsb.org/structure/1H1R). protein•ligand complex.
 ## Derive Force Field Parameters for the Ligands with Antechamber
 [Force fields](https://ambermd.org/AmberModels.php) are specific for the molecule type. Ligands, usually organic molecules, are modelled in Amber with [GAFF](https://ambermd.org/AmberModels_organic.php) or some modification of this. For this tutorial, GAFF2 parameters will be derived with Antechamber. Each new ligand will require derivation of parameters with Antechamber unless they exist elsewhere in the literature.
 4. Add hydrogens to the ligand using a tool such a VMD or Chimera X. Using VMD to examine the ligand structure and add H atoms back to the ligand if missing. (The *Extension -> Modeling -> Molefacture* function can be used to add missing atoms). Save the structure.   ![ligand](https://ambermd.org/tutorials/advanced/tutorial39/include/1h1q_lig.png) Figure 3. Ligand from 1H1Q (NU6094) after hydrogen atoms added with *Molefacture*.
 5. Name the residue *LIG* Using a text editor, rename the residue name (columns 18 to 20) of the ligand from “2A6” to “LIG”




 The residue name for all ligands should be *“LIG”* in order for FE-Workflow to recognize the ligands.
  6. Make a new file for just the ligand. Using a text editor, copy the ligand to a new file named `1h1q_0.pdb`. Simply copy-paste. Do not cut-paste because you don't want to remove it from the original protein•ligand complex coordinates.
 7. Run *antechamber* to get a `mol2` file. This run will use GAFF2 parameters and set net charge = 0 because LIG is a neutral ligand. `ndiis_attempts=700` ⇒ Forces antechamber to run enough times to generate a force field. Here, {PDB} can be changed to "1h1q".
```
antechamber -i ${PDB}_0.pdb -fi pdb -o ${PDB}_0.mol2 -fo mol2 \
-at "gaff2" -pf y -an y -c 2 -nc 0 -rn LIG -ek "ndiis_attempts=700"

```
  8. Run *parmchk2* to get an `frcmod` file. Here, {PDB} can be changed to "1h1q".
```
parmchk2 -i ${PDB}_0.mol2 -o ${PDB}_0.frcmod -f mol2 -s 1

```
  9. Run *tleap* to generate an `lib` file. Here, {PDB} can be changed to "1h1q". `tleap` commands
```
loadamberparams ${PDB}_0.frcmod
LIG = loadmol2 ${PDB}_0.mol2
saveoff LIG ${PDB}_0.lib
quit

```
  10. Repeat step 4-9 for other ligands. For example, in this case, you want to repeat steps 4-9 for 1H1R where H is instead Cl.
 # Simulation Setup
 In this tutorial, we’re going to use `FE-Workflow` to set up the TI calculations, which will employ:
 - a second-order smoothstep softcore potential.
- enhanced sampling `git_add_sc` functionalities.
- Hamiltonian REMD with lambda TI as implemented in the ACES code.
 ## Setup directories and the `input` file for FE-Workflow
 `FE-Workflow` uses the text file called `input` (provided below) to specify parameters. It assumes the following: - A single place to put all the `.pdb` files and `.frcmod` files.
- Use the `input` file to control different steps of the TI calculation (ex: setup or analysis; in this first case setup).
 11. Make folders for the simulation. - Make a folder called `ti`.
- In the `ti` folder, make two folders, one called `rbfe` and another called `initial`.
 `FE-Workflow` will read files found in `initial/CDK2` and put results in the rbfe folder. In the `rbfe` folder, `FE-Workflow` will also make folders for each trial and lambda window, and then place TI input files, scripts, results in the correct folders. The `current` folder will contain the current coordinate file that is being simulated during production runtime. Figure 4 contains a list of files and/or folders the user needs to make (blue) and the ones the `FE-Workflow` will make (green).
 ![](https://ambermd.org/tutorials/advanced/tutorial39/include/file_setup.png) Figure 4. File setup for `FE-Workflow`.
 12. Edit the `input` file for FE-Workflow - Download the [`input`](https://ambermd.org/tutorials/advanced/tutorial39/include/rbfe/input) file or copy from the opened link.
- Place the `input` file in your `ti` folder.
 For this tutorial, we are going to use the [ff14SB protein force field](https://ambermd.org/AmberModels_proteins.php) and [gaff2 for ligand force field](https://ambermd.org/AmberModels_organic.php). The water model is set to tip4pew and the buffer size is 16 Å for complex system and 20 Å for solution system. We'll use `tleap` to add a water box and ions.
 `FE-Workflow` will build an extensive set of input and scripting files necessary to run a robust equilbiration, including solvent relaxation, heating and production run for 2 ns at each lambda step. This is all contained in the `run_alltrials.slurm` file. Default values are used in `FE-Workflow`. If you want to modify the relaxation, heating or production runs, you can manually modify the ` md.in` files for each step after `FE-Workflow` has made all the files.
 In this tutorial, we are running 11 lambda windows in parallel for 20 steps. After 20 steps, an exhange is attempted with a neighboring lambda window. At the exchange step, will try to exchange if successful, will run for another 20 steps and then another exchange. The flag `numexchg` controls how many exchanges will be attempted and `nstlim` specifies how many steps between exchange attempts.
 There is a detailed explanation for the flags inside [the sample `input` file](https://ambermd.org/tutorials/advanced/tutorial39/include/rbfe/input)and [Amber Reference Manual](https://ambermd.org/Manuals.php).
 ## Run FE-Workflow
 13. Go to the folder with `input`, type the command ```
setup_fe

```
 Workflow should begin to generate the files for you.
 ![setup](https://ambermd.org/tutorials/advanced/tutorial39/include/setup_fe1.png) In your transformation folder, `FE-Workflow` will make a set of folders, one for the ligand solvated in water ( `aq`) and for the protein•ligand complex ( `com`). It will copy over needed `.parm7` and `.rst7` files, as well as make `.sh` and slurm files.
 You will find simulation topology, coordinate, md input and slurm files under CDK2/1h1q~1h1r/aq and com folders respectively (see Figure 4). Files for the set of the [transformation 1h1q to 1h1r are provided](https://ambermd.org/tutorials/advanced/tutorial39/include/tutorial.mdout.tar.gz) for your reference.
 In the beginning, the trial folders (t1, t2 and t3) only contain `.rst7` files for end states (0 and 1). Later it will later write `.mdout`, `.nc` files as well as `.rst7` files for propagated lambda windows from equilibration.
## Running the TI Calculation
 14. The user should edit the slurm script to customize the number of GPUs requested. ### GPU Usage
 This TI implementation uses ACES, a replica exchange method, that requires multiple GPUs to run simultaneoulsy.
You might ask how many GPUs you need. You do not need the same number of GPUs as replicas. The number of replicas per GPU is determined by the GPU memory size. If your GPU has a lot of memory (ex: 32 GB), you can put 20 replicas on 1 GPU. A typical system with 100,000 atoms, will take 1 GB GPU memory. Modern GPUs typically have 8 GB memory. Therefore, about 8 replicas can go on 1 GPU simultaneously. The code will of course go faster if you have more GPUs. The slurm script will automatically distribute the replicas to the specified number of GPUs available. For instance, if you tell slurm to grab 2 GPUs, will put 5 in 1 GPU and 6 replicas in the other. If you have 4 GPUs, it will distribute the 11 replicas evenly. **NOTE: All the GPUs should be on the same node.** If you are using multiple nodes, it will slow down the exchange attempt due to I/O communication issues.
15. Submit the slurm script.
 Now transfer the files to any HPC where Amber is installed. Confirm the “#SBATCH --” lines are set appropriately. Submit the slurm script and wait for the simulations to finish. The simulation will go through minimization, equilibration and production stages. Each step will generate the output file and coordinate. Here is an example production run input.
 ```
&cntrl
imin            = 0
nstlim          = 20
numexchg        = 100000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 1
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = 10
iwrap           = 1

ntb             = 2
temp0           = 298.
ntt             = 3
gamma_ln        = 2.
tautp           = 1
ntp             = 1
barostat        = 2
pres0           = 1.01325
taup            = 5.0

ifsc            = 1
icfe            = 1

ifmbar          = 1
bar_intervall   = 1
mbar_states     = 11
mbar_lambda(1)     = 0.00000000
mbar_lambda(2)     = 0.10000000
mbar_lambda(3)     = 0.20000000
mbar_lambda(4)     = 0.30000000
mbar_lambda(5)     = 0.40000000
mbar_lambda(6)     = 0.50000000
mbar_lambda(7)     = 0.60000000
mbar_lambda(8)     = 0.70000000
mbar_lambda(9)     = 0.80000000
mbar_lambda(10)     = 0.90000000
mbar_lambda(11)     = 1.00000000

noshakemask     = '@1-90'

clambda         = 1.00000000
timask1         = '@1-45'
timask2         = '@46-90'
crgmask         = ""
scmask1         = '@41'
scmask2         = '@90'
scalpha         = 0.5
scbeta          = 1.0

gti_cut         = 1
gti_output      = 1
gti_add_sc      = 25
gti_scale_beta  = 1
gti_cut_sc_on   = 8
gti_cut_sc_off  = 10
gti_lam_sch     = 1
gti_ele_sc      = 1
gti_vdw_sc      = 1
gti_cut_sc      = 2
gti_ele_exp     = 2
gti_vdw_exp     = 2

gremd_acyc      = 1
/
&ewald
/

```
 # Analysis of the Simulation
 This tutorial uses `edgembar` (part of `FE-Toolkit`) within `AmberTools25` to analyze the data automatically.
If simulations finish without errors, you will find a series of `*ti.mdout` files in the t1/ folder. For the free energy results, we will extract energy values from these `*ti.mdout` files.
 ### Look at the trajectories of the endstates in VMD
 If you want to take a close look at the conformation or the structure of the ligand, load in VMD the unisc.parm7 and 0.00000000_ti.nc for real state ligand 1h1q (resid :1), 1.00000000_ti.nc for real state ligand 1h1r (resid :2).
 ### Perform Analysis
 - Transfer *ti.mdout files in both com/ and aq/ folders back to your local computer.
- Modify the “stage=analysis” and “path_to_data” to the folder containing the output files.
- Make sure you’ve loaded the `AmberTools` or `FE-Toolkit` module.
- Again, type the command
 ```
setup_fe
```
  For a more detailed discussion, please refer to [edgembar documentation](https://rutgerslbsr.gitlab.io/fe-toolkit/edgembar/index.html). Essentially, the analysis uses the `edgembar` program, and there are more advanced options for different features you may want to explore. An example of `edgembar` is provided [here](https://rutgerslbsr.gitlab.io/fe-toolkit/quickstart/#EDGES). A general usage is consisted of three steps:
 ```
# extract dvdl_*.dat and efep_*.dat from *_ti.mdout files to the ./dats folder
edgembar-amber2dats.py -o ./dats ./*_ti.mdout  

```
 ```
OOMP_NUM_THREADS=4 edgembar_omp --halves --fwdrev --no-auto 1h1q~1h1r.xml

```
 ```
python 1h1q~1h1r.py

```
 ```
edgembar-WriteGraphHtml.py -o Graph.html 1h1q~1h1r.py

```
 The `FE-Workflow` serves as an interface to connect `edgembar` and pass some variables for analysis.
 The experimental value to the binding free energy is not passed to the analysis mode in this example, but you can always apply the experimental restraints by providing a file consisting of two columns, first are the names of the complex, and second the experimental value in unit of kcal/mol. Let’s look at the graph.html file. Since we only have two nodes, no cycle is detected.  ![ligand](https://ambermd.org/tutorials/advanced/tutorial39/include/result.png) In the "Edges" table, it generates the relative binding free energy from 1h1q to 1h1r based on our one trial simulation. The ddG value is **-0.549 kcal/mol**, which agrees with reported value in ref[ **[4](#References)**]. The experimental value for such alchemical transformation is ~0.5 kcal/mol, therefore the calculated value roughly falls within the ±1 kcal/mol error range.
 ## References
 1. Davies, T. G., Bentley, J., Arris, C. E., Boyle, F. T., Curtin, N. J., Endicott, J. A., Gibson, A. E., Golding, B. T., Griffin, R. J., Hardcastle, I. R., Jewsbury, P., Johnson, L. N., Mesguiche, V., Newell, D. R., Noble, M. E. M., Tucker, J. A., Wang, L.,and Whitfield, H. J. [*Nature Structural and Molecular Biology*, 2002, *9*, 745-749.](https://www.nature.com/articles/nsb842)
2. Cao, Y., Jiang, T., and Girke, T. A maximum common substructure-based algorithm for searching and predicting drug-like compounds. [*Bioinformatics*, 2008, *24*, i366–i374.](https://academic.oup.com/bioinformatics/article/24/13/i366/236402)
3. Lee, T. S., Cerutti, D. S., Mermelstein, D., Lin, C., LeGrand, S., Giese, T.J., Roitberg, A., Case, D. A., Walker, R. C., and York, D. M. GPU-Accelerated Molecular Dynamics and Free Energy Methods in Amber18: Performance Enhancements and New Features. [*Journal of Chemical Information and Modeling*, 2018, *58*, 2043-2050.](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00462)
4. Lee, T. S., Tsai, H. C., Ganguly, A., and York, D. M. ACES: Optimized Alchemically Enhanced Sampling. [*Journal of Chemical Theory and Computation*, 2023, *19*, 472-487.](https://pubs.acs.org/doi/10.1021/acs.jctc.2c00697)
5. Steinbrecher, T.; Joung, I.; Case, D. A. Soft-core potentials in thermodynamic integration: Comparing one- and two-step transformations [*J. Comput. Chem.*, 2011, *32*, 3253– 3263.](https://onlinelibrary.wiley.com/doi/10.1002/jcc.21909)
6. Lee, T.S., Lin, Z., Allen, B.K., Lin, C., Radak, B.K., Tao, Y., Tsai, H.-C., Sherman, W. and York, D.M. Improved Alchemical Free Energy Calculations with Optimized Smoothstep Softcore Potentials [*J. Chem. Theory Comput.*, 2020, *16*, 5512–5525.](https://pubs.acs.org/doi/10.1021/acs.jctc.0c00237)
7. Tsai, H. C., Lee, T. S., Ganguly, A., Giese, T. J., Ebert, M. C., Labute, P., Merz Jr., K. M., and York, D. M. AMBER Free Energy Tools: A New Framework for the Design of Optimized Alchemical Transformation Pathways. [*Journal of Chemical Theory and Computation*, 2023, *19*, 640-658.](https://pubs.acs.org/doi/10.1021/acs.jctc.2c00725)
8. Ganguly, A., Tsai, H. C., Pendás-Fernández, M., Lee, T. S., Giese, T. J., and York, D. M. AMBER Drug Discovery Boost Tools: Automated Workflow for Production Free-Energy Simulation Setup and Analysis (ProFESSA). [*Journal of Chemical Information and Modeling*, 2022, *62*, 6069-6083.](https://pubs.acs.org/doi/10.1021/acs.jcim.2c00879)

