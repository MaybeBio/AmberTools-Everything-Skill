

# Simulating the Green Fluorescent Protein and Building a Modified Amino Acid Residue
 Formerly known as TUTORIAL B5
 **By Jason Swails, Dave Case, and Tai-Sung Lee**


   ![DNA](https://ambermd.org/tutorials/basic/tutorial5/files/gfp_img.jpg) The starting structure of the green fluorescent protein with the ligand molecule highlighted with a sticks representation.
  # Introduction
 This tutorial demonstrates how you would go about preparing a system with a modified amino acid in order to run a simulation using the AMBER molecular dynamics engines. Note that this tutorial differs from other examples that parametrize a small organic compound because the custom residue is part of a polymer chain, so there are additional complications involved with deriving the charges and preparing the files to appropriately insert the custom residue into the peptide chain.
 We assume that you are already familiar with running basic simulations of systems that contain only standard residues (familiarity with the [antechamber Sustiva tutorial](https://ambermd.org/tutorials/basic/tutorial4b/index.html) is suggested as well), so that the focus here will be on the steps needed to parametrize a custom residue. The tutorial will, however, follow through on the full workflow of a molecular dynamics investigation of GFP.
 The basic outline of this tutorial is as follows:
1. Preparing the Protein Data Bank (PDB) file for use with AMBER
2. Computing partial charges and atom types of the custom residue
3. Preparing the residue library and force field parameters for use with LEaP
4. Creating the topology and coordinate files for the simulation
5. Minimizing, heating, equilibrating, and running production dynamics using the generated files
 Let's get started (below)!
  # Section 1: Preparing the PDB file
 We start with the PDB file [1EMA.pdb](http://www.rcsb.org/pdb/explore/explore.do?structureId=1ema), which contains the GFP with a fluorophore named CRO. So go ahead and download the PDB file from the RCSB website (or [click here to download a copy](https://ambermd.org/tutorials/basic/tutorial5/files/1EMA.pdb)). You can see the CRO residue in line 896 of the PDB file as residue 66 (residues 65 and 67 joined with residue 66 in the creation of the fluorophore, as described in the header of the PDB file).
 The first stage is to "massage" the PDB file to prepare it for use in tleap. Amber includes a script called pdb4amber that will modify the PDB file for use with tleap. The --help flag will show you all of the available options. We will use the --dry flag to remove crystallographic waters and the --reduce flag to add hydrogen atoms in their optimal locations:
 pdb4amber -i 1EMA.pdb -o gfp.pdb --dry --reduce
 Pay careful attention to the output from pdb4amber so that you understand what has been done. The resulting PDB file should now contain all of the hydrogens for each of the residues (and there are no disulfide bonds, but if there were pdb4amber would have added them correctly for you). After you run *pdb4amber*, you need to modify the resulting gfp.pdb file to convert the MSE residues (selenomethionine) into the more standard MET (we intend to simulate the wild-type variant, not the mutated version needed to ensure proper crystallization). You can do this with *sed* or your favorite editor by doing a global replace of MSE with MET. Furthermore, you need to convert the selenium atom into the sulfur atom (named SD) and change the element from SE to S. You are encouraged to run this command yourself, but our copy of [gfp.pdb can be downloaded here](https://ambermd.org/tutorials/basic/tutorial5/files/gfp.pdb) for comparison.
  # Section 2: Computing partial charges and atom types for CRO
 So we know from the header of 1EMA.pdb that CRO is not a standard residue (it is, in fact, a cyclization of three amino acids). As a result, the standard Amber residue libraries do not contain atom types or charges for this residue. These are needed for any molecular mechanical simulation. This section will focus on deriving charges and determining the atom types of the CRO residue.
 We will use *antechamber* to derive charges using the *bcc* charge scheme that you will be familiar with if you ran the Sustiva tutorial mentioned in the introduction above. Note that this choice is not unique, nor necessarily optimal for every application. Other tools, such as [R.E.D. Tools](http://q4md-forcefieldtools.org/) provide a more rigorous and reproducible set of charges, but antechamber will suffice for our needs in this tutorial. If you wish to use R.E.D., they have a number of tutorials that you can use to generate the charges and atom types, after which you can rejoin this tutorial in Section 4.
 The CRO template employed herein is from the very useful [Chemical Component Dictionary (CCD)](https://www.wwpdb.org/data/ccd) that contains observed and idealized structures for all residue and small molecule components that appear in any PDB entry (components.cif, i.e. the CCD, can be [downloaded](https://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz) as a multi megabyte gzip file). Using the CCD avoids the coordinate construction task and ensures that the atom and residue names will match those in the PDB. If you search for "CRO" in the PDB then it will bring you to a site where you can download the CRO component definition file in CIF format, CRO.cif. ( [You can also download a copy here](https://ambermd.org/tutorials/basic/tutorial5/files/CRO.cif).)
 The antechamber program can read the component CIF (ccif) files and generate charges and atom types. You should consult the "Antechamber and GAFF" chapter of the Amber reference manual for details and examples. You can also type antechamber -help at the command-line for information about the various options. Here we will ask for Amber atom types (-at amber), since this is a modified amino acid, and resembles the molecules for which the Amber atom types have been parametrized. (For more general organic molecules, it is usually better to use gaff atom types.) We use the following command to run *antechamber* (which should complete within a couple minutes):
 antechamber -fi ccif -i CRO.cif -bk CRO -fo ac -o cro.ac -c bcc -at amber
 Among many other files, you will see cro.ac, which vaguely resembles a PDB file, but contains a listing of bonds as well as partial atomic charges and atom types. Note that antechamber is mostly set up to use gaff atom types, so it occasionally makes a mistake with the Amber atom types. In this example, we will fix up the beginning *N* atom, which should have the same type as an amide nitrogen (N). You can simply change the atom type from NT to N by hand. (Or for practice, try to use *sed*, *awk*, or *perl* to change it.) You can [download a copy of cro.ac](https://ambermd.org/tutorials/basic/tutorial5/files/cro.ac) that we created to compare to yours.
 We now have the partial atomic charges and atom types for the atoms within the antechamber molecule file. The next section will discuss how to create the residue library file and prepare it for use with LEaP.
 NOTE: you would generally use *antechamber* to create a mol2 file to be read into LEaP. However, since we need to process the resulting file with *prepgen* later, we need to use the *ac* format which is all that *prepgen* supports.
  # Section 3: Preparing the residue library and force field parameters for use with LEaP
 The "CRO" component we downloaded from the PDB is a complete molecule *with all of its hydrogen atoms included*, which is what antechamber (or *any* quantum mechanical program) needs in order to compute charges. But we need to strip off the atoms at the beginning and the end to make an "amino acid"-like residue that is ready to be connected to other amino acids at its N- and C-termini.
 This can be done using the *prepgen* program with an "mc" ( *m*ain *c*hain) file that identifies the atoms to be removed as well as which atoms are part of the main chain (i.e., backbone atoms). You generally need to create this file yourself. Here are the contents for this particular example:
 HEAD_NAME N1


 TAIL_NAME C3 

 MAIN_CHAIN CA1


 MAIN_CHAIN C1


 MAIN_CHAIN N3


 MAIN_CHAIN CA3


 OMIT_NAME H2


 OMIT_NAME HN11


 OMIT_NAME OXT


 OMIT_NAME HXT


 PRE_HEAD_TYPE C


 POST_TAIL_TYPE N


 CHARGE 0.0
 The HEAD_NAME and TAIL_NAME lines identify the atoms that will connect to the previous and following amino acids, respectively. The MAIN_CHAIN lines list the atoms along the chain that connect the *head* and the *tail* atoms. The OMIT_NAME lines list the atoms in the CRO residue that should be removed from the final structure, as they are not present in the intact protein. The PRE_HEAD_TYPE and POST_TAIL_TYPE lines let prepgen know what atom types in the surrounding protein will be used for the covalent connection. (Remember that *prepgen* can be used for polymers other than just proteins and peptides.) The CHARGE line gives the total charge on the residue; *prepgen* will ensure that the charges of the "omitted" atoms are redistributed among the remaining atoms so that the total charge is correct (i.e., 0 in this case).
 Again, you can use the the -help flag with *prepgen* to get a listing of all available options for the program. The following command runs *prepgen* to strip out the unneeded atoms at the N- and C-termini, assuming that you named the above mainchain file cro.mc:
 prepgen -i cro.ac -o cro.prepin -m cro.mc -rn CRO
 After you run this command, you should have a cro.prepin file that contains the definition of the CRO residue (you can [download our copy of the file here](https://ambermd.org/tutorials/basic/tutorial5/files/cro.prepin) for comparison).
 At this point, we have our residue library containing the charges for the atoms in our modified amino acid . We next need to check that the covalent parameters (bonds, angles, and dihedrals) that will be needed are available. The *parmchk2* program figures out what parameters will be needed and checks to see if they are in the standard files. If not, it tries to make educated guesses, and puts these new parameters into a file we are calling "frcmod.cro" here.
 To run *parmchk2*, use the following command (again, the -help flag will list all available options):
 parmchk2 -i cro.prepin -f prepi -o frcmod.cro -a Y \


 -p $AMBERHOME/dat/leap/parm/parm10.dat
 Here we specify the parm10.dat file because it is the main parameter database for the force field we plan to use, *ff14SB* (this information can be found in the leaprc file for the force field you plan to use). If we omit this flag, then *parmchk2* will match to parameters in the *gaff* database, which is not what we wanted when we chose the *amber* atom types when we ran *antechamber*. We also asked for *all* parameters to be printed (even the ones with perfect matches in parm10.dat), for reasons that will be apparent soon.
 After this step is done, look at frcmod.cro (you can compare [to our copy of frcmod.cro](https://ambermd.org/tutorials/basic/tutorial5/files/frcmod.cro)). You should immediately see problems with lines containing parameters of 0's labeled ATTN, need revision. This means that *parmchk2* could not find suitably similar parameters in the parm10.dat database. There are also numerous other parameters that were selected having a high penalty (which indicates a poor fit as determined by *parmchk2*).
 A simple solution that is likely to yield sensible parameters is to request that *parmchk2* searches through *gaff.dat* to "fill in" the parameters that were either not represented at all, or poorly represented, by the options in parm10.dat. To do this, we will delete the parameters that say ATTN: need revision in frcmod.cro and tell *parmchk2* to look in gaff.dat (default file) for parameters. The parameters from parm10.dat that we *want* to keep instead of the gaff parameters must be specified in the frcmod1.cro file (this is why we needed to print *all* of them with the -a Y flag). This is done by the set of commands below:
 grep -v "ATTN" frcmod.cro > frcmod1.cro # Strip out ATTN lines


 parmchk2 -i cro.prepin -f prepi -o frcmod2.cro
 At this point, we have two frcmod files; one with the parameters extracted from the Amber parameter database (frcmod1.cro, which you can [download here to compare to your file](https://ambermd.org/tutorials/basic/tutorial5/files/frcmod1.cro)) and one with the parameters extracted by matching to gaff atom types (frcmod2.cro, which you can [download here to compare to your file](https://ambermd.org/tutorials/basic/tutorial5/files/frcmod2.cro)). Between these two files, we have all of the parameters that we need. If you want to, you can combine these into a single frcmod file by extracting the 7 missing parameters (and any additional parameters with high penalties) from frcmod2.cro and adding them to frcmod1.cro. However, if you are careful to follow the order in which these parameters are loaded in the next section, this step is not needed. Let's proceed on to our step of creating the prmtop and inpcrd files.
 As always, the parameters generated here--particularly the ones supplied by gaff, should be viewed as a starting point that should be validated against available experimental or high-level QM data. For the purposes of this tutorial, we will simply proceed with what we have generated here.
  # Section 4: Creating the topology and coordinate files for simulation
 We now have all of the files we need to create the topology and coordinate files for sander or pmemd! We just need to load these files into LEaP to create these files. If you used R.E.D. or some other tool to do the charge derivation and atom typing, welcome back.
 For this example, we choose the ff14SB force field and the igb=8 implicit solvent model, as described in the Amber reference manual. Setting the default PBRadii to the "mbondi3" set is a part of specifying this force field.
 We then load the cro.prepin and the parameters created above. In this case, we need to load frcmod2.cro *first*, followed by *frcmod1.cro* to make sure that all of the gaff parameters are overwritten where we would rather use the *ff14SB* parameters instead. Then we pull in the 1EMA PDB file we prepared (as gfp.pdb), and save the results. This is shown in the script below, which we will save to tleap.in:
 source leaprc.protein.ff14SB


 set default PBRadii mbondi3


 loadAmberPrep cro.prepin


 loadAmberParams frcmod2.cro


 loadAmberParams frcmod1.cro


 x = loadPDB gfp.pdb


 saveAmberParm x gfp.parm7 gfp.rst7


 quit
 We can then run tleap with the following command:
 tleap -f tleap.in
 After this step, you should have your finished topology and coordinate files, and you are ready to run simulations! You can [download the versions we created here](https://ambermd.org/tutorials/basic/tutorial5/files/gfp_parmcrd.tgz) to compare with your own. On to section 5!
  # Section 5: Simulations; minimization, heating, equilibration, and production
 Since the purpose of this tutorial was to walk through parametrizing a modified polymer "link", this section will be brief. In your own project, you are of course free to choose an explicit solvent model and a more careful minimization, heating, and equilibration routine than what we use here (perhaps utilizing positional restraints to prevent the structure from distorting).
 ### Minimization
 The input file we use here for minimization is shown below. We will name this file min.in:
 simple generalized Born minimization script


 &cntrl


 imin=1, ntb=0, maxcyc=100, ntpr=10, cut=1000., igb=8, 

 /
 We can run the minimization with *sander* using the command:
 sander -O -i min.in -p gfp.parm7 -c gfp.rst7 -o min1.out -r min1.rst7
 As always, we suggest that you visualize the resulting structure (min1.rst7) to make sure that nothing obviously bad happened, as well as the output file to make sure that everything looks OK (e.g., that the structure remained intact and that the total energy and maximum gradient steadily decreased during the minimization). The output files we created can be downloaded as part of a tarball containing most of the files generated during the calculations at the end of this section.
 ### Heating
 The input file we use here for heating is shown below. We will name this file heat.in, and it will vary the target temperature linearly over the course of 200 ps from 10 K to 300 K.
 Implicit solvent initial heating mdin


 &cntrl


 imin=0, irest=0, ntx=1,


 ntpr=1000, ntwx=1000, nstlim=100000,


 dt=0.002, ntt=3, tempi=10,


 temp0=300, gamma_ln=1.0, ig=-1,


 ntp=0, ntc=2, ntf=2, cut=1000,


 ntb=0, igb=8, ioutfm=1, nmropt=1,


 /


 &wt


 TYPE='TEMP0', ISTEP1=1, ISTEP2=100000,


 VALUE1=10.0, VALUE2=300.0,


 /


 &wt TYPE='END' /
 Remember, we want to use the structure we generated in the last step (minimization) when we heat the structure. So the command to heat using *sander* will look something like this:
 sander -O -i heat.in -p gfp.parm7 -c min1.rst7 -o heat.mdout \


 -x heat.nc -r heat.rst7
 Note, using *sander* in serial may take a long time. We used *pmemd.cuda* to run our simulation, which finishes much faster when run on our GTX 680. And as always, check the resulting structure and trajectory with your favorite visualization tool to make sure nothing obviously bad happened. As before, the files we generated will be available in a compressed tarball at the end of this section
 ### Production
 We have successfully heated up our structure! Now we are ready for our production simulation. Note that what many people call "equilibration" is really just a portion of the production simulation that you neglect during analysis, either because you used restraints to stabilize the structure or you are letting the system migrate toward an equilibrium configuration. For the purposes of this tutorial, we will not distinguish between those two phases, as we would use the same input for both.
 The input file we will use is shown below, which we will call md.in:
 Implicit solvent molecular dynamics


 &cntrl


 imin=0, irest=1, ntx=5,


 ntpr=1000, ntwx=1000, nstlim=500000,


 dt=0.002, ntt=3, tempi=300,


 temp0=300, gamma_ln=1.0, ig=-1,


 ntp=0, ntc=2, ntf=2, cut=1000,


 ntb=0, igb=8, ioutfm=1,


 /
 You can use the following command to run the simulation (again, we will use *pmemd.cuda*:
 sander -O -i md.in -p gfp.parm7 -c heat.rst7 -o md1.mdout \


 -x md1.nc -r md1.rst7
 And that's it! Of course, you still need to analyze your simulation in order to test whatever hypothesis drove you to run the calculation in the first place. But you now know how to parametrize a new, modified monomeric unit within a polymer chain (such as a modified nucleotide or amino acid). The strategy employed here can be applied to *any* polymeric unit you wish to use.


