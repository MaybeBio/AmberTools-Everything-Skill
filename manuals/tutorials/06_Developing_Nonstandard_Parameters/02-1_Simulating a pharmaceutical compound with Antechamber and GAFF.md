

# Simulating a pharmaceutical compound using antechamber and the Generalized Amber Force Field
 **By Ross Walker and Sishi Tang** 

 *Updated for AMBER 18 by Michael Barton and Tyler Luchko*
**![Sustiva](https://ambermd.org/tutorials/basic/tutorial4b/figures/sustiva.png)**
## Learning Outcomes
 Using Antechamber to Create LEaP Input Files for SimulatingSustiva (efavirenz)-RT complex using the General Amber Force Field
 

 ## Introduction
 In this tutorial we will make use of the Antechamber package, which ships with AmberTools, to create prmtop and inpcrd files for the simulations of a protein-ligand complex, and perform a short GB simulation on it.  Antechamber is designed to be used with the "general AMBER force field (GAFF)"1. This force field has been specifically designed to cover most pharmaceutical molecules and is compatible with the traditional AMBER force fields in such a way that the two can be mixed during a simulation. Like the traditional AMBER force fields, GAFF uses a simple harmonic function form for bonds and angles but unlike the traditional protein and DNA orientated AMBER force fields the atom types used in GAFF are much more general such that they cover most of the organic chemical space. The current implementation of the GAFF force field consists of basic atom types and special atom types. The charge methods used can be HF/6-31G* RESP2 or AM1-BCC3.
 By design, GAFF, is a complete force field (so that missing parameters rarely occur), it covers almost all the organic chemical space that is made up of C, N, O, S, P, H, F, Cl, Br and I. Moreover, since GAFF is totally compatible with the AMBER macromolecular force fields it should prove to be a useful molecular mechanical tool for rational drug design. Especially in binding free energy calculations and molecular docking studies.
 The Antechamber tool set is designed to allow the rapid generation of topology files for use with the AMBER simulation programs. This is useful in situations where you want to automatically screen a large number of compounds. Allowing antechamber to calculate charges and atom types automatically for us using GAFF allows it to be included in shell scripts that process a large number of compounds. However, as with any automated system it is not perfect and if you plan on focusing on a single system you should consider manually assigning atom types and carefully validating things. With Antechamber, one may solve the following problems:
 1. Automatically identify bond and atom types
2. Judge atomic equivalence
3. Generate residue topology files
4. Find missing force field parameters and supply reasonable suggestions
 Remember though that Antechamber is not a replacement for due diligence. You should always closely examine the atom types that Antechamber assigns and verify to yourself that the choices are reasonable. You should never use scientific software in a "Black Box" approach!
 In this tutorial we shall use the Antechamber tools with LEaP to create topology and coordinate files for the prescription drug Sustiva (Efavirenz). [Efavirenz](http://en.wikipedia.org/wiki/Efavirenz)) is a human immunodeficiency virus type 1 (HIV-1) specific, non-nucleoside, reverse transcriptase (RT) inhibitor marketed by Bristol Myers Squibb for controlling the progression of HIV infection in humans. The chemical name for Sustiva is (S)-6-chloro-(cyclopropylethynyl)-1,4-dihydro-4-(trifluoromethyl)-2H-3,1-benzoxazin-2-one. Its empirical formula is C14H9ClF3NO2 and it's 2D structure is:
![sustiva 2d](https://ambermd.org/tutorials/basic/tutorial4b/files/sustiva_2d.gif)
 Here is a basic 3-dimensional geometry for Sustiva from which we will start to build our topology and coordinate files: [sustiva.pdb](https://ambermd.org/tutorials/basic/tutorial4b/files/sustiva.pdb). Our Sustiva pdb file is extracted from the RT-sustiva complex pdb file (PDB ID [1FKO](http://www.rcsb.org/pdb/explore/explore.do?structureId=1FKO) ). The coordinates of sustiva are associated with a residue called "EFZ" (Efavirenz). 

By all means open it up in VMD and take a look at it.
We shall use Antechamber to assign atom types to this molecule and also calculate a set of point charges. Antechamber is the most important program within the set of Antechamber tools. It can perform many file conversions and can also assign atomic charges and atom types. Depending on its inputs, antechamber executes the following programs (all provided with AmberTools): *sqm*, *atomtype*, *am1bcc*, *bondtype*, *espgen*, *respgen* and *prepgen*. It will also generate a series of intermediate files (all in capital letters).
 ## Create parameter and coordinate files for Sustiva
 First, let's run *reduce* to add all the hydrogen atoms to the pdb file. The hydrogenated sustiva coordinates can be found in [sustiva_h.pdb](https://ambermd.org/tutorials/basic/tutorial4b/files/sustiva_h.pdb).
 reduce sustiva.pdb > sustiva_h.pdb
 (A note: *reduce* relies on a database of known ligands that are found in the PDB. If you have a molecule that is not already in the PDB component database, and you need to add hydrogens, give *OpenBabel* a try:
 obabel sustiva.pdb -O sustiva_h.pdb -h
 Be sure to check that either *reduce* or *obabel* have made the correct changes. And note that *OpenBabel* is not distributed with Amber.)
 To be consistent with the name of pdb, we will change the name of the residue from "EFZ" to "SUS", and create a new pdb file: [sustiva_new.pdb](https://ambermd.org/tutorials/basic/tutorial4b/files/sustiva_new.pdb).
Now let's try using antechamber on our sustiva pdb file. To create the "mol2" file, required to define a new unit in LEaP, we simply run the following command:
antechamber -i sustiva_new.pdb -fi pdb -o sustiva.mol2 -fo mol2 -c bcc -s 2
 Here the *-i sustiva.pdb* specifies the name of the 3D structure file and the *-fi pdb* tells antechamber that this is a pdb format file (we could easily have used any number of other supported formats including Gaussian Z-Matrix [gzmat], Gaussian Output [gout], MDL [mdl], amber Restart [rst], Sybyl Mol2 [mol2]). The *-o sustiva.mol2* specifies the name of our output file and the *-fo mol2* states that we want the output file to be of Tripos Mol2 format (this is a format supported by LEaP). The *-c bcc* option tells antechamber to use the AM1-BCC charge model in order to calculate the atomic point charges while the *-s 2* option defines the verbosity of the status information provided by antechamber. In this case we have selected verbose output (2).
 So, go ahead and run the above command. The screen output should be as follows:
```
Running: /usr/local/amber10/bin/bondtype -j full -i ANTECHAMBER_BOND_TYPE.AC0 -o ANTECHAMBER_BOND_TYPE.AC -f acRunning: /usr/local/amber10/bin/atomtype -i ANTECHAMBER_AC.AC0 -o ANTECHAMBER_AC.AC -p gaffTotal number of electrons: 160; net charge: 0Running: /usr/local/amber10/bin/sqm.shRunning: /usr/local/amber10/bin/am1bcc -i ANTECHAMBER_AM1BCC_PRE.AC -o ANTECHAMBER_AM1BCC.AC -f ac-p /usr/local/amber10/dat/antechamber/BCCPARM.DAT -s 2 -j 1Running: /usr/local/amber10/bin/atomtype -f ac -p bcc -o ANTECHAMBER_AM1BCC.AC -i ANTECHAMBER_AM1BCC_PRE.AC
``` You should also get a whole series of files written to your directory.
```
ANTECHAMBER_AC.AC      ANTECHAMBER_AM1BCC_PRE.AC  ATOMTYPE.INF  sqm.out  sustiva.mol2ANTECHAMBER_AC.AC0     ANTECHAMBER_BOND_TYPE.AC0  divcon.pdb    sqm.pdb  sustiva.pdbANTECHAMBER_AM1BCC.AC  ANTECHAMBER_BOND_TYPE.AC   sqm.in     
``` The files in CAPITALS are all intermediate files used by antechamber and are not required here. You can safely delete them. These files are not deleted by default since they may be of interest if things didn't work correctly. The sqm.xxx files are input and output from the sqm quantum mechanics code used by Antechamber to calculate the atomic point charges. We are not interested in the data here except to check that the sqm calculation completed successfully:
 [sqm.out](https://ambermd.org/tutorials/basic/tutorial4b/files/mopac.out)[](https://ambermd.org/tutorials/basic/tutorial4b/files/sqm.out) The last line of your sqm.out file should read:


 ************* Calculation Completed **************
 The file that we are really interested in, and the reason we ran Antechamber in the first place, is the [sustiva.mol2](https://ambermd.org/tutorials/basic/tutorial4b/files/sustiva.mol2) file. This contains the definition of our sustiva residue including all of the charges and atom types that we will load into LEaP to when creating our prmtop and rst7 files. Let's take a quick look at the file:
 [sustiva.mol2](https://ambermd.org/tutorials/basic/tutorial4b/files/sustiva.mol2) @<TRIPOS>MOLECULE


SUS


30 32 1 0 0


SMALL


bcc






@<TRIPOS>ATOM


1 CL -4.6850 -32.7250 25.2220 cl 999 SUS -0.073100


2 F1 -0.7550 -36.6320 25.6970 f 999 SUS -0.231400


3 F2 1.0780 -37.0430 24.6720 f 999 SUS -0.221400


4 F3 -0.7840 -37.1770 23.6260 f 999 SUS -0.216800


5 O1 1.5240 -34.9340 20.9100 o 999 SUS -0.573600


6 O2 0.9890 -34.8800 23.0580 os 999 SUS -0.371900


7 N -0.6810 -34.9710 21.4340 n 999 SUS -0.459500


8 C1 -1.6620 -34.4140 22.3130 ca 999 SUS 0.084700


9 C2 -2.9150 -33.9470 21.8430 ca 999 SUS -0.167000


10 C3 -3.8380 -33.4230 22.7710 ca 999 SUS -0.069500


11 C4 -3.5330 -33.3730 24.1190 ca 999 SUS -0.025500


12 C5 -2.3100 -33.8290 24.5930 ca 999 SUS -0.040200
 As you can see this file contains the 3 dimensional structure of our sustiva molecule as well as the charge on each atom, final column, the atom number (column 1), its name (column 2) and it's atom type (column 6). It also specifies the bonding at the end of the file. This file does not, however, contain any parameters. The GAFF parameters are all defined in $AMBERHOME/dat/leap/parm/gaff.dat. The other thing you should notice here is that all of the GAFF atom types are in lower case. This is the mechanism by which the GAFF force field is kept independent of the macromolecular AMBER force fields. All of the traditional AMBER force fields use uppercase atom types. In this way the GAFF and traditional force fields can be mixed in the same calculation.
 While the most likely combinations of bond, angle and dihedral parameters are defined in the parameter file it is possible that our molecule might contain combinations of atom types for bonds, angles or dihedrals that have not been parameterised. If this is the case then we will have to specify any missing parameters before we can create our prmtop and rst7 files in LEaP.
 We can use the utility parmchk2 to test if all the parameters we require are available.
 parmchk2 -i sustiva.mol2 -f mol2 -o sustiva.frcmod
 Run this command now and it will produce a file called [sustiva.frcmod](https://ambermd.org/tutorials/basic/tutorial4b/files/sustiva.frcmod). This is a parameter file that can be loaded into LEaP in order to add missing parameters. Here it will contain all of the missing parameters. If it can antechamber will fill in these missing parameters by analogy to a similar parameter. You should check these parameters carefully before running a simulation. If antechamber can't empirically calculate a value or has no analogy it will either add a default value that it thinks is reasonable or alternatively insert a place holder (with zeros everywhere) and the comment "ATTN: needs revision". In this case you will have to manually parameterise this yourself. It is the hope that as GAFF is developed, the number of missing parameters will decrease. Let's look at our frcmod file:
 [sustiva.frcmod](https://ambermd.org/tutorials/basic/tutorial4b/files/sustiva.frcmod) ```
remark goes hereMASSBONDANGLEca-c3-c1   64.784     110.735   Calculated with empirical approachc1-c1-cx   56.400     177.990   same as c1-c1-c3c1-cx-hc   48.300     109.750   same as c1-c3-hcc1-cx-cx   64.200     111.590   same as c1-c3-c3DIHEIMPROPERca-ca-ca-ha         1.1          180.0         2.0          General improper torsional angle (2 general atom types)n -o -c -os        10.5          180.0         2.0          General improper torsional angle (2 general atom types)c -ca-n -hn         1.1          180.0         2.0          General improper torsional angle (2 general atom types)ca-ca-ca-n          1.1          180.0         2.0          Using default valueNONBON
```
 We can see that there were a total of 4 missing angle parameters and 4 missing improper dihedrals. For the purposes of this tutorial we shall assume that the parameters Antechamber has suggested for us are acceptable. Ideally you should really test these parameters (by comparing to *ab initio* calculations for example) to ensure they are reasonable. If you see any parameters listed with the comment "ATTN: NEEDS REVISION" then it means that Antechamber could not determine suitable parameters and so you must manually provide these before you can proceed with the simulation. By default Antechamber will have set the values to zero.
 We now have everything we need to load sustiva as a unit in LEaP. We just need to run tleap and ensure the GAFF force field is available.
$tleap -f oldff/leaprc.ff99SB
 Once tleap is up and running we also need to ensure that it knows about the GAFF force field. There is a script in $AMBERHOME/dat/leap/cmd/ that will do this for us. We can load it into tleap with:
 source leaprc.gaff
 Our tleap console should now look something like this:
Welcome to LEaP!


(no leaprc in search path)


Sourcing: /usr/local/amber10/dat/leap/cmd/oldff/leaprc.ff99SB


Log file: ./leap.log


Loading parameters: /usr/local/amber10/dat/leap/parm/parm99.dat


Reading title:


PARM99 for DNA,RNA,AA, organic molecules, TIP3P wat. Polariz.& LP incl.02/04/99


Loading parameters: /usr/local/amber10/dat/leap/parm/frcmod.ff99SB


Reading force field modification type file (frcmod)


Reading title:


Modification/update of parm99.dat (Hornak & Simmerling)


Loading library: /usr/local/amber10/dat/leap/lib/all_nucleic94.lib


Loading library: /usr/local/amber10/dat/leap/lib/all_amino94.lib


Loading library: /usr/local/amber10/dat/leap/lib/all_aminoct94.lib


Loading library: /usr/local/amber10/dat/leap/lib/all_aminont94.lib


Loading library: /usr/local/amber10/dat/leap/lib/ions94.lib


Loading library: /usr/local/amber10/dat/leap/lib/solvents.lib


> source leaprc.gaff 

----- Source: /usr/local/amber10/dat/leap/cmd/leaprc.gaff


----- Source of /usr/local/amber10/dat/leap/cmd/leaprc.gaff done


Log file: ./leap.log


Loading parameters: /usr/local/amber10/dat/leap/parm/gaff.dat


Reading title:


AMBER General Force Field for organic mol., add. info. at the end (June, 2003)


>
Now we can load our sustiva unit (sustiva.mol2):
SUS = loadmol2 sustiva.mol2
If you now type list in tleap you should see the new SUS unit (highlighted in bold):
> SUS = loadmol2 sustiva.mol2


Loading Mol2 file: ./sustiva.mol2


Reading MOLECULE named SUS


> list


ACE ALA ARG ASH ASN ASP CALA CARG


CASN CASP CCYS CCYX CGLN CGLU CGLY CHCL3BOX


CHID CHIE CHIP CHIS CILE CIO CLEU CLYS


CMET CPHE CPRO CSER CTHR CTRP CTYR CVAL


CYM CYS CYX Cl- Cs+ DA DA3 DA5


DAN DC DC3 DC4 DC5 DCN DG DG3


DG5 DGN DT DT3 DT5 DTN GLH GLN


GLU GLY HID HIE HIP HIS HOH IB


ILE K+ LEU LYN LYS Li+ MEOHBOX MET


MG2 NALA NARG NASN NASP NCYS NCYX NGLN


NGLU NGLY NHE NHID NHIE NHIP NHIS NILE


NLEU NLYS NMABOX NME NMET NPHE NPRO NSER


NTHR NTRP NTYR NVAL Na+ PHE PL3 POL3BOX


PRO QSPCFWBOX RA RA3 RA5 RAN RC RC3


RC5 RCN RG RG3 RG5 RGN RU RU3


RU5 RUN Rb+ SER SPC SPCBOX SPCFWBOX SPF


SPG SUS T4E THR TIP3PBOX TIP3PFBOX TIP4PBOX TIP4PEWBOX


TP3 TP4 TP5 TPF TRP TYR VAL WAT frcmod99SBgaff parm99
At this point we haven't loaded the frcmod file that parmchk2 gave us. Thus if we check our SUS unit we should find that there are 4 missing angle type parameters.
check SUS
> check SUS


Checking 'SUS'....


Checking parameters for unit 'SUS'.


Checking for bond parameters.


Checking for angle parameters.


Could not find angle parameter: ca - c3 - c1


Could not find angle parameter: c1 - c1 - cx


Could not find angle parameter: c1 - cx - hc


Could not find angle parameter: c1 - cx - cx


Could not find angle parameter: c1 - cx - cx


There are missing parameters.


Unit is OK.
Our missing angle type parameters were ca-c3-c1, c1-c1-cx, c1-cx-hc and c1-cx-cx. These correspond to the propyl ring and the c-c triple bond. This is what we would expect since this type of system is fairly rare in organic molecules. We can now load our frcmod file in order to tell tleap the parameters for these missing angle types.
 loadamberparams sustiva.frcmod
 If we now check out SUS unit we should find that there are no missing parameters:
> loadamberparams sustiva.frcmod loading parameters: ./sustiva.frcmod


Reading force field modification type file (frcmod)


Reading title:


remark goes here


> check SUS


Checking 'SUS'....


Checking parameters for unit 'SUS'.


Checking for bond parameters.


Checking for angle parameters.


Unit is OK. We can now create the library file for sustiva ( [sus.lib](https://ambermd.org/tutorials/basic/tutorial4b/files/sus.lib)), as well as the prmtop and rst7 files ( [sustiva.prmtop](https://ambermd.org/tutorials/basic/tutorial4b/files/sustiva.prmtop), [sustiva.rst7](https://ambermd.org/tutorials/basic/tutorial4b/files/sustiva.rst7)).
 saveoff SUS sus.lib
saveamberparm SUS sustiva.prmtop sustiva.rst7
 The output from tleap shows a few warnings, which can be safely ignored (in this case!) due to the triangular bond geometry of sustiva:
> saveoff SUS sus.lib


Building topology.


Building atom parameters.


>


> saveamberparm SUS sustiva.prmtop sustiva.inpcrd


Checking Unit.


Building topology.


Building atom parameters.


Building bond parameters.


Building angle parameters.


Building proper torsion parameters.


1-4: angle 7 12 duplicates bond ('triangular' bond) or angle ('square' bond)




1-4: angle 7 9 duplicates bond ('triangular' bond) or angle ('square' bond)




1-4: angle 9 12 duplicates bond ('triangular' bond) or angle ('square' bond)




Building improper torsion parameters.


total 8 improper torsions applied


Building H-Bond parameters.


Not Marking per-residue atom chain types.


Marking per-residue atom chain types.


(Residues lacking connect0/connect1 -


these don't have chain types marked:




res total affected




SUS 1


)


(no restraints)


>
Instead of typing everything in the tleap console, the list of tleap commands mentioned above can also be used to create a tleap input file ( [tleap.in](https://ambermd.org/tutorials/basic/tutorial4b/files/tleap.in)) and generate all the required files:
tleap -f tleap.in
 ## Creating topology and coordinate files for Sustiva-RT complex
 Since we can mix the traditional AMBER force fields with GAFF, we can at this point load a fragment of the reverse transcriptase (RT) from HIV virus and treat this using the ff99SB force field while treating the Sustiva molecule using the GAFF force field. We will need to use the Sustiva library file ( [sus.lib](https://ambermd.org/tutorials/basic/tutorial4b/files/sus.lib)) that was created in the previous step.
The RT-Sustiva complex can be found in the RCSB protein data bank (pdb code 1FKO). The corresponding pdb file is [1FKO.pdb](https://ambermd.org/tutorials/basic/tutorial4b/files/1FKO.pdb).
HIV reverse transcriptase is a heterodimer composed of the p51 and p66 subunits. It is a large protein with a molecular mass of 117 kDa. For the purpose of this tutorial, we will use a truncated system in close proximity of Sustiva, including the finger and palm domains of the p66 subunit. The truncated pdb file is: [1FKO_trunc.pdb](https://ambermd.org/tutorials/basic/tutorial4b/files/1FKO_trunc.pdb).
In order to have the 1FKO pdb file recognized by tleap once the Sustiva library file is loaded, we need to change the residue name in the 1FKO pdb file from EFZ to SUS. Since the Sustiva libary file includes the same atom names as the Sustiva molecule in the 1FKO pdb file, no further modification is necessary. In other cases it is always a good idea to check the pdb file against the library file to make sure they have matching atom names as well as residue name. The modified pdb file is now: [1FKO_trunc_sus.pdb](https://ambermd.org/tutorials/basic/tutorial4b/files/1FKO_trunc_sus.pdb) Now we are able to load the pdb file into tleap.
First, we start tleap just like what we did in the previous section:
tleap -f oldff/leaprc.ff99SB
>source leaprc.gaff
>loadamberparams sustiva.frcmod
Now we load the Sustiva library file (sus.lib), followed by the complex pdb file 1FKO_trunc_sus.pdb.
>loadoff sus.lib
>complex = loadpdb [1FKO_trunc_sus.pdb](https://ambermd.org/tutorials/basic/tutorial4b/files/1FKO_trunc_sus.pdb)
Finally, we are ready to create our topology and coordinate files of the truncated RT-sustiva complex. 



>saveamberparm complex 1FKO_sus.prmtop 1FKO_sus.rst7


>savepdb complex 1FKO_sus.pdb


>quit




You can take a look at the truncated RT-Sustiva complex structure ( [1FKO_sus.pdb](https://ambermd.org/tutorials/basic/tutorial4b/files/1FKO_sus.pdb)) in VMD.


![RT-Sustiva Complex](https://ambermd.org/tutorials/basic/tutorial4b/figures/RT.png)


Again we can create a tleap input file ( [tleap2.in](https://ambermd.org/tutorials/basic/tutorial4b/files/tleap2.in)) and generate all the files for the Sustiva-RT complex:
tleap -f tleap2.in
 ## Minimize and Equilibrate the Sustiva-RT complex
 Once we have the topology and coordinate files of RT-Sustiva complex, we are ready to run short GB simulations on Sustiva-RT. Note, the procedure given here is very short in order to make the simulations compatible with the timescale of this tutorial. In a real "production" simulation you would typically run a much longer simulation (ns) in order to obtain good statistical convergence.
First, we will minimize our complex to remove any possible bad contacts. Here's our input file [min.in](https://ambermd.org/tutorials/basic/tutorial4b/files/min.in). We will do a total of 200 steps of minimization (MAXCYC) with the first 50 being steepest descent (NCYC), the remainder will be conjugate gradient (MAXCYC-NCYC). We will use a reasonably large cut off of 16 Angstroms since this is not going to be a periodic simulation and we want to deal with our electrostatics accurately (NTB=0,CUT=16). For our implicit solvent model we will use the GB model of Hawkins, Cramer and Truhlar, see the AMBER manual for a full description (IGB=1):


[min.in](https://ambermd.org/tutorials/basic/tutorial4b/files/min.in)
Initial minimisation of sustiva-RT complex


 &cntrl


 imin=1, maxcyc=200, ncyc=50,


 cut=16, ntb=0, igb=1,


 &end


Let's run our minimization: 

sander -O -i min.in -o 1FKO_sus_min.out -p 1FKO_sus.prmtop-c 1FKO_sus.rst7 -r 1FKO_sus_min.ncrst &


Now we wait for it to run. (Takes about 3 minutes on my local machine). If you can't wait the output files are here: [1FKO_sus_min.out,](https://ambermd.org/tutorials/basic/tutorial4b/files/1FKO_sus_min.out) [1FKO_sus_min.ncrst](https://ambermd.org/tutorials/basic/tutorial4b/files/1FKO_sus_min.ncrst).
You can use ambpdb to generate a pdb of the final minimized structures if you want:
 ambpdb -p 1FKO_sus.prmtop -c 1FKO_sus_min.ncrst > [1FKO_sus_min.pdb](https://ambermd.org/tutorials/basic/tutorial4b/files/1FKO_sus_min.pdb)
It is always a good idea to visualize your results since the human eye is very good at spotting anomalies.
 The next step if to heat the RT-Sustiva complex. For speed we will do this very rapidly over 1ps. Ideally you should do this for much longer.
 Here's our input file, we will run MD (imin=0) and this is not a restart (irest=0). In this example we will not use shake since it is possible that the hydrogen motion may effect the binding energy (probably not, but it serves as an example here). As we are not using shake we will need a time step smaller than normal. Here I will use a time step of 1 fs and run for 1000 steps [2 ps] (dt = 0.001, nstlim=1000, ntc=1). We will also write to our output file every 20 steps and to our trajectory [mdcrd] file every 20 steps (ntpr=20,ntwx=20). For temperature control we will use a Langevin dynamics approach with a collision frequency of 1 ps^-1. We will start our system at 0K and we want a target temperature of 300K (ntt=3, gamma_ln=1.0, tempi=0.0, temp0=300.0). And here's the input file: eq.in
[eq.in](https://ambermd.org/tutorials/basic/tutorial4b/files/eq.in)


Initial MD equilibration


 &cntrl


 imin=0, irest=0,


 nstlim=1000,dt=0.001, ntc=1,


 ntpr=20, ntwx=20,


 cut=16, ntb=0, igb=1,


 ntt=3, gamma_ln=1.0,


 tempi=0.0, temp0=300.0,


 &end


Now we run:


sander -O -i eq.in -o 1FKO_sus_eq.out-p 1FKO_sus.prmtop-c 1FKO_sus_min.ncrst -r 1FKO_sus_eq.rst -x 1FKO_sus_eq.nc &
The heating trajectory and restart coordinates are saved in [1FKO_sus_eq.mdcrd](https://ambermd.org/tutorials/basic/tutorial4b/files/1FKO_sus_eq.mdcrd) and [1FKO_sus_eq.rst](https://ambermd.org/tutorials/basic/tutorial4b/files/1FKO_sus_eq.rst), respectively. Now that the sustiva-RT complex is minimized and heated, you can take a look at the snapshot at 300K. This structure can be used as the starting point for further equilibration.




 ## References
 1Wang, J., Wolf, R.M., Caldwell, J.W., Kollman, P.A., Case, D.A. "Development and Testing of a General Amber Force Field", J. Comp. Chem., 2004, **25**, 1157 - 1173.
 2Bayly, C.I., Cieplak, P., Cornell, W.D., Kollman, P.A. "A Well-Behaved Electrostatic Potential Based Method Using Charge Restraints for Deriving Atomic Charges : The RESP Model", J. Phys. Chem, 1993, 10269-10280.
 3Jakalian, A., Bush, B.L., Jack, B.D., Bayly, C.I., "Fast, Efficient Generation of High-Quality Atomic Charges. AM1-BCC Model: I. Method.", J. Comp. Chem., 2000, **21**, 132-146.

