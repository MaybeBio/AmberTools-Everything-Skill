

# Building a Peptide Sequence
 

 

 ![chignolin](https://ambermd.org/tutorials/basic/tutorial10/include/chignolin.png)
 Implementation of methods outlined in: 

 Nguyen, H., Maier, J., Huang, H., Perrone, V. and C. Simmerling, Folding Simulations for Proteins with Diverse Topologies are Accessible in Days with Physics-Based Force Field and Implicit Solvent, [*J. Am. Chem. Soc.*, **2014,** *136(40),* 13959–13962.](https://pubs.acs.org/doi/full/10.1021/ja5032776)
 ## Learning Outcomes
 - Load an appropriate force field
- Specify radii for an implicit solvent GB calculation.
- Utilize the *sequence* command within LEaP to build chignolin, a 10 amino acid synthetic peptide.
- Save a pdb file in tleap.
- Save a prmtop and inpcrd file tleap.
 ## Overview
 If you want to start a simulation from an extended structure and not based upon an experimentally-determined structure, then the *sequence* command in LEaP will make this starting structure for you. This short tutorial shows you how to build a sequence, choose the correct force field and save an inputs for sander or pmemd. Note that if you want to run implicit solvent simulations, you will want to specify your radii before saving your topology and input coordinate files. This tutorial uses the 10 amino acid peptide [chignolin](https://www.cell.com/structure/fulltext/S0969-2126(04)00242-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0969212604002424%3Fshowall%3Dtrue).
 ## Process
 1. Start tLEaP with the "tleap" command ```
>$ tleap-I: Adding /amber20/dat/leap/prep to search path.-I: Adding /amber20/dat/leap/lib to search path.-I: Adding /amber20/dat/leap/parm to search path.-I: Adding /amber20/dat/leap/cmd to search path.Welcome to LEaP!>(no leaprc in search path)
```
 2. Load an appropriate force field. 

 

 A MD force field is the potential energy function that describes the intermolecular and intramolecular forces governing the molecules in the system. Chignolin is a small peptide, made of amino acids. Therefore, a protein force field will need to be chosen to represent the atoms in chignolin. Each type of molecule (protein, DNA, RNA, carbohydrate, lipid, water, ions, co-factors) requires a different set of parameters. To learn more about force fields, please look through section **3.1.1** on page **36** of the [Amber 2020 Reference Manual](https://ambermd.org/doc12/Amber20.pdf#page=36).
In Amber, the standard is that there is both an ***lib*** file and a ***dat*** or ***parm*** file that describes each class of molecule. Lib files can be found under:
 ```
$AMBERHOME/dat/leap/lib
```
 Dat and parm files can be found under:
 ```
$AMBERHOME/dat/leap/parm
```
 Some premade ***leaprc*** files which load a set of lib and parm files for different biomolecules can be found under:
 ```
$AMBERHOME/dat/leap/cmd
```
 In Nguygen *et al.,* the ff14SBonlysc force field is used. This will be used to describe the potential energy of chignolin. ```
> source leaprc.protein.ff14SBonlysc
----- Source: /amber20/dat/leap/cmd/leaprc.protein.ff14SBonlysc
----- Source of /amber20/dat/leap/cmd/leaprc.protein.ff14SBonlysc done
Log file: ./leap.log
Loading parameters: /amber20/dat/leap/parm/parm10.dat
Reading title:
PARM99 + frcmod.ff99SB + frcmod.parmbsc0 + OL3 for RNA
Loading parameters: /amber20/dat/leap/parm/frcmod.ff14SB
Reading force field modification type file (frcmod)
Reading title:
ff14SB protein backbone and sidechain parameters
Loading parameters: /amber20/dat/leap/parm/frcmod.ff99SB14
Reading force field modification type file (frcmod)
Reading title:
ff99SB backbone parameters (Hornak & Simmerling) with ff14SB atom types
Loading library: /amber20/dat/leap/lib/amino12.lib
Loading library: /amber20/dat/leap/lib/aminoct12.lib
Loading library: /amber20/dat/leap/lib/aminont12.lib


```
 3. Use GBneck2 radii for implicit solvent calculation. Implicit solvent calculations use an electric field to mimic the solvent and the Generalized Born model is a common implimentation. For a thorough discussion of the theory and different parameters, see [Cramer and Truhlar.](https://link.springer.com/chapter/10.1007/0-306-46931-6_1) Nguyen *et al.* employed the [GBneck2](https://pubs.acs.org/doi/10.1021/ct3010485") parameters.
```
set default PBradii mbondi3
```
 You should get the tleap message
 ```
Using ArgH and AspGluO modified Bondi2 radii
```
 

 4. Build the chignolin molecule with the "sequence" command 

 

 Once we have loaded the appropriate force field into tLEaP, the "sequence" command allows us to connect amino acid "units" together to build a protein, stored here in the **"chig"** unit. Leap will automatically recognize that GLY1 is the first residues and should be the N-terminus with an NH3+ and that GLY10 should be the C-terminal end with a carboxylate. To learn about other ways to build proteins, please look through section **13.4.2** on pages **222** and **223** of the [Amber 2020 Reference Manual](https://ambermd.org/doc12/Amber20.pdf#page=222). ```
> chig = sequence { GLY TYR ASP PRO GLU THR GLY THR TRP GLY }
```
 You now have an chignolin molecule stored in the unit **chig**. 

 

 5. Save the chig unit as a pdb file. 

 

 To save as a pdb file, type
```
> savepdb chig chignolin_ext.pdb
```
 You should get
 ```
Writing pdb file: chignolin_ext.pdb
```
 6. Save the chig unit as parm7 topology files and inpcrd files. A topology file and a starting set of coordinates are required to run anything with Amber (minimization, molecular dynamics, etc.). This step will save your UNIT (chig) with the appropriate force field parameters (parm7) and starting coordinates (inpcrd file). 

 

To save an Amber parm7 topology file and inpcrd file, type
```
>saveamberparm chig chig_ext.parm7 chig_ext.crd
```
 You should get
 ```
>Checking Unit.

Warning: The unperturbed charge of the unit (-2.000000) is not zero.

Note: Ignoring the warning from Unit Checking.

Building topology.
Building atom parameters.
Building bond parameters.
Building angle parameters.
Building proper torsion parameters.
Building improper torsion parameters.
 total 35 improper torsions applied
Building H-Bond parameters.
Incorporating Non-Bonded adjustments.
Not Marking per-residue atom chain types.
Marking per-residue atom chain types.
 (no restraints)

```
 In general when running *tLEap*, you should run it with the goal of no Warnings or Errors. Anytime an atom is created or a Warning or Error appears, you should review it and make sure it's what you thought should happen. If atoms were added that you don't think should be there, then you did something wrong. Once you quit *tLEap*, you can also review the leap.log file that is made. Note that the leap.log file is a concatenated log of all your *tLEap* runs in that folder. So make sure you look in the most recent run.
 Above, there's a Warning say that the overall charge is nonzero. This is fine for an implicit solvent calculation. This will not be fine for an explicit solvent calculation. [Tutorial 1.4 Building a Protein System](https://ambermd.org/tutorials/basic/tutorial7/index.php) walks you through how to add counterions, among other things.
 7. Quit tleap To quit tleap, just type 'quit'
```
> quit
	Quit

Exiting LEaP: Errors = 0; Warnings = 1; Notes = 1.


```
 You now have an extended starting set of coordinates ( [chig_ext.crd](https://ambermd.org/tutorials/basic/tutorial10/include/chig_ext.crd)) and a parm7 topology file ( [chig_ext.top](https://ambermd.org/tutorials/basic/tutorial10/include/chig_ext.parm7)) that you can use to run an implicit solvent calculation. If you don't want to use GBneck2, then skip step 3 or substitute with the preferred radii. It should be noted that because of the limits of the Ergodic hypothesis, lack of sampling, that most people in the field begin from an experimental structure.
 by Jan Ziembicki and Maria Nagan


