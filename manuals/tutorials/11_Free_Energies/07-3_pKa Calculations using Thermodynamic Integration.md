

(Note: These tutorials are meant to provide illustrative examples of how to use the AMBER software suite to carry out simulations that can be run on a simple workstation in a reasonable period of time. They do not necessarily provide the optimal choice of parameters or methods for the particular application area.)


 Copyright Ross Walker 2006
 **pKa Calculations using Thermodynamic Integration**
 **By [Ross Walker](mailto:amber_tutorial_query@rosswalker.co.uk) & Mike Crowley**
 ![](https://ambermd.org/tutorials/advanced/tutorial6/images/energy_plot.gif)
 This tutorial reproduces the calculation of the pKa value of the ASP residue in the protein thioredoxin as described in the following paper:
 > Simonson, T., Carlsson, J., Case, D.A., "Proton Binding to Proteins: pKa Calculations with Explicit and Implicit Solvent Models", JACS 2004, 126, pp4167-4180.
 You should make sure that you read this paper before attempting this tutorial. If you are in the Amber 2006 advanced workshop this was provided as part of the pre-requisite reading material.
 The authors of the paper carried out a number of simulations using both implicit and explicit solvent and with the AMBER and CHARMM force fields. For the purposes of this tutorial we will just go over setting up the implicit solvent simulation for a single aspartate residue (ASP 26) in the protein thioredoxin. You are welcome to try the explicit solvent case of your own accord and/or the other aspartate residues or a different protein.
 The calculation of the pKa value is based on the following equation:
 ![](https://ambermd.org/tutorials/advanced/tutorial6/images/pka_equation1.gif)
 The value of pKa,model is the experimental pKa for a model system. In this case ACE/NME blocked aspartate. We are actually interested in calculating the pKa shift of the aspartate residue when part of the protein matrix. In this case we don't need to know the experimental pKa of the model system. We are interested in obtaining the value of the second half of the equation above involving delta-delta G. This represents the pKa shift and is what we will calculate. We will do this using the thermodynamic cycle shown in the image at the top of this page.
 To do this we need to use thermodynamic integration to calculate the deltaG for the conversion of our model AspH to Asp- and the deltaG for the conversion of thioredoxin-AspH to thioredoxin-Asp-. We will start by doing the calculation of deltaG for the model compound.
 This tutorial consists of a total of 4 steps as follows:
 1. [section1](https://ambermd.org/tutorials/advanced/tutorial6/section1.php) : Setup the prmtop and inpcrd files for the AspH and Asp- models.
1. [section2](https://ambermd.org/tutorials/advanced/tutorial6/section2.php) : Run thermodynamic integration on the model compound to find delta-G(model).
2. [section3](https://ambermd.org/tutorials/advanced/tutorial6/section3.php) : Setup the prmtop and inpcrd files for the AspH and Asp- thioredoxin protein.
3. [section4](https://ambermd.org/tutorials/advanced/tutorial6/section4.php) : Run thermodynamic integration on the thioredoxin protein. Calculate deltadeltaG and pKa shift.
  [CLICK HERE TO GO TO SECTION 1](https://ambermd.org/tutorials/advanced/tutorial6/section1.php)


