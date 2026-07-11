

# Constant pH MD Example
 ### Calculating pKas for titratable side chains in HEWL
 ![](https://ambermd.org/tutorials/advanced/tutorial18/images/main_pic.jpg)
 In this tutorial, we will learn how to use software in Amber and AmberTools to carry out molecular dynamics simulations at constant pH (CpHMD) on the hen egg white lysozyme (HEWL). HEWL is a common benchmark for simulations at constant pH. The CpHMD method was implemented in sander by John Mongan, and is described in the corresponding journal article. [[1]](https://livecomsjournal.org/index.php/livecoms/article/view/v4i1e1563/1389) This tutorial will assume familiarity with Ref. [[1]](https://livecomsjournal.org/index.php/livecoms/article/view/v4i1e1563/1389), and will outline only the practical steps for running CpHMD simulations.
 The CpHMD method described by Mongan, et al. works only in Generalized Born implicit solvent. This method is available in *sander*, *pmemd* and *pmemd.cuda*.
 This tutorial consists of four sections:
 > 1) [section1](https://ambermd.org/tutorials/advanced/tutorial18/section1.php): Creating the initial structure and input files.
>  2) [section2](https://ambermd.org/tutorials/advanced/tutorial18/section2.php): Preparing the system by relaxing bad initial contacts, followed by heating and equilibration stages.
>  3) [section3](https://ambermd.org/tutorials/advanced/tutorial18/section3.php): Running the production simulations at different pH environments.
>  4) [section4](https://ambermd.org/tutorials/advanced/tutorial18/section4.php): Analyzing the results.
  [CLICK HERE TO GO TO SECTION 1](https://ambermd.org/tutorials/advanced/tutorial18/section1.php)
  References
 **[1]** John Mongan, David A. Case, and J. Andrew McCammon, "Constant pH molecular dynamics in generalized Born implicit solvent", *J. Comput. Chem.*, **2004**, 25 (16), pp. 2038-2048
 by Jason Swails & T. Dwight McGee Jr 2013


