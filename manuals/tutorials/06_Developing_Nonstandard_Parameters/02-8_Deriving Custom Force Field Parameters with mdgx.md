

# Deriving bonded parameters with `mdgx`: Introduction
 Parameter development need not be a black box. The key is to have tools that speak in terms of coordinates and energies, recognize that the molecular mechanics energy is merely a relative quantity, treat the physical phases of matter correctly, and understand that the Newtonian parameters are prone to over-fitting. There's some finesse to the process, but machines can take away a huge amount of the tedium and handle much bigger problems (in sheer terms of the number of parameters and size of the data sets) than a human being fitting one parameter at a time. This is the right way to approach the problem: parameters do not exist in isolation‐not at the level of approximation we are taking wth molecular mechanics. Moreover, a sophisticated apparatus for collecting and interpreting the data can be enhanced by a rich set of post-processing tools to interpret the resulting parameters, place them in the context of the data, probe their vulnerabilities to prevent over-fitting, and ensure that the models will function in the intended conditions. This is what `mdgx` aims to offer.
 [Stage 1: The molecule in many, many guises](https://ambermd.org/tutorials/advanced/tutorial32/Part1.php)


 [Stage 2: Fitting parameters](https://ambermd.org/tutorials/advanced/tutorial32/Part2.php)


 [Stage 3: Iterate](https://ambermd.org/tutorials/advanced/tutorial32/Part3.php)


 [Stage 4: Use the new parameters](https://ambermd.org/tutorials/advanced/tutorial32/Part4.php)
 Proper credit should be given to `ParamFit` and its creators Robin Betz and Ross Walker. A number of the features of the `mdgx` force field tools are inspired by that project, particularly the methods for generating conformations of the molecules of interest and dealing with data in batches. In building the `mdgx` force field tools, I have tried to extend and refine a lot of the process in [the ParamFit tutorial](http://ambermd.org/tutorials/advanced/tutorial23/). I hope that the following tutorial will convince you that `mdgx` sets a high standard in this field, and that the wisdom of numerous experts has been distilled into the tools.
 This bonded parameter development is ideally used after fitting charges with the [IPolQ scheme](http://ambermd.org/tutorials/advanced/tutorial28/index.php), but if you are happy with the charges you've got it can also be taken standalone. Like the IPolQ that preceded it, this is an iterative process, but in this case the iterations are meant to relieve human operators from the burden of having to intercede to eliminate spurious parameters. We feel that the automated approach is in fact better for this purpose in the way that it fixes what can be widespread and subtle errors: humans cannot easily identify these errors, much less correct them, but new parameter sets are rife with them. The automated tuning is also just as effective at taking care of the obvious catastrophic problems that human intuition and meticulous care have solved in the past.
 **The parameter development scheme is detailed in these papers:**


 - D.S. Cerutti, J.E. Rice, W.C. Swope, and D.A. Case. (2013) "Derivation of Fixed Partial Charges for Amino Acids Accommodating a Specific Water Model and Implicit Polarization." *J. Phys. Chem. B* 117: 2328-2338. [link](http://pubs.acs.org/doi/suppl/10.1021/jp311851r)


 - D.S. Cerutti, W.C. Swope, J.E. Rice, and D.A. Case. (2014) "ff14ipq: A Self-Consistent Force Field for Condensed-Phase Simulations of Proteins." *J. Chem. Theory Comput.* 10: 4515-4534. [link](http://pubs.acs.org/doi/abs/10.1021/ct500643c)


 - K.T. Debiec, D.S. Cerutti, L.R. Baker, A.M. Gronenborn, D.A. Case, and L.T. Chong. (2016) "Further along the Road Less Traveled: AMBER ff15ipq, an Original Protein Force Field Built on a Self-Consistent Physical Model." *J. Chem. Theory Comput.* 12: 3926-3947. [link](http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00567)


