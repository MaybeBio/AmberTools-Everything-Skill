Title: Overview Amber Tutorials

URL Source: https://ambermd.org/tutorials/Overview.php

Markdown Content:
Prior to molecular dynamics simulations, it can be helpful to familiarize yourself with the [tools](https://ambermd.org/tutorials/Tools.php) listed here. Many of the tutorials listed below depend on this prerequisite knowledge.

## For the impatient

All users _should_ peruse the sections below. But if you just want to quickly see how Amber works, try using this [this tutorial](https://ambermd.org/tutorials/basic/tutorial0/index.php) as a starting point.

## Introduction

If you are new to simulations, here are the general steps involved:

### [1. Building a System](https://ambermd.org/tutorials/BuildingSystems.php)

    
This is the most complicated step. Before you start a simulation you will need to decide a number of items:

    *   What types of ions/molecules do you have in your system?

    *   Do you have [a force field](https://ambermd.org/AmberModels.php) for each type of ion/molecule? If you don't [how will you obtain or derive it](https://ambermd.org/tutorials/ForceField.php)?

    *   Will you run your simlulation with explicit or implicit solvent? Which model for explicit solvent will you employ?

[](https://ambermd.org/tutorials/Relaxation.php)
### [2. Relaxing the System](https://ambermd.org/tutorials/Relaxation.php)

    
Once you have built the system, then you will want to get rid of bad contacts in the solute and if you have explicit solvent and/or ions, you will want to relax them around the solute to form a stable system.

### [3. Collecting Trajectories under Production Run Conditions](https://ambermd.org/tutorials/Production.php)

    
Amber has a number of production run MD engines including sander, sander.MPI, pmemd, pmemd.cuda, pmemd.cuda.MPI and others. For a full list and the differences [please see the manual](https://ambermd.org/Manuals.php).

### [4. Analyzing Data](https://ambermd.org/tutorials/TrajectoryAnalysis.php)

    
Once the trajectory is collected, you'll want to analyze the data. This is often the most time consuming part of a project and one that will need to be customized to your scientific question.

The Amber tutorials are organized with these items in mind in a modular fashion so that you can learn about each step as you need it.

## Case Studies

Case studies have been developed that cover a lot of these steps collectively. Some of these tutorials are meant to provide illustrative examples of how to use the AMBER software suite. These can be run on a simple workstation in a reasonable period of time. They do not necessarily provide the optimal choice of parameters or methods for the particular application area. These should not be used as a standard for excellent science but rather a guide for using the software.

Please refer to published papers for better guidance regarding scientific methods and the [Amber Reference Manual](https://ambermd.org/Manuals.php) for implementation of specific flags.

We recommend using this [this tutorial](https://ambermd.org/tutorials/basic/tutorial0/index.php) as a starting point.

For a full list of case studies, please see [this page](https://ambermd.org/tutorials/Introductory.php).

## More Advanced Methods

For information on more advanced methods, please see the specific area of the [Amber Reference Manual.](https://ambermd.org/Manuals.php)

*   [Sampling Configuration Space Tutorials](https://ambermd.org/tutorials/EnhancedSampling.php)

*   [Calculating Free Energies Tutorials](https://ambermd.org/tutorials/FreeEnergy.php)

*   [Maintaining Chemical Equilibria Tutorials](https://ambermd.org/tutorials/ChemicalReactions.php)

*   [NMR Refinement Tutorial](https://ambermd.org/tutorials/Refinement.php)




