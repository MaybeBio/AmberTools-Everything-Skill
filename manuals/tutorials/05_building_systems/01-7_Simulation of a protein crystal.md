

# Simulation of a Protein Crystal
Formerly known as "Simulating Crystals with the AMBER Molecular Dynamics Package"
By David Cerutti, Dave Case; last updated Oct. 15, 2015
**Note:** This tutorial is written for users who are already familiar with the basic operation of *LEaP, sander* and *pmemd.* The simulation of molecular crystals is in principle no different than the simulation of other periodic systems, but setting up the simulations requires some auxiliary programs and iterative procedures that most solution-phase simulations do not.
 Nearly all simulations of biomolecules begin with an X-ray or NMR structure of the system of interest, a starting point for which the coordinates of most atoms in the system are known with certainty. Nearly all simulations of biomolecules are also done in solution, with the biomolecular system immersed in a bath that approximates water. For direct comparison between crystallographic or NMR structures and molecular simulations, the conditions of the original experiment must be reproduced in the computer model. The solution conditions of NMR structures are fairly straightforward to reproduce, and can often be approximated by a bath of pure water or even an implicit solvent such as Generalized Born. The environment of a biomolecule in a protein crystal lattice is somewhat more difficult to model: in addition to the numerous symmetry-related copies of the biomolecule itself, there typically exists a heterogeneous solvent environment which is mostly aqueous but may also contain precipitating agents or cryoprotectants. This tutorial is geared for advanced users of the AMBER biomolecular simulations package, provides somewhat depthy narrative on the challenge of adequately modeling the crystalline environment, and incorporates several auxiliary programs to aid in the preparation of crystal simulations.
  ## Summary of Protocol
 1. [Constructing a lattice: You'll need some more tools for this!](https://ambermd.org/tutorials/advanced/tutorial13/Tarball.php)
[](https://ambermd.org/tutorials/advanced/tutorial13/Tarball.php) [](https://ambermd.org/tutorials/advanced/tutorial13/Selection.php)
2. [System selection: Not all crystal structures are equal!](https://ambermd.org/tutorials/advanced/tutorial13/Selection.php)
[](https://ambermd.org/tutorials/advanced/tutorial13/Selection.php) [](https://ambermd.org/tutorials/advanced/tutorial13/Construction.php)
3. [Constructing the unit cell: Where crystal simulations depart from ordinary MD](https://ambermd.org/tutorials/advanced/tutorial13/Construction.php)
[](https://ambermd.org/tutorials/advanced/tutorial13/Construction.php) [](https://ambermd.org/tutorials/advanced/tutorial13/Solvation.php)
4. [Solvating the unit cell: Filling in the gaps](https://ambermd.org/tutorials/advanced/tutorial13/Solvation.php)
[](https://ambermd.org/tutorials/advanced/tutorial13/Solvation.php) [](https://ambermd.org/tutorials/advanced/tutorial13/Topology.php)
5. [Generating a topology: Working with what we have](https://ambermd.org/tutorials/advanced/tutorial13/Topology.php)
[](https://ambermd.org/tutorials/advanced/tutorial13/Topology.php) [](https://ambermd.org/tutorials/advanced/tutorial13/References.php)
6. [Suggested reading on crystallographic simulations](https://ambermd.org/tutorials/advanced/tutorial13/References.php)
[](https://ambermd.org/tutorials/advanced/tutorial13/References.php)

