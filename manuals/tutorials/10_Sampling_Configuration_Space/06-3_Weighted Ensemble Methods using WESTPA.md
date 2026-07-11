

# Westpa Tutorials
 ## Learning Outcomes
 - Understand how Westpa methods can enhance sampling of rare events.
- Understand how to implement Westpa methods to sample rare events.
 ## Introduction
  Figure 1. Sampling of rare events.
 The problem with standard simulations (blue spheres in Figure 1) in sampling long-timescale processes such as those that involve large protein conformational changes is that most of the time is spent “waiting” in the initial stable state for a “lucky” transition over the free energy barrier.[1](#References) However, these transitions are relatively fast once they occur. To more efficiently simulate these interesting, functional transitions, the weighted ensemble strategy focuses the computing power on the actual transitions rather than the stable states themselves.
 By focusing computing power on the transitions between stable states (Figure 1), the weighted ensemble strategy can be orders of magnitude more efficient than standard simulations in generating unbiased, continuous pathways and rate constants for rare events (e.g., large-scale protein conformational transitions and protein binding). In addition, the strategy can be applied to any stochastic process[2,3,4](#References) and a target state does not need to be strictly defined in advanced.[5](#References)
 To run our weighted ensemble simulations, we used the [open-source WESTPA software](https://github.com/westpa/westpa).[6](#References)
 Here are set of tutorials for WESTPA developed by the WESTPA developers.
 - [Live Coms Tutorials](https://livecomsjournal.org/index.php/livecoms/article/view/v1i2e10607)
- [Protein Folding of Chignolin](https://github.com/westpa/westpa/wiki/Protein-Folding-%28Chignolin%29-with-AMBER-16)
- [Github page for WESTPA Tutorials](https://github.com/westpa/westpa/wiki/Tutorials)
  ## References
 1. Zwier, M. C. & Chong, L. T. Reaching biological timescales with all-atom molecular dynamics simulations. Current Opinion in Pharmacology vol. 10 745–752 (2010).
2. Huber, G. A. & Kim, S. Weighted-ensemble Brownian dynamics simulations for protein association reactions. Biophysical Journal vol. 70 97–110 (1996).
3. Zhang, B. W., Jasnow, D. & Zuckerman, D. M. The ‘weighted ensemble’ path sampling method is statistically exact for a broad class of stochastic processes and binning procedures. J. Chem. Phys. 132, 054107 (2010).
4. Zuckerman, D. M. & Chong, L. T. Weighted Ensemble Simulation: Review of Methodology, Applications, and Software. Annu. Rev. Biophys. 46, 43–57 (2017).
5. Suárez, E. et al. Simultaneous Computation of Dynamical and Equilibrium Information Using a Weighted Ensemble of Trajectories. J. Chem. Theory Comput. 10, 2658–2667 (2014).
6. Zwier, M. C. et al. WESTPA: an interoperable, highly scalable software package for weighted ensemble simulation and analysis. J. Chem. Theory Comput. 11, 800–809 (2015).
  By Jeremy Leung and Darian Yang



