# AmberTools Workflow Patterns

## Standard MD Pipeline

**When to use**: Any protein/DNA/RNA simulation in explicit solvent.

**How**:
```
pdb4amber --in raw.pdb --out clean.pdb
tleap -f leaprc.protein.ff19SB
  > source leaprc.water.opc
  > mol = loadpdb clean.pdb
  > solvateoct mol OPCBOX 10.0
  > addions2 mol Na+ 0
  > saveamberparm mol prmtop inpcrd

# Minimization: imin=1, ncyc=500, maxcyc=1000, ntr=1, restraint_wt=10.0
# Heating: imin=0, ntt=3, temp0=300, ntc=2, ntf=2, dt=0.002, nstlim=25000
# Equilibration: ntp=1, ntb=2, nstlim=50000
# Production: nstlim=5000000, ioutfm=1, ntwx=5000, ntwr=5000

cpptraj -p prmtop << EOF
trajin prod.nc
autoimage
rms first :1-300@CA
strip :WAT,Na+,Cl-
trajout stripped.nc
go
EOF
```

**Trade-offs**: OPC water recommended for FF19SB. Cubic box is simplest; truncated octahedron saves ~25% water. 4 fs timestep requires HMR.

---

## Explicit Water Solvation

**When to use**: Any simulation needing atomistic solvent.

**How**:
```
solvateoct mol OPCBOX 10.0   # 10 A buffer, truncated octahedron
solvatebox mol TIP3PBOX 12.0  # 12 A buffer, cubic box
addions2 mol Na+ 0            # Neutralize (0 = auto)
addions2 mol Cl- 0            # Neutralize
addions2 mol Na+ Cl- 50.0     # Add 50 mM NaCl
```

**Trade-offs**: Truncated octahedron preferred for globular proteins. Minimum buffer: 8 A for PME, 10 A safer. TIP3P is fastest but less accurate than OPC.

---

## Implicit Solvent (GB)

**When to use**: Fast sampling, large systems, or when explicit water is too expensive.

**How**:
```
 &cntrl igb=8, saltcon=0.1, ntb=0, cut=999.0, /
```
igb=1 (GB^OBC I), igb=2 (GB^OBC II), igb=5 (GB^OBC II + mbondi2), igb=7 (GB^neck2), igb=8 (GB^neck2 + mbondi3).

**Trade-offs**: igb=8 with mbondi3 radii and FF19SB is current best. igb=5 is standard for FF14SB. Much faster than explicit water but less accurate for solvent-mediated interactions.

---

## Small Molecule Parameterization (antechamber + GAFF)

**When to use**: Any non-standard small molecule (ligand, cofactor, drug).

**How**:
```
antechamber -i lig.mol2 -fi mol2 -o lig.mol2 -fo mol2 -c bcc -s 2 -nc 0
parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod -s 2
```
In LEaP: `lig = loadmol2 lig.mol2` then `loadamberparams lig.frcmod`.

**Trade-offs**: AM1-BCC is fast but less accurate than RESP. RESP requires Gaussian03 or equivalent. For charged ligands, RESP recommended.

---

## Modified Residue

**When to use**: Non-standard amino acids, post-translational modifications, fluorescent dyes.

**How**: Build residue, parameterize like a small molecule, create mol2/prep file, load in LEaP with `loadamberprep` or `loadmol2`.

**Trade-offs**: Requires careful charge derivation. For modified amino acids, consider ff15ipq-m force field.

---

## HMR (Hydrogen Mass Repartitioning)

**When to use**: Enable 4 fs timestep in any production MD.

**How**:
```
parmed -p original.parm7 << EOF
hmassrepartition
outparm hmr.parm7
go
EOF
```
Run with `dt=0.004`, `ntc=2`, `ntf=2`.

**Trade-offs**: 2x speedup for PME simulations. Requires SHAKE. Works with all force fields. Not recommended for QM/MM or very high temperature.

---

## Protein-Ligand Binding Free Energy (MM-PBSA)

**When to use**: Rank-order binding affinity estimates for congeneric series.

**How**:
```
MMPBSA.py -O -i mmpbsa.in -o FINAL_RESULTS.dat \
  -sp complex.prmtop -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop \
  -y complex_traj.nc
```
mmpbsa.in: `&general startframe=500, endframe=1000, interval=5, / &gb igb=5, saltcon=0.15, /`

**Trade-offs**: 1-trajectory protocol is faster but may miss conformational changes. 3-trajectory is more rigorous. GB is faster than PB. Use 100+ frames for convergence.

---

## Alchemical Free Energy (TI)

**When to use**: Rigorous relative binding free energies between similar ligands.

**How**:
```
&cntrl icfe=1, ifsc=1, scalpha=0.5, scbeta=12.0, clambda=0.0,
        timask1=':LIG', timask2=':LIG', /
```
Run 12-24 lambda windows (0.0 to 1.0), analyze with BAR or MBAR.

**Trade-offs**: Softcore potentials essential (ifsc=1). scalpha=0.5, scbeta=12.0 are standard. ACES for enhanced sampling. Expensive but most accurate.

---

## Umbrella Sampling PMF

**When to use**: Free energy along a reaction coordinate (conformational change, binding path).

**How**: Generate windows along coordinate, restrain with harmonic potential, analyze with WHAM.

**Trade-offs**: Requires careful window spacing. NFE toolkit provides automated umbrella sampling via pmd module.

---

## Membrane Protein Setup

**When to use**: Any membrane-bound protein simulation.

**How**:
```
# PACKMOL-Memgen
packmol-memgen --pdb protein.pdb --lipids POPC --ratio 1 --salt --saltconc 0.15
```
Or use CHARMM-GUI Membrane Builder, then convert with `charmmlipid2amber.py`.

**Trade-offs**: PACKMOL-Memgen is simpler. CHARMM-GUI gives more control. Lipid21 is the recommended force field.

---

## Metal Ion Parameterization (Bonded Model)

**When to use**: Catalytic metal sites with strong coordination.

**How**:
```
MCPB.py -i input.in -s 1  # Step 1: identify metal site
MCPB.py -i input.in -s 2  # Step 2: QM calculation
MCPB.py -i input.in -s 3  # Step 3: force constant derivation
MCPB.py -i input.in -s 4  # Step 4: tleap input generation
```

**Trade-offs**: Bonded model captures specific coordination geometry. Nonbonded 12-6-4 model is simpler for structural ions. pyMSMT provides alternative.

---

## Constant pH MD

**When to use**: Systems where protonation states matter (enzyme mechanism, pH-dependent binding).

**How**:
```
&cntrl icnstph=1, ntcnstph=100, solvph=7.0, /
```
Use `pdb4amber --constantph` for residue naming.

**Trade-offs**: CpHMD captures protonation dynamics. Requires GB implicit solvent (explicit solvent CpHMD under development). C(E)pHMD for continuous pH.

---

## REMD for Enhanced Sampling

**When to use**: Overcome energy barriers, sample conformational space.

**How**:
```
mpirun -np 32 pmemd.cuda.MPI -rem 3 -remlog rem.log -ng 32 -groupfile groupfile
```
Typical: 32-64 replicas spanning 300-450 K with geometric temperature distribution.

**Trade-offs**: More replicas = better exchange but more expensive. GPU T-REMD is efficient. H-REMD for Hamiltonian exchange.

---

## Multi-GPU Parallel MD

**When to use**: Large systems (>500K atoms) or long timescales (>1 us).

**How**:
```
mpirun -np 4 pmemd.cuda.MPI -O -i mdin -o mdout -p prmtop -c inpcrd
```

**Trade-offs**: CUDA_VISIBLE_DEVICES for GPU assignment. Near-linear scaling to 4-8 GPUs for typical systems. Use pmemd.cuda.MPI, not sander.MPI for production.

---

## QM/MM Simulation

**When to use**: Bond breaking/forming, electronic transitions, metal reactivity.

**How**:
```
&cntrl ifqnt=1, /
&qmmm qmmask=':1-30', qmcharge=0, qm_theory='DFTB3',
       qmcut=8.0, writepdb=1, /
```
Supported QM: SCC-DFTB, PM6, AM1, PM3, DFTB3, external (Gaussian, ORCA, TeraChem, QUICK, xTB, DFTB+).

**Trade-offs**: sander only. No GPU for QM portion. Adaptive QM/MM for exchanging solvent. DPRc for ML corrections.

---

## GIST Water Analysis

**When to use**: Identify favorable/displaceable water sites in binding pockets.

**How**:
```
cpptraj -p prmtop << EOF
trajin prod.nc 1 last 10
gist doorder doeij gridcntr 10 15 5 griddim 40 40 40 \
     gridspacn 0.5 refdens 0.0334 temp 300.0 prefix gist
go
EOF
```

**Trade-offs**: Requires 10-20 ns restrained-solute simulation. gridspacn=0.5 A for detail, 0.75 A for faster convergence. Strip ions before GIST.

---

## FEW Automated Free Energy

**When to use**: Batch processing of multiple ligands against one receptor.

**How**: Prepare command files, run `perl FEW.pl MMPBSA command_file` for MD setup, then `perl FEW.pl MMPBSA analysis_file` for free energy analysis.

**Trade-offs**: Automates the entire pipeline. Uses ff12SB by default. Requires Perl and specific modules. Only cubic water boxes supported.

---

## Lipid21 Membrane Simulation

**When to use**: Simulations with biological membranes.

**How**:
```
source leaprc.lipid21
source leaprc.protein.ff19SB
source leaprc.water.opc
```
Build with PACKMOL-Memgen or CHARMM-GUI. Use `ntp=2,3` for semi-isotropic pressure coupling (membrane plane vs. normal).

**Trade-offs**: Lipid21 is current recommended. ntp=2 for semi-isotropic membrane pressure. Use ntp=3 for surface tension control. Membrane area per lipid should be monitored.