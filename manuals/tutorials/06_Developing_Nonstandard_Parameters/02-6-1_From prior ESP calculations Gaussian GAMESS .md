

# RESP Charge Derivation with PyPE_RESP
 This tutorial is compatible with AmberTools25 or later versions.
 ## Learning Outcomes
 - Understand the role of PyPE_RESP as an automated pipeline for RESP charge fitting.
- Learn how to write a PyPE_RESP input file to define molecules and run a complete RESP fit with a single command.
- Learn how to define charge constraint groups using atom indices and SMARTS patterns.
- Learn how to set up multi-molecule RESP fits with intermolecular constraints and atom equivalencies.
- Understand how to analyze PyPE_RESP output files: charge summaries, mol2 files, and run logs.
- Appreciate the difference between the manual pyresp_gen.py/py_resp.py workflow and the PyPE_RESP pipeline.
 ## Introduction
 In this tutorial we will use PyPE_RESP to derive RESP charges from pre-computed quantum mechanical (QM) electrostatic potential (ESP) data. PyPE_RESP is a command-line pipeline that wraps the pyresp_gen.py and py_resp.py programs shipped with AmberTools, automating the entire two-stage RESP fitting procedure through a single, human-readable input file.
 Traditional RESP charge fitting with py_resp.py requires the user to (1) manually run `espgen` to convert QM output, (2) invoke `pyresp_gen.py` to build first- and second-stage input files, (3) execute `py_resp.py` twice, and (4) extract the final charges. PyPE_RESP replaces this multi-step scripting with a single command:
 ```
pype-resp.py -i my_input.in

```
 PyPE_RESP reads Gaussian output files directly, handles format conversion internally, builds the PyRESP inputs, runs both RESP stages, and produces mol2 files with fitted charges and a charge summary report. It also supports advanced features such as SMARTS-based charge constraints, intermolecular constraint groups, and cross-molecule atom equivalencies.
 We will use capped alanine (ACE-ALA-NME) as our example molecule throughout this tutorial. This is a common system in AMBER force field development: the acetyl (ACE) and N-methyl (NME) capping groups are used to mimic the peptide backbone environment, and their charges must be constrained to zero so that the central alanine residue carries the correct total charge.
 For background on the RESP algorithm itself, see the [RESP tutorial](https://ambermd.org/tutorials/advanced/tutorial1/index.php) and the [PyRESP tutorial](https://ambermd.org/tutorials/basic/tutorial19/index.php). The [PyPE_RESP paper](https://doi.org/10.1021/acs.jcim.5c00041) describes the program in detail.
 ## Prerequisites
 This tutorial assumes you have already performed QM geometry optimization and ESP calculations using Gaussian. If you are unfamiliar with these steps, refer to steps 1–2 of the [PyRESP tutorial](https://ambermd.org/tutorials/basic/tutorial19/index.php). The Gaussian route line for ESP generation should include:
 ```
#HF/6-31G* Pop=MK IOp(6/33=2,6/42=6,6/43=20)

```
 **Note:** It is advisable to run the geometry optimization and the ESP calculation as *separate* Gaussian jobs rather than combining them into a single input. PyPE_RESP can currently struggle to parse the ESP data correctly from a combined optimization-plus-ESP output; running the ESP calculation on its own (starting from the already-optimized geometry) avoids this. This limitation will be addressed in a future release.
 If you would like PyPE_RESP to handle the QM calculations for you using the QUICK quantum chemistry program, see the companion tutorial: [Automated RESP Charge Derivation with PyPE_RESP and QUICK](https://ambermd.org/tutorials/basic/tutorial23/index.php).
 ## Getting Help
 PyPE_RESP provides two useful flags for quick reference without running a calculation:
 ```
pype-resp.py -ih           # Print detailed input file documentation
pype-resp.py -t            # Generate a template input file (PyPE_RESP_template.in)

```
 The `-ih` flag displays a complete description of every keyword. The `-t` flag writes a commented template file to the current directory that you can copy and edit.
 ## Process
 1. ### Setting Up the Working Directory
 Create a working directory and place your Gaussian output files in it. For this tutorial we use two conformers of ACE-ALA-NME. Each Gaussian output file must contain the ESP data (produced by the `Pop=MK IOp(6/33=2,...)` keywords).
 ```
mkdir ala_resp
cp ace_ala_nme_conf1.log ace_ala_nme_conf2.log ala_resp/
cd ala_resp

```
 PyPE_RESP identifies molecules by pattern-matching filenames. In this case both files contain the substring `ace_ala_nme`, so a single `MOL_1 = ace_ala_nme` definition will capture both conformers automatically.
 **Where the conformers come from:** the two ESP calculations start from distinct conformations of ACE-ALA-NME. To generate a small set of diverse starting geometries from a single structure, we used a minimal RDKit-based helper script, [generate_conformers.py](https://ambermd.org/tutorials/basic/tutorial22/include/generate_conformers.py) (it takes a PDB and writes RMSD-pruned conformer PDBs; it requires [RDKit](https://www.rdkit.org/) to be installed). Those conformers were then optimized and submitted to Gaussian for the ESP calculation. The QM steps are beyond the scope of this tutorial, but the starting structures ( [ace_ala_nme.pdb](https://ambermd.org/tutorials/basic/tutorial22/include/ace_ala_nme.pdb), [ace_gly_nme.pdb](https://ambermd.org/tutorials/basic/tutorial22/include/ace_gly_nme.pdb)) and the resulting Gaussian outputs are provided with the downloadable examples below.
2. ### Writing a Minimal Input File
 Create a file called `ala_dipeptide.in` with the following content:
 ```
# PyPE_RESP input for ACE-ALA-NME RESP charge fitting
# Two conformers from Gaussian output

MOL_1 = ace_ala_nme

FORMAT    = gaussian

```
 This is the simplest possible PyPE_RESP input. It defines one molecule ( `MOL_1`) whose files match the pattern `ace_ala_nme` and specifies Gaussian output format. Every other option takes its default value, including the two we deliberately leave out here:
 - `WORK_DIR` defaults to the directory in which `pype-resp.py` is executed, so there is no need to set it when you run from inside the working directory (as we do above).
- `ATOM_TYPE` defaults to GAFF2, so it only needs to be specified when you want a different atom type set.
 The remaining defaults give a two-stage RESP fit with standard restraint weights (0.0005 for stage 1, 0.001 for stage 2).
 Run PyPE_RESP:
 ```
pype-resp.py -i ala_dipeptide.in

```
 **Supported formats:** The `FORMAT` keyword accepts `gaussian`, `gamess`, `xyz`, and `mol2`. When using `gaussian` or `gamess`, PyPE_RESP reads the ESP data directly from the QM output files. When using `xyz` or `mol2`, no ESP data is present in the coordinate files and PyPE_RESP will use the QUICK quantum chemistry program to calculate the ESP (see the [companion QUICK tutorial](https://ambermd.org/tutorials/basic/tutorial23/index.php)).
3. ### Understanding the Output
 After a successful run, PyPE_RESP creates the following directory structure:
 ```
ala_resp/                                  # your working directory
  ├── ace_ala_nme_conf1.log                # Gaussian ESP outputs (your QM inputs)
  ├── ace_ala_nme_conf2.log
  ├── _PyPE_RESP.log                       # full run log
  └── _PyPE_RESP/                          # PyPE_RESP working directory
      ├── ESP_FILES/                       # ESP data extracted for PyRESP
      │   ├── ace_ala_nme_conf1.esp
      │   └── ace_ala_nme_conf2.esp
      ├── RESP_RUN_FILES/                  # generated PyRESP inputs + raw outputs
      │   ├── First_stage.in               # generated 1st-stage PyRESP input
      │   ├── First_stage.out              # 1st-stage PyRESP output
      │   ├── First_stage.chg              # 1st-stage fitted charges
      │   ├── First_stage.espout           # 1st-stage ESP comparison
      │   ├── Second_stage.in              # generated 2nd-stage PyRESP input
      │   ├── Second_stage.out             # 2nd-stage PyRESP output (fit statistics)
      │   ├── Second_stage.chg             # final fitted charges
      │   ├── Second_stage.espout          # 2nd-stage ESP comparison
      │   └── Concatenated.esp             # combined ESP data from all conformers
      ├── Mol3_files/
      │   └── ace_ala_nme_conf1_AFTER_RESP.mol2   # final charges, mol2 form
      └── Run_statistics/
          ├── RESP_charge_summary.log      # fit/charge summary (text)
          └── RESP_charge_summary.csv      # fit/charge summary (CSV)

```
 All generated files are collected in the `_PyPE_RESP/` working directory (the leading underscore marks it as a PyPE_RESP-managed folder), alongside the top-level `_PyPE_RESP.log`. The mol2 file is named after the first conformer, and the two summary files take the name set by `OUTPUT_NAME` (defaulting to `RESP_charge_summary`; in the complete example below they become `ala_charges.log`/ `.csv`). The key output files are:
 - **Mol3_files/<mol>_AFTER_RESP.mol2** — the main deliverable: a mol2 file carrying the final *per-atom* RESP charges, ready for use with LEaP or other force field tools.
- **RESP_RUN_FILES/Second_stage.chg** — the same final per-atom charges in PyRESP's native format, together with the equivalencing ( `ivary`) information.
- **Run_statistics/<output_name>.log / .csv** — a *summary* of the fit, not the per-atom charges: the total molecular charge, total dipole moment, the fit-quality statistics, and the net charge of each constraint group.
- **RESP_RUN_FILES/Second_stage.out** — the raw 2nd-stage PyRESP output, including the detailed fitting-statistics section.
- **ESP_FILES/** — the per-conformer ESP data extracted from the QM output.
- **_PyPE_RESP.log** — complete log of the run, including parsed options, ESP point counts, and fitting statistics.
 We will look at the actual contents of these files in the [complete single-molecule example](#complete-single) below, which is the run provided for download.
4. ### Defining Charge Constraints with Atom Indices
 When parameterizing a capped amino acid like ACE-ALA-NME, the capping groups should each sum to zero so that the central alanine residue carries the full molecular charge. This ensures the derived charges are transferable to a peptide chain where the caps are absent.
 Constraints are defined with the `MOL_X_CON_Y` keyword, where `X` is the molecule index and `Y` is the constraint group index. Atom indices are 1-based and space-separated. The target charge follows a colon; if omitted it defaults to 0.0.
 For ACE-ALA-NME with the atom ordering CH3–CO–NH–CH(CH3)–CO–NH–CH3, the ACE cap (atoms 1–6: 3H, C, C, O) and NME cap (atoms 17–22: N, H, C, 3H) can be constrained as follows:
 ```
# ACE cap (CH3-CO) constrained to 0.0
MOL_1_CON_1 = 1 2 3 4 5 6 : 0.0

# NME cap (NH-CH3) constrained to 0.0
MOL_1_CON_2 = 17 18 19 20 21 22 : 0.0

```
 If no charge is specified, the default is 0.0:
 ```
# Equivalent to MOL_1_CON_1 = 1 2 3 4 5 6 : 0.0
MOL_1_CON_1 = 1 2 3 4 5 6

```
5. ### Defining Charge Constraints with SMARTS Patterns
 For complex molecules, specifying atom indices by hand can be tedious and error-prone. PyPE_RESP supports SMARTS patterns as an alternative. The SMARTS pattern must match a unique set of atoms in the molecule (i.e., it must be unambiguous).
 For ACE-ALA-NME, the ACE methyl group and the NME methyl group can be distinguished by their chemical context using SMARTS:
 ```
# ACE methyl: CH3 attached to a carbonyl carbon
MOL_1_CON_1 = C(H)(H)(H)C(=O)N : 0.0

# NME methyl: CH3 attached to the NME nitrogen
MOL_1_CON_2 = C(H)(H)(H)|NC(=O)| : 0.0

```
 **SMARTS slicing:** When a SMARTS pattern alone would match multiple locations, you can add context atoms to make the match unique, then exclude those context atoms from the constraint using pipe ( `|`) delimiters. Atoms between the pipes are used for matching but are not included in the constraint group.
 ```
# Match hydroxyl attached to aliphatic C attached to aromatic C,
# but only constrain the HO part (exclude the Cc context):
MOL_1_CON_1 = HO|Cc|

```
 In this example, `HO|Cc|` matches the full substructure H–O–C(aliphatic)–C(aromatic), but only the H and O atoms end up in the constraint group. This is particularly useful for distinguishing between otherwise identical functional groups in different chemical environments.
 **Note:** Each constraint group ( `MOL_X_CON_Y`) is interpreted as *either* atom indices or a SMARTS pattern, depending on its content—so a single group cannot itself mix the two styles. You can, however, freely combine index-based and SMARTS-based constraint groups within the same molecule (e.g. `MOL_1_CON_1` by indices and `MOL_1_CON_2` by SMARTS).
 **Testing SMARTS patterns:** Before running a full calculation, you can verify that your SMARTS patterns match the intended atoms using the `TEST_SMARTS` keyword:
 ```
# Add to your input file:
TEST_SMARTS = True

```
 This will report which atoms each SMARTS pattern matches and exit without performing the RESP fit.
6. ### A Complete Single-Molecule Example
 Putting it all together, here is a complete input file for ACE-ALA-NME with cap constraints. This is the run provided for download (in the `ala_resp` working directory):
 ```
# PyPE_RESP input: ACE-ALA-NME with capping group constraints
# Two conformers from Gaussian ESP calculations

MOL_1 = ace_ala_nme

FORMAT    = gaussian

# Constrain ACE cap (atoms 1-6) to neutral
MOL_1_CON_1 = 1 2 3 4 5 6 : 0.0

# Constrain NME cap (atoms 17-22) to neutral
MOL_1_CON_2 = 17 18 19 20 21 22 : 0.0

OUTPUT_NAME = ala_charges
VERBOSE     = True

```
 ```
pype-resp.py -i ala_dipeptide.in

```
 **Downloads:** the complete worked run — input, both Gaussian outputs, and the full `_PyPE_RESP/` tree — is packaged as [ala_resp.tar.gz](https://ambermd.org/tutorials/basic/tutorial22/include/ala_resp.tar.gz). The key files are also available individually: input [ala_dipeptide.in](https://ambermd.org/tutorials/basic/tutorial22/include/ala_resp/ala_dipeptide.in); final charges [ace_ala_nme_conf1_AFTER_RESP.mol2](https://ambermd.org/tutorials/basic/tutorial22/include/ala_resp/ace_ala_nme_conf1_AFTER_RESP.mol2) and [Second_stage.chg](https://ambermd.org/tutorials/basic/tutorial22/include/ala_resp/Second_stage.chg); fit output [Second_stage.out](https://ambermd.org/tutorials/basic/tutorial22/include/ala_resp/Second_stage.out); and summary [ala_charges.log](https://ambermd.org/tutorials/basic/tutorial22/include/ala_resp/ala_charges.log) / [ala_charges.csv](https://ambermd.org/tutorials/basic/tutorial22/include/ala_resp/ala_charges.csv).
 **The fitted charges.** The final per-atom charges are written into the mol2 file (and into `Second_stage.chg`). Here are the first six atoms — the ACE cap — from `ace_ala_nme_conf1_AFTER_RESP.mol2`:
 ```
@<TRIPOS>ATOM
      1 H        -4.00324    0.08443    0.18500 hc     1 MOL      0.073671865
      2 C        -3.32052   -0.75157    0.44462 c3     1 MOL     -0.174801630
      3 H        -3.22139   -0.81241    1.54918 hc     1 MOL      0.073671865
      4 H        -3.75594   -1.70112    0.06799 hc     1 MOL      0.073671865
      5 C        -1.97787   -0.50910   -0.16273 c      1 MOL      0.479125900
      6 O        -1.45853   -1.38941   -0.90167 o      1 MOL     -0.525339860

```
 These six charges sum to 0.0 — exactly the ACE-cap constraint we requested — so the cap is electrostatically neutral and the parameters are transferable to a full peptide where the cap is absent.
 **The summary file.** Rather than per-atom charges, `ala_charges.log` reports the fit quality and the net charge of each constraint group, confirming both caps came out neutral:
 ```
 Fitting Statistics Summary
 ...
 The relative RMSE           (RRMSE = sqrt(RSS/ssvpot))         1.0424322E-01

CONSTRAINED SUMMARY
...
ace_ala_nme_conf1:
The sum of charges over all atoms of ace_ala_nme_conf1 is: 9e-09
The total dipole moment (Debye) is 8.06302 (X: 1.50293 Y: 1.79172 Z: -7.71642)

Constraints:
1. INTRAMOLECULAR
    Ids: 1 2 3 4 5 6
    Names: (H C H H C O)
    Charge: 0.000000005
2. INTRAMOLECULAR
    Ids: 17 18 19 20 21 22
    Names: (N H C H H H)
    Charge: 0.000000010

```
 The relative RMSE of about 0.10 indicates a good fit to the QM ESP, the total molecular charge is zero to within rounding, and each cap constraint group is satisfied to ~10−8 e. The remaining alanine-residue atoms (7–16) carry the full molecular charge.
7. ### Multi-Molecule RESP Fitting
 PyPE_RESP can fit charges for multiple molecules simultaneously. This is essential when parameterizing molecules that share a chemical interface, such as polymer building blocks or host-guest systems. Each molecule is defined with its own `MOL_X` keyword:
 ```
# Two molecules, identified by filename substrings
MOL_1 = ace_ala_nme
MOL_2 = ace_gly_nme

FORMAT    = gaussian

```
 PyPE_RESP will search the working directory for all Gaussian output files whose names contain `ace_ala_nme` and assign them to molecule 1 (treating multiple matches as conformers), and likewise for `ace_gly_nme` and molecule 2.
8. ### Intermolecular Charge Constraints
 When fitting multiple molecules, you may need to constrain atoms across molecule boundaries to a combined charge. For example, when parameterizing two peptide building blocks simultaneously, the backbone atoms at the junction should sum to zero across both molecules to ensure charge neutrality at the interface.
 First, define intramolecular constraint groups for each molecule, then link them with an `INTER_CON_GROUP` definition:
 ```
# Define the C-terminal backbone atoms in ALA
MOL_1_CON_1 = 13 14 15 16

# Define the N-terminal backbone atoms in GLY
MOL_2_CON_1 = 7 8 9

# Link them: backbone atoms at the junction must sum to 0.0
INTER_CON_GROUP_1 = (MOL_1_CON_1 MOL_2_CON_1) : 0.0

```
 When a constraint group is part of an intermolecular definition, any charge specified on the individual `MOL_X_CON_Y` line is ignored — the charge on the `INTER_CON_GROUP` line takes precedence. Constraint groups that are *not* referenced by an `INTER_CON_GROUP` remain intramolecular constraints.
 A maximum of two constraint groups can be combined in a single `INTER_CON_GROUP` definition.
9. ### Atom Equivalencies Across Molecules
 By default, PyPE_RESP automatically determines conformers of each molecule and equivalences corresponding atoms across conformers. However, if atoms in *different* molecules should receive identical charges, use the `EQUATE` keyword:
 ```
# The carbonyl oxygen in ALA must equal the carbonyl oxygen in GLY
EQUATE = MOL_1 14 ; MOL_2 8

```
 Each `EQUATE` statement specifies one atom per molecule, separated by semicolons. You can have multiple `EQUATE` lines for different equivalence relationships. SMARTS patterns are also supported:
 ```
# The amide nitrogen in MOL_1 equals the amide nitrogen in MOL_2
EQUATE = MOL_1 [NH] ; MOL_2 [NH]

```
 **Note:** Unlike constraints, atom indices and SMARTS patterns are *not* mutually exclusive for equivalencies — you can mix both styles in the same input file.
10. ### A Complete Multi-Molecule Example
 Below is a complete input file demonstrating most PyPE_RESP features for a simultaneous fit of capped alanine and capped glycine:
 ```
# PyPE_RESP input: simultaneous ALA + GLY parameterization
# Each molecule has two conformers in the working directory

MOL_1 = ace_ala_nme
MOL_2 = ace_gly_nme

FORMAT    = gaussian

# -- Cap constraints for both molecules --
# ACE caps constrained to 0
MOL_1_CON_1 = 1 2 3 4 5 6 : 0.0
MOL_2_CON_1 = 1 2 3 4 5 6 : 0.0

# NME caps constrained to 0
MOL_1_CON_2 = 17 18 19 20 21 22 : 0.0
MOL_2_CON_2 = 14 15 16 17 18 19 : 0.0

# -- C-terminal backbone atoms for intermolecular constraint --
MOL_1_CON_3 = 13 14 15 16
MOL_2_CON_3 = 10 11 12 13

# -- Intermolecular constraint: backbone junction sums to 0 --
INTER_CON_GROUP_1 = (MOL_1_CON_3 MOL_2_CON_3) : 0.0

# -- Cross-molecule equivalency: carbonyl oxygens --
EQUATE = MOL_1 14 ; MOL_2 11

# -- Optional settings --
OUTPUT_NAME = dipeptide_charges
VERBOSE     = True

```
 **Running it.** Set up a separate working directory containing the Gaussian outputs for *both* molecules — here two conformers of ACE-ALA-NME and five of ACE-GLY-NME — save the input above as `mult_mol_resp.in`, and run a single command:
 ```
mkdir mult_mol
cp ace_ala_nme_conf*.log ace_gly_nme_conf*.log mult_mol/
cd mult_mol
pype-resp.py -i mult_mol_resp.in

```
 PyPE_RESP groups the conformers by filename pattern, fits both molecules simultaneously, and writes one mol2 per molecule ( `ace_ala_nme_conf1_AFTER_RESP.mol2` and `ace_gly_nme_conf1_AFTER_RESP.mol2`).
 **Downloads:** the complete run is packaged as [mult_mol.tar.gz](https://ambermd.org/tutorials/basic/tutorial22/include/mult_mol.tar.gz). Individually: input [mult_mol_resp.in](https://ambermd.org/tutorials/basic/tutorial22/include/mult_mol/mult_mol_resp.in); final charges [ace_ala_nme_conf1_AFTER_RESP.mol2](https://ambermd.org/tutorials/basic/tutorial22/include/mult_mol/ace_ala_nme_conf1_AFTER_RESP.mol2), [ace_gly_nme_conf1_AFTER_RESP.mol2](https://ambermd.org/tutorials/basic/tutorial22/include/mult_mol/ace_gly_nme_conf1_AFTER_RESP.mol2), [Second_stage.chg](https://ambermd.org/tutorials/basic/tutorial22/include/mult_mol/Second_stage.chg); fit output [Second_stage.out](https://ambermd.org/tutorials/basic/tutorial22/include/mult_mol/Second_stage.out); and summary [dipeptide_charges.log](https://ambermd.org/tutorials/basic/tutorial22/include/mult_mol/dipeptide_charges.log) / [dipeptide_charges.csv](https://ambermd.org/tutorials/basic/tutorial22/include/mult_mol/dipeptide_charges.csv).
 **Reading the summary.** Because this run combines intramolecular cap constraints, an intermolecular junction constraint, and a cross-molecule equivalency, the `dipeptide_charges.log` summary is the clearest way to confirm everything was applied. Its constrained summary lists, for each molecule, every constraint group and its net charge — including the intermolecular link between the two backbone groups:
 ```
ace_ala_nme_conf1:
...
Constraints:
1. INTRAMOLECULAR
    Ids: 1 2 3 4 5 6
    Names: (H C H H C O)
    Charge: -0.000000004
2. INTRAMOLECULAR
    Ids: 17 18 19 20 21 22
    Names: (N H C H H H)
    Charge: 0.000000010
3. INTERMOLECULAR with 10 11 12 13 (H H C O) of ace_gly_nme_conf1
    Ids: 13 14 15 16
    Names: (H H C O)
    Charge: 0.000000002
...
ace_gly_nme_conf1:
...
3. INTERMOLECULAR (see ace_ala_nme_conf1 above)
    Ids: 10 11 12 13

```
 Group 3 of each molecule is reported as a single `INTERMOLECULAR` constraint whose combined charge sums to zero across the junction, while the four cap groups are each neutral — exactly the constraint scheme defined in the input.
11. ### Additional Options
 PyPE_RESP exposes several options for controlling the RESP fitting procedure:
 KeywordDefaultDescription `CHARGE_MODEL`chgElectrostatic model: `chg` (point charge/RESP), `ind` (induced dipole/pGM-ind), or `perm` (permanent dipole/pGM-perm) `POLAR_FILE`NonePath to atomic polarizability file; required if CHARGE_MODEL is `ind` or `perm` `CHARGE_WEIGHT_1`0.0005Charge restraint weight for 1st RESP stage `CHARGE_WEIGHT_2`0.001Charge restraint weight for 2nd RESP stage `FIRST_STAGE_ONLY`FalseRun only the 1st RESP stage (skip 2nd stage) `SEARCH_DEPTH`3Maximum graph traversal depth for automatic equivalent atom detection `MOL_NAME`MOLMolecule name inserted into mol2 files (max 3 characters) `NCPUS`1Number of CPUs for parallel execution `NORESP`FalsePrepare all PyRESP input files but do not run the fit (useful for manual inspection) `MOL3_ONLY`FalseRegenerate the final charge-containing mol2 files from a previous run's output (e.g. after a `NORESP` run) without repeating the RESP fit For a complete list of all options, run `pype-resp.py -ih` or generate a template with `pype-resp.py -t`.
12. ### Comparison with the Manual Workflow
 To illustrate the convenience of PyPE_RESP, consider the steps required to fit RESP charges for ACE-ALA-NME with two conformers using py_resp.py directly (as shown in the [PyRESP tutorial](https://ambermd.org/tutorials/basic/tutorial19/index.php)):
 ```
# Manual workflow:
espgen -i ace_ala_nme_conf1.log -o conf1.dat -p 1
espgen -i ace_ala_nme_conf2.log -o conf2.dat -p 1
cat conf1.dat conf2.dat > combined.dat
pyresp_gen.py -i combined.dat -f1 1st.in -f2 2nd.in -p chg -q 0

#  --> MANUAL STEP: open 1st.in and 2nd.in in a text editor and hand-edit
#      the charge-constraint section to neutralize the ACE and NME caps
#      (correct atom indices, charge group counts, and the relevant flags)

py_resp.py -O -i 1st.in -o 1st.out -e combined.dat -s 1st.esp -t 1st.chg
py_resp.py -O -i 2nd.in -o 2nd.out -e combined.dat -s 2nd.esp -t 2nd.chg -q 1st.chg

#  --> MANUAL STEP: read the final charges out of 2nd.out (or 2nd.chg) by hand
#      and assemble them into a mol2/charge file for downstream use

```
 The two `MANUAL STEP` comments above are the easy-to-overlook part. `pyresp_gen.py` only produces an *unconstrained* pair of input files: to enforce the cap constraints you must open `1st.in` and `2nd.in` and edit the charge-constraint block by hand—getting the atom indices, group counts, and control flags exactly right, and keeping them consistent between the two stages. And `py_resp.py` writes its charges into a verbose text output, so you then have to locate and transcribe the final per-atom charges yourself before they can be used. Both steps are manual, error-prone, and have to be repeated for every molecule and every change.
 The equivalent PyPE_RESP input file and command:
 ```
# ala_dipeptide.in
MOL_1     = ace_ala_nme
FORMAT    = gaussian

```
 ```
pype-resp.py -i ala_dipeptide.in

```
 PyPE_RESP handles format conversion, concatenation, input generation, both fitting stages, and post-processing in a single invocation. Crucially, it also absorbs the two manual steps above: constraints are declared once with `MOL_X_CON_Y` keywords and written into *both* stage inputs automatically, and the fitted charges are extracted and written straight into ready-to-use mol2 files and a charge summary—no hand-editing of input files and no manual transcription of charges. For multi-molecule fits with constraints and equivalencies, the manual approach becomes significantly more complex and error-prone while the PyPE_RESP input file remains straightforward.
 ## Next Steps
 - If you do not have pre-computed QM data and want PyPE_RESP to handle the ESP calculation, see [Automated RESP Charge Derivation with PyPE_RESP and QUICK](https://ambermd.org/tutorials/basic/tutorial23/index.php).
- For background on the RESP algorithm and the pyresp_gen.py/py_resp.py programs, see the [PyRESP tutorial](https://ambermd.org/tutorials/basic/tutorial19/index.php).
- For details on manually running multi-conformational RESP fits, see [Setting up a DNA-Ligand System](https://ambermd.org/tutorials/advanced/tutorial1/index.php).
 by Marco Lapsien

