# Chapter 25: ProPrep -- Protein Preparation

## Core Commands & Syntax

```bash
proprep [options]
```

### Input Options
| Option | Description |
|--------|-------------|
| `-pdbid ID` | Download structure by PDB ID, create `{ID}_results/` |
| `-pdbfile PATH` | Load local PDB file, create `{filename}_results/` |
| `-output-dir DIR` | Specify output directory instead of being prompted |

### Interface Mode Options
| Option | Description |
|--------|-------------|
| `-workflow` | Start in guided workflow mode |
| `-full-menu` | Start in full menu mode |
| `-set-default {workflow\|full-menu\|ask}` | Set default interface mode |

### Session Recording Options
| Option | Description |
|--------|-------------|
| `-no-auto-record` (`-no-session`) | Disable automatic session recording |
| `-session-dir DIR` | Directory for session files |
| `-session-description TEXT` | Session metadata description |
| `-resume-session FILE` | Resume a previous session |
| `-force-replay FILE` | Replay session without prompting |
| `-demo-delay SECONDS` | Delay between replayed actions (default 0.0) |

### Template and Batch Processing
| Option | Description |
|--------|-------------|
| `-create-template SESSION_FILE` | Convert session to reusable template |
| `-validate-template TEMPLATE_FILE` | Validate template without executing |
| `-template-info TEMPLATE_FILE` | Display template variables and values |
| `-batch-replay TEMPLATE_FILE` | Batch process with template |
| `-batch-list INPUT_FILE` | File of identifiers and variable substitutions |
| `-batch-continue` | Continue batch after entry failure |

### Shortcut Options
| Option | Description |
|--------|-------------|
| `-analysis` | Jump directly to simulation analysis browser |
| `-pdbview ID_OR_FILE` | Launch interactive structure viewer then exit |

### General Options
| Option | Description |
|--------|-------------|
| `-debug` | Enable debug mode with full stack traces |
| `-verbose` (`-v`) | Print detailed internal logging |
| `-version` | Show version information |

## Interface Modes

### Guided Workflow Mode (5 sequential stages)
```
Load -> Detect -> Compare -> Fix -> Simulate
```

**Load:** Structure Loader, Biological Assembly Generator
**Detect:** Redox Site Detector, Forcefield Explorer, Forcefield Parameterizer
**Compare:** Homology Searcher, Structure Aligner
**Fix:** PDB Filter, Amino Acid Mutator, Structure Fixer, Redox Site Preparer, Protonation State Analyzer, Membrane Builder, Structure Orientator
**Simulate:** TLEaP Topology Generator, Molecular Dynamics Manager, ONIOM QM/MM Preparator, ORCA QM/MM Preparator

Navigation: `n` (next), `b` (previous), or stage keys: `l` (Load), `d` (Detect), `c` (Compare), `f` (Fix), `s` (Simulate). `v` opens structure viewer, `u` opens utilities, `a` switches to full menu.

### Full Menu Mode
All modules exposed simultaneously in numbered sections. Select by alphanumeric key (e.g., `1a` for Structure Loader, `4c` for Structure Fixer).

## Key Modules

### Structure Loader
Four sources: RCSB PDB (by ID or keyword search), AlphaFold DB (by UniProt ID or gene name), AlphaFill DB (enriched predictions), local PDB files. Supports batch download, biological assembly download, mmCIF loading (inspection only -- redox requires PDB format). Pre-existing AMBER files (.prmtop, .parm7, .rst7, .inpcrd) can be loaded directly.

### Redox Site Detector
Six-phase interactive detection:
1. **Configuration:** select metal elements, non-metal elements, search radius (default 4.0 A), S-S distance threshold (default 2.5 A)
2. **Inventory scan:** organometallic cofactors, isolated metal ions, redox-active amino acids (TYR, TRP, PHE, MET), disulfide bonds
3. **Selection and grouping:** group by chain/residue name/number; assign site type labels
4. **Site type classification:** labels enable template-based automation
5. **Site refinement:** distance cutoff search (fixed or adaptive radius) or residue count cutoff search; interactive bond definition
6. **Review, finalize, and export** to JSON

### Forcefield Explorer
Interactive browser for AMBER forcefield parameters. Features: overview statistics, browse atom types, browse residues, five parameter search modes (bond, angle, dihedral, VDW/nonbonded, mass), reverse lookup, structure atom typing with automatic terminal detection and protonation state resolution.

### Forcefield Parameterizer
Six-level residue classification cascade:
1. User-defined, 2. Redox site, 3. Metal ion, 4. Modified amino acid, 5. Small molecule, 6. Unknown

**Small molecule workflow:** Gaussian input generation (B3LYP/6-31+G(d) opt, HF/6-31G(d) ESP) -> antechamber GAFF2 typing + RESP -> parmchk2 -> optional Seminario method refinement.

**Modified amino acid workflow:** ACE-XXX-NME tripeptide -> Gaussian -> RESP with backbone restraints -> residuegen (.lib) + parmchk2.

**Metal site workflow:** Preprocessing (9 steps) -> atom typing/model building (small model + large model) -> Seminario bonded parameters -> RESP charge fitting -> mol2/frcmod generation -> forcefield integration.

### PDB Filter
Chain interface analysis (buried surface area heatmap, topology analysis). Component selection with redox site awareness. Water analysis: metal proximity, H-bond analysis, B-factor, burial, interface proximity, multi-radius burial profiling, directional burial, water network analysis.

### Amino Acid Mutator
Standard (`A:123:ALA->GLY`) and non-standard (`A:154:PHE->CNF`) mutations. Non-standard mutations specify atom retention (backbone only, backbone+CB, all, custom). MODELLER required for standard mutations.

### Structure Fixer
Detects missing residues (REMARK 465, SEQRES vs. coordinates, FASTA comparison, sequence gaps), missing atoms (Chemical Component Dictionary, REMARK 470), alternate conformations. Fill missing segments with MODELLER.

### Interactive Structure Viewer
Browser-based NGL Viewer (WebGL 3D). Local HTTP server on port 8765 (auto-fallback). 11 representation styles, 6 color schemes, distance and angle measurements, view controls (center, reset, fullscreen, screenshot). Non-blocking.

## Session Recording and Templates
Every interaction recorded to timestamped JSON session files. Session editor commands: `e N` (edit), `d N` (delete), `u N` (undo), `v N` (template variable), `save`, `quit`. Templates use named variables for batch processing. Batch input supports plain text (one identifier per line) or CSV with headers.

## Optional Dependencies
- **MODELLER:** required for structure repair and mutagenesis (free for academic use)
- **Gaussian:** required for QM calculations in parameterization and QM/MM workflows (commercial)

## Key Takeaways

1. **ProPrep** is an interactive protein preparation framework with 5-stage guided workflow or full menu mode.
2. **Session recording and replay** enable reproducible workflows; templates support batch processing of multiple structures.
3. **Redox site detection** is fully interactive with 6-phase workflow; site types enable template-based automation.
4. **Forcefield Parameterizer** handles small molecules (GAFF2), modified amino acids (capped peptide RESP), and metal sites (MCPB workflow).
5. **All state is workspace-based**; modules auto-detect prerequisites and display status indicators.

## Connects To

- Chapter 2: LEaP (topology generation)
- Chapter 21: antechamber (GAFF atom typing)
- Chapter 23: sqm (semi-empirical QM)
- Chapter 27: MCPB (metal center parameterization)
- Chapter 31: External library interface (external QM programs)