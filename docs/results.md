# Results & Progress Notes

## Week 1 Checkpoint (2026-11-24)
- Curated 21 CYP2D6 star-allele variants from PharmVar/CPIC (15 mappable missense substitutions and 6 non-mappable/in-frame alleles) with verified South Asian/Indian population frequencies attached for 15 rows and 6 left blank where reliable data were unavailable.
- Automated PDB (`2F9Q`) and UniProt (`P10635`) downloads (`01_fetch_data.py`) and parsed variants (`02_build_variant_table.py`), confirming 100% wild-type amino-acid agreement with UniProt `P10635`.
- Computed minimum heavy-atom distance to the catalytic heme iron (`HEM:FE`) and Shrake–Rupley relative SASA on `2F9Q` Chain A (`03_structural_features.py`); `[X]` missense variants mapped inside the crystal structure (`wt_matches = True` for all) and `[Y]` N-terminal variants were flagged `in_structure = False`.
- Calculated CYP2D-family Shannon entropy across `[N]` reviewed mammalian CYP2D sequences aligned with EBI Clustal Omega (`04_conservation.py`) and merged all metrics into `data/processed/variants_final.csv`.
- Ready for Week 2: PyMOL overview/close-up rendering, scripted figure automation, analytical summary plots, and full interpretation.

## Day 9: Local Environment & Illustrative Mutagenesis Notes (2026-11-25)

Local 5.0 A heavy-atom environments and candidate polar contacts (N/O pairs <= 3.5 A) extracted from PDB `2F9Q` Chain A for four showcase variants spanning no-function, decreased-function, and normal-function classes:

### `R440H` (`*31`, no function)
- **Quantitative metrics (`2F9Q` Chain A):** Heme Fe distance = `11.31 A` | Relative SASA = `0.4` | CYP2D-family entropy = `0.53`
- **Residues within 5.0 A:** GLY439, ARG441, ARG132, TYR124, TRP128, ALA442, ARG129, LEU91, GLU96, GLY95, ALA99, ALA438
- **Candidate polar contacts (<= 3.5 A):** O-ARG132:NH1 (2.84 A); NH1-TYR124:OH (3.24 A)
- **Structural interpretation:** Sits in the heme-binding loop just upstream of the heme-coordinating cysteine (C443) on the proximal face of the enzyme; substituting histidine shortens the basic side chain and alters local electrostatic/polar packing around the heme-binding region.

### `F120I` (`*49`, decreased function)
- **Quantitative metrics (`2F9Q` Chain A):** Heme Fe distance = `5.91 A` | Relative SASA = `0.18` | CYP2D-family entropy = `1.77`
- **Residues within 5.0 A:** LEU121, VAL119, ARG101, ASP301, ARG441, ALA122, HEM600, MET374, GLY118, ALA305
- **Candidate polar contacts (<= 3.5 A):** O-ARG101:NE (2.77 A); N-ASP301:OD1 (2.89 A); O-ARG441:NH1 (3.27 A); O-ARG101:NH1 (3.37 A)
- **Structural interpretation:** Sits on the B-C loop directly lining the hydrophobic active-site pocket over the heme; replacing the planar aromatic ring of phenylalanine with branched isoleucine alters active-site geometry and removes aromatic stacking interactions.

### `H324P` (`*7`, no function)
- **Quantitative metrics (`2F9Q` Chain A):** Heme Fe distance = `25.04 A` | Relative SASA = `0.24` | CYP2D-family entropy = `0.832`
- **Residues within 5.0 A:** LEU323, PRO325, LEU320, VAL327, ASP326, GLN328, ILE322, PHE415, MET321, ARG474, LEU319
- **Candidate polar contacts (<= 3.5 A):** N-LEU320:O (2.68 A); O-VAL327:N (2.98 A); O-ASP326:N (3.12 A); ND1-LEU323:O (3.22 A); O-GLN328:N (3.30 A); N-MET321:O (3.37 A); O-MET321:O (3.45 A)
- **Structural interpretation:** Positioned on the I-helix; introducing a rigid proline removes the backbone amide H-bond donor and induces helical distortion that destabilises local packing.

### `R296C` (`*2`, normal function)
- **Quantitative metrics (`2F9Q` Chain A):** Heme Fe distance = `15.75 A` | Relative SASA = `0.21` | CYP2D-family entropy = `0.957`
- **Residues within 5.0 A:** ILE297, LEU295, ALA300, ASN294, VAL298, GLU293, ASP292, VAL299, LEU251, LEU248, LEU255, ASP252, LEU110, ASP301
- **Candidate polar contacts (<= 3.5 A):** O-ALA300:N (2.91 A); N-ASN294:O (3.33 A)
- **Structural interpretation:** Positioned near the N-terminal end of the I-helix, distal to the catalytic heme iron; serves as a tolerated normal-function backbone control shared across multiple star alleles.

> **Illustrative Mutagenesis Caveat:** Side-by-side wild-type vs. mutant figures (`figures/closeups/mutagenesis_R440H.png`, `mutagenesis_F120I.png`, `mutagenesis_H324P.png`) were generated using PyMOL's Mutagenesis Wizard (backbone-dependent rotamer selection with minimal steric clashes). These views are strictly illustrative—no energy minimization, backbone relaxation, or folding free-energy ($\Delta\Delta G$) calculation was performed.
