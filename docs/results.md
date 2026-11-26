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

## Day 11: Quantitative Structural & Conservation Summary (2026-11-26)

Quantitative metrics were summarized across the **12 missense variants resolved in PDB `2F9Q` Chain A** (`figures/summary_structural_metrics.png`), alongside South Asian allele frequency context contrasting 3D-mappable and non-mappable star alleles (`figures/summary_sas_frequencies.png`). Three missense variants—`V11M` (`*35`, normal), `R28C` (`*22`, unknown), and `G42R` (`*12`, no function)—fall in regions unresolved in `2F9Q` Chain A and were excluded from 3D coordinate calculations while retaining their sequence-conservation scores.

### Summary Metrics by CPIC Clinical Function Class (`2F9Q` Chain A, n = 12)

| CPIC Function Class | Mapped n | Heme Fe Distance Median [Range] (A) | Relative SASA Median [Range] | CYP2D-Family Entropy Median [Range] (bits) |
| :--- | :---: | :---: | :---: | :---: |
| No function | 2 | 18.18 [11.31–25.04] | 0.32 [0.24–0.40] | 0.681 [0.530–0.832] |
| Decreased function | 6 | 22.43 [5.91–32.52] | 0.24 [0.18–0.38] | 0.439 [0.000–2.084] |
| Normal function | 3 | 15.75 [15.10–25.69] | 0.21 [0.19–0.30] | 1.509 [0.957–1.986] |
| Unknown / uncertain (`*25` R343G) | 1 | 26.52 [26.52–26.52] | 0.32 [0.32–0.32] | 0.439 [0.439–0.439] |

### Key Observations
- **Evolutionary constraint separates clinical classes more cleanly than 3D distance alone:** Normal-function variants occupy more variable positions across the mammalian CYP2D family (median entropy `1.509` bits), whereas decreased-function (`0.439` bits, including a 100% conserved position at `0.000` bits) and no-function (`0.681` bits) variants sit at more conserved sites.
- **Impaired function arises via multiple structural mechanisms, not only heme contact:** While active-site (`5.91 A`) and proximal heme-loop (`11.31 A`, `R440H`) variants lie near the cofactor, other loss- and reduced-function variants lie >20 A away (e.g., `H324P` at `25.04 A`), where helical disruption or folding instability rather than direct heme contact drives functional loss. Accordingly, median heme distance does not separate clinical classes monotonically.
- **Solvent accessibility overlaps across classes:** Median `rel_sasa` values are closely clustered (`0.21–0.32`), reinforcing that local stereochemistry (such as proline insertion or proximal-loop electrostatics) and sequence conservation provide stronger mechanistic context than simple burial depth.
- **Uncertain-function variant note (`*25` `R343G`):** Although distal to the heme (`26.52 A`, `rel_sasa = 0.32`), `R343G` occurs at a conserved CYP2D-family position (`0.439` bits). This is noted solely as a structural observation and does not alter its CPIC uncertain-function classification.
- **Population & scope takeaway:** In the SG10K_Health Indian cohort (`N = 328`), major clinically actionable reduced- and no-function alleles include non-mappable splicing and deletion variants—`*41` (`12.20%`), `*4` (`6.40%`), and `*5` (`6.10%`)—confirming why structural missense atlases must be paired with explicit documentation of non-mappable alleles (`data/processed/variants_not_mappable.csv`).
- **Statistical caveat:** Sample sizes per CPIC class (`n = 1` to `6`) are intentionally small in this curated atlas; findings are reported descriptively using medians and ranges without underpowered hypothesis testing.