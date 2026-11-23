# Structure Selection & Inspection Decisions

**Date:** 2026-11-22  
**Target Protein:** Human Cytochrome P450 2D6 (UniProt: `P10635`, 497 amino acids)  
**Selected Structure:** PDB `2F9Q` (Chain A)

## 1. Inspection of PDB `2F9Q`
* **Experimental method & resolution:** X-ray diffraction at 3.00 Å (Rowland et al., *J Biol Chem* 2006).
* **Chains present:** `[fill from get_chains — e.g., 4 chains: A, B, C, D]`. Chain A will be used as the single monomer for all measurements and visualisations.
* **Modeled residues:** `[fill total CA count]` CA atoms across all chains; `[fill Chain A CA count]` residues modeled in Chain A (spanning residues `[start]` to `[end]`).
* **Missing regions:**
  * N-terminal transmembrane anchor (~residues 1–`[start - 1]`), truncated for soluble expression in *E. coli*.
  * Internal unresolved loop gaps on Chain A: `[fill from PyMOL output, or "none" if empty]`.
  * C-terminal tail beyond residue `[end]` (if applicable).
* **Engineered/construct substitutions noted on RCSB:** Truncated N-terminus plus solubility mutations (`L230D`, `L231R`) and polymorphic `V374M`.
* **UniProt numbering check (`P10635`):** Checked positions 296 (`ARG`) and 486 (`SER`) in PyMOL; PDB residue numbering matches UniProt `P10635` directly (`OFFSET = 0`).
* **Active-site cofactor and ligands:** Contains catalytic heme (`HEM`, iron atom `FE` verified in Chain A at residue `[fill HEM resi]`). No drug or substrate ligand is bound in the active site (apo conformation); only crystallization sulfate ions (`SO4`) appear on Chain B.

## 2. Comparison with Other Human CYP2D6 Structures (`P10635`)
Searched RCSB PDB for alternative human CYP2D6 structures:
* **`3TBG` (thioridazine complex):** Bound to a psychiatric antipsychotic drug; useful for future ligand-contact analysis, but the pocket reflects ligand-induced conformational changes.
* **`3QM4` (prinomastat complex, 2.85 Å)** and **`4WNU` (quinidine complex, 2.26 Å):** Higher-resolution inhibitor-bound structures.

## 3. Final Decision
Proceeding with **`2F9Q` (Chain A)** as the primary reference structure for the atlas because it provides an unbiased apo conformational baseline for measuring variant distances to the catalytic heme iron (`HEM:FE`) and residue solvent exposure. Any variants falling in the missing N-terminal anchor or unresolved loop gaps will be flagged (`in_structure = False`) rather than modeled. Higher-resolution holo structures (`3TBG`, `4WNU`) are noted as natural extensions for drug-interaction or docking follow-up studies.