# Results & Progress Notes

## Week 1 Checkpoint (2026-11-24)
- Curated 21 CYP2D6 star-allele variants from PharmVar/CPIC (15 mappable missense substitutions and 6 non-mappable/in-frame alleles) with verified South Asian/Indian population frequencies attached for 15 rows and 6 left blank where reliable data were unavailable.
- Automated PDB (`2F9Q`) and UniProt (`P10635`) downloads (`01_fetch_data.py`) and parsed variants (`02_build_variant_table.py`), confirming 100% wild-type amino-acid agreement with UniProt `P10635`.
- Computed minimum heavy-atom distance to the catalytic heme iron (`HEM:FE`) and Shrake–Rupley relative SASA on `2F9Q` Chain A (`03_structural_features.py`); `[X]` missense variants mapped inside the crystal structure (`wt_matches = True` for all) and `[Y]` N-terminal variants were flagged `in_structure = False`.
- Calculated CYP2D-family Shannon entropy across `[N]` reviewed mammalian CYP2D sequences aligned with EBI Clustal Omega (`04_conservation.py`) and merged all metrics into `data/processed/variants_final.csv`.
- Ready for Week 2: PyMOL overview/close-up rendering, scripted figure automation, analytical summary plots, and full interpretation.