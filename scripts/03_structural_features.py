import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, ShrakeRupley
from Bio.Data.IUPACData import protein_letters_3to1

PDB_FILE = "structures/2F9Q.pdb"
CHAIN = "A"   # confirmed in PyMOL on Day 2
OFFSET = 0    # PDB residue number minus UniProt P10635 number

# Theoretical maximum ASA per amino acid (Tien et al., 2013, PLoS ONE)
MAX_ASA = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167,
    "Q": 225, "E": 223, "G": 104, "H": 224, "I": 197,
    "L": 201, "K": 236, "M": 224, "F": 240, "P": 159,
    "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}

model = PDBParser(QUIET=True).get_structure("cyp2d6", PDB_FILE)[0]

# Remove crystallographic water molecules (hetflag == "W")
for ch in model:
    for rid in [r.id for r in ch if r.id[0] == "W"]:
        ch.detach_child(rid)

chain = model[CHAIN]

# Compute residue-level SASA on the single Chain A monomer (+ its heme)
ShrakeRupley().compute(chain, level="R")

# Locate the catalytic iron (FE) atom of the heme (HEM) cofactor in Chain A
heme_fe = next(r["FE"] for r in chain if r.get_resname().strip() == "HEM")

rows = []
for _, v in pd.read_csv("data/processed/variants_mappable.csv").iterrows():
    d = v.to_dict()
    # Standard Biopython residue ID format: (" ", resseq, " ")
    key = (" ", int(v["position"]) + OFFSET, " ")
    if key not in chain:
        rows.append({**d, "in_structure": False})
        continue

    res = chain[key]
    wt = protein_letters_3to1.get(res.get_resname().strip().capitalize(), "?")
    dist = min(np.linalg.norm(a.coord - heme_fe.coord) for a in res)

    rows.append({
        **d,
        "in_structure": True,
        "pdb_residue": wt,
        "wt_matches": wt == v["wt"],
        "min_dist_heme_Fe_A": round(float(dist), 2),
        "rel_sasa": round(res.sasa / MAX_ASA.get(wt, 200), 2),
    })

out = pd.DataFrame(rows)
out.to_csv("data/processed/variants_structural.csv", index=False)

# Print summary and mandatory QC check
print(f"Saved data/processed/variants_structural.csv ({len(out)} rows)")
in_struct = out[out["in_structure"] == True]
missing = out[out["in_structure"] == False]
print(f"In structure: {len(in_struct)} | Unresolved in 2F9Q Chain A: {len(missing)}")

mismatches = in_struct[in_struct["wt_matches"] == False]
if len(mismatches) > 0:
    print("WARNING - wt_matches is False for:", mismatches[["allele", "protein_change", "pdb_residue"]].to_dict("records"))
else:
    print("QC PASSED: wt_matches is True for every residue present in 2F9Q Chain A!")

print()
cols = ["allele", "protein_change", "function", "in_structure", "pdb_residue", "wt_matches", "min_dist_heme_Fe_A", "rel_sasa"]
print(out[cols].to_string(index=False))