from Bio import AlignIO
import numpy as np
import pandas as pd

aln = AlignIO.read("data/raw/cyp2d_family_aln.fasta", "fasta")
human = next(r for r in aln if "CP2D6_HUMAN" in r.id or "P10635" in r.id)


def entropy(col):
    aa = [c for c in col if c != "-"]
    if len(aa) < 0.7 * len(col):
        return np.nan  # too many gaps (>30%)
    p = pd.Series(aa).value_counts(normalize=True)
    return round(float(-(p * np.log2(p)).sum()), 3)


rows, pos = [], 0
for i, h in enumerate(human.seq):
    if h == "-":
        continue
    pos += 1
    rows.append({
        "position": pos,
        "human_aa": h,
        "entropy": entropy(aln[:, i]),
    })

cons = pd.DataFrame(rows)
cons.to_csv("data/processed/conservation.csv", index=False)

struct = pd.read_csv("data/processed/variants_structural.csv")
final = struct.merge(cons[["position", "entropy"]], on="position", how="left")
final.to_csv("data/processed/variants_final.csv", index=False)

print(f"Alignment: {len(aln)} sequences, {aln.get_alignment_length()} columns")
print(f"Saved data/processed/conservation.csv ({len(cons)} human positions)")
print(f"Saved data/processed/variants_final.csv ({len(final)} rows)\n")

cols = [
    "allele", "protein_change", "function",
    "in_structure", "min_dist_heme_Fe_A", "rel_sasa", "entropy",
]
print(final[cols].to_string(index=False))