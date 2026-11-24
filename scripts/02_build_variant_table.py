import pathlib
import re
import pandas as pd

raw = pd.read_csv("data/raw/cyp2d6_variants_curated.csv")

# Matches single-letter missense changes like R296C or p.R296C
pattern = re.compile(r"^(?:p\.)?([A-Z])(\d+)([A-Z])$")


def parse(pc):
    m = pattern.match(str(pc).strip())
    return m.groups() if m else (None, None, None)


raw[["wt", "position", "mut"]] = raw["protein_change"].apply(
    lambda s: pd.Series(parse(s))
)

mappable = raw[raw["variant_type"].eq("missense") & raw["position"].notna()].copy()
mappable["position"] = mappable["position"].astype(int)
unmapped = raw.drop(mappable.index)

# Ensure data/processed/ exists on clean clones
out_dir = pathlib.Path("data/processed")
out_dir.mkdir(parents=True, exist_ok=True)

mappable.to_csv(out_dir / "variants_mappable.csv", index=False)
unmapped.to_csv(out_dir / "variants_not_mappable.csv", index=False)

print(f"{len(mappable)} mappable; {len(unmapped)} not mappable")

# Cross-check every mappable wild-type residue against UniProt P10635
fasta_path = pathlib.Path("data/raw/P10635.fasta")
if fasta_path.exists():
    seq = "".join(
        line.strip()
        for line in fasta_path.read_text().splitlines()
        if not line.startswith(">")
    )
    mismatches = []
    for _, row in mappable.iterrows():
        pos = row["position"]
        uniprot_aa = seq[pos - 1] if 1 <= pos <= len(seq) else "?"
        if uniprot_aa != row["wt"]:
            mismatches.append((row["allele"], row["protein_change"], uniprot_aa))
    if mismatches:
        print("WARNING - UniProt P10635 mismatches:", mismatches)
    else:
        print(
            f"UniProt P10635 check passed: all {len(mappable)} wild-type residues match "
            f"(reference length {len(seq)} aa)."
        )