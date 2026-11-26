import pathlib
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

PDB_FILE = "structures/2F9Q.pdb"
CHAIN = "A"
FINAL_CSV = "data/processed/variants_final.csv"
RESULTS_MD = pathlib.Path("docs/results.md")

TARGETS = [
    (
        "R440H", "*31", 440, "no function",
        "Sits in the heme-binding loop just upstream of the heme-coordinating cysteine (C443) on the proximal face of the enzyme; substituting histidine shortens the basic side chain and alters local electrostatic/polar packing around the heme-binding region."
    ),
    (
        "F120I", "*49", 120, "decreased function",
        "Sits on the B-C loop directly lining the hydrophobic active-site pocket over the heme; replacing the planar aromatic ring of phenylalanine with branched isoleucine alters active-site geometry and removes aromatic stacking interactions."
    ),
    (
        "H324P", "*7", 324, "no function",
        "Positioned on the I-helix; introducing a rigid proline removes the backbone amide H-bond donor and induces helical distortion that destabilises local packing."
    ),
    (
        "R296C", "*2", 296, "normal function",
        "Positioned near the N-terminal end of the I-helix, distal to the catalytic heme iron; serves as a tolerated normal-function backbone control shared across multiple star alleles."
    ),
]

# Load 2F9Q Chain A without crystallographic waters
model = PDBParser(QUIET=True).get_structure("cyp2d6", PDB_FILE)[0]
for ch in model:
    for rid in [r.id for r in ch if r.id[0] == "W"]:
        ch.detach_child(rid)
chain = model[CHAIN]

df = pd.read_csv(FINAL_CSV).set_index("protein_change")


def inspect_site(pos):
    res = chain[(" ", pos, " ")]
    t_atoms = [a for a in res.get_atoms() if a.element != "H"]
    neighbors = []
    polar = []

    for other in chain:
        if other == res:
            continue
        resname = other.get_resname().strip()
        if other.id[0] != " " and resname != "HEM":
            continue

        o_atoms = [a for a in other.get_atoms() if a.element != "H"]
        if not o_atoms:
            continue

        min_d = min(
            float(np.linalg.norm(a1.coord - a2.coord))
            for a1 in t_atoms for a2 in o_atoms
        )
        if min_d <= 5.0:
            neighbors.append((min_d, f"{resname}{other.id[1]}"))

        # Candidate polar contacts (N/O heavy atoms within 2.2-3.5 A)
        for a1 in t_atoms:
            if a1.element not in {"N", "O"}:
                continue
            for a2 in o_atoms:
                if a2.element not in {"N", "O"}:
                    continue
                # Skip immediate peptide-bond backbone N/O neighbors (i +/- 1)
                if (
                    other.id[0] == " "
                    and abs(other.id[1] - pos) == 1
                    and a1.get_name() in {"N", "O"}
                    and a2.get_name() in {"N", "O"}
                ):
                    continue
                d = float(np.linalg.norm(a1.coord - a2.coord))
                if 2.2 <= d <= 3.5:
                    polar.append((d, f"{a1.get_name()}-{resname}{other.id[1]}:{a2.get_name()} ({d:.2f} A)"))

    neighbors.sort()
    polar.sort()
    nbr_str = ", ".join(label for _, label in neighbors) if neighbors else "None within 5.0 A"
    pol_str = "; ".join(label for _, label in polar[:8]) if polar else "No N/O contacts <= 3.5 A (primarily hydrophobic/van der Waals packing)"
    return nbr_str, pol_str


lines = [
    "",
    "## Day 9: Local Environment & Illustrative Mutagenesis Notes (2026-11-25)",
    "",
    "Local 5.0 A heavy-atom environments and candidate polar contacts (N/O pairs <= 3.5 A) extracted from PDB `2F9Q` Chain A for four showcase variants spanning no-function, decreased-function, and normal-function classes:",
    "",
]

for pc, allele, pos, fn, context in TARGETS:
    nbr_str, pol_str = inspect_site(pos)
    row = df.loc[pc]
    lines.extend([
        f"### `{pc}` (`{allele}`, {fn})",
        f"- **Quantitative metrics (`2F9Q` Chain A):** Heme Fe distance = `{row['min_dist_heme_Fe_A']} A` | Relative SASA = `{row['rel_sasa']}` | CYP2D-family entropy = `{row['entropy']}`",
        f"- **Residues within 5.0 A:** {nbr_str}",
        f"- **Candidate polar contacts (<= 3.5 A):** {pol_str}",
        f"- **Structural interpretation:** {context}",
        "",
    ])

lines.extend([
    "> **Illustrative Mutagenesis Caveat:** Side-by-side wild-type vs. mutant figures (`figures/closeups/mutagenesis_R440H.png`, `mutagenesis_F120I.png`, `mutagenesis_H324P.png`) were generated using PyMOL's Mutagenesis Wizard (backbone-dependent rotamer selection with minimal steric clashes). These views are strictly illustrative—no energy minimization, backbone relaxation, or folding free-energy ($\Delta\Delta G$) calculation was performed.",
    "",
])

section = "\n".join(lines)
existing = RESULTS_MD.read_text(encoding="utf-8") if RESULTS_MD.exists() else ""

if "## Day 9:" in existing:
    print("Day 9 section already exists in docs/results.md — skipping append to avoid duplicates.")
else:
    RESULTS_MD.parent.mkdir(parents=True, exist_ok=True)
    RESULTS_MD.write_text(existing.rstrip() + "\n" + section, encoding="utf-8")
    print("Appended verified Day 9 Local Environment Notes to docs/results.md!\n")

print(section)