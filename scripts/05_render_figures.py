import pathlib
import pandas as pd
from pymol import cmd, util

PDB_FILE = "structures/2F9Q.pdb"
VARIANTS_CSV = "data/processed/variants_final.csv"
FIG_DIR = pathlib.Path("figures")
CLOSEUP_DIR = FIG_DIR / "closeups"
PYMOL_DIR = pathlib.Path("pymol")

COLOR_MAP = {
    "normal function": "marine",
    "decreased function": "orange",
    "no function": "red",
    "unknown function": "gray50",
}


def setup_scaffold():
    cmd.reinitialize()
    cmd.load(PDB_FILE, "cyp2d6")
    cmd.remove("solvent")
    cmd.remove("not chain A")
    cmd.hide("everything")
    cmd.show("cartoon", "polymer")
    cmd.bg_color("white")
    cmd.color("gray80", "polymer")

    cmd.select("heme", "resn HEM")
    cmd.show("sticks", "heme")
    cmd.color("yellow", "heme")

    # Publication styling + Day 8 label-clipping fix
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.set("sphere_scale", 0.65)
    cmd.set("label_color", "black")
    cmd.set("label_outline_color", "white")
    cmd.set("label_size", 16)
    cmd.set("label_position", [0, 0, 4])


def render_all():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    CLOSEUP_DIR.mkdir(parents=True, exist_ok=True)
    PYMOL_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(VARIANTS_CSV)
    mapped = df[df["in_structure"].astype(str).str.lower() == "true"].copy()
    mapped["position"] = mapped["position"].astype(float).astype(int)

    setup_scaffold()

    # 1. Map all resolved missense variants onto the overview scaffold
    for _, row in mapped.iterrows():
        pc = str(row["protein_change"]).strip()
        pos = int(row["position"])
        fn = str(row["function"]).strip().lower()
        color = COLOR_MAP.get(fn, "gray50")
        sel = f"v_{pc}"

        cmd.select(sel, f"polymer and chain A and resi {pos}")
        cmd.show("spheres", f"{sel} and name CA")
        cmd.color(color, sel)
        cmd.label(f"{sel} and name CA", f'"{pc}"')

    cmd.orient("cyp2d6")

    # Save final session and front/back overview renders
    pse_path = PYMOL_DIR / "cyp2d6_atlas_final.pse"
    cmd.save(str(pse_path))
    print(f"Saved session: {pse_path}")

    front_path = FIG_DIR / "overview_front.png"
    print(f"Rendering {front_path} ...")
    cmd.png(str(front_path), width=2000, height=1500, dpi=300, ray=1)

    cmd.turn("y", 180)
    back_path = FIG_DIR / "overview_back.png"
    print(f"Rendering {back_path} ...")
    cmd.png(str(back_path), width=2000, height=1500, dpi=300, ray=1)
    cmd.turn("y", 180)

    # 2. Switch cleanly to stick close-ups (hide overview spheres & labels)
    cmd.hide("spheres", "all")
    cmd.label("all", "")
    cmd.color("gray80", "polymer")
    cmd.set("label_size", 14)
    cmd.set("label_position", [0, 0, 2.5])
    cmd.set("dash_gap", 0.3)

    total = len(mapped)
    for idx, (_, row) in enumerate(mapped.iterrows(), start=1):
        pc = str(row["protein_change"]).strip()
        allele = str(row["allele"]).strip()
        pos = int(row["position"])
        fn = str(row["function"]).strip().lower()
        color = COLOR_MAP.get(fn, "gray50")
        sel = f"v_{pc}"

        cmd.select("nbrs", f"byres ({sel} around 5) and (polymer or heme)")
        cmd.show("sticks", "nbrs")
        cmd.show("sticks", sel)
        cmd.color(color, sel)
        util.cnc("nbrs and not heme")
        util.cnc(sel)
        cmd.distance("hb", sel, f"nbrs and not {sel}", mode=2)
        cmd.label(f"{sel} and name CA", f'"{pc} ({allele})"')

        cmd.zoom(sel, 8)
        out_png = CLOSEUP_DIR / f"closeup_{pc}.png"
        print(f"[{idx}/{total}] Rendering {out_png} ...")
        cmd.png(str(out_png), width=1400, height=1100, dpi=300, ray=1)

        # Clean reset before next variant
        cmd.delete("hb")
        cmd.hide("sticks", f"(nbrs or {sel}) and not heme")
        cmd.label(f"{sel} and name CA", "")
        cmd.color("gray80", "polymer")

    print(f"\nDone! Rendered 2 overview PNGs + {total} close-up PNGs.")


# Supports both `pymol -cq scripts/05_render_figures.py` and `python scripts/05_render_figures.py`
if __name__ == "__main__":
    import pymol
    pymol.finish_launching(["pymol", "-cq"])

render_all()