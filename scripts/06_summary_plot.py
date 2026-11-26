import pathlib
import matplotlib
matplotlib.use("Agg")  # Headless rendering without popping up GUI windows
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

FINAL_CSV = "data/processed/variants_final.csv"
RAW_CSV = "data/raw/cyp2d6_variants_curated.csv"
FIG_DIR = pathlib.Path("figures")

# Consistent color palette matching PyMOL renders (Days 8 & 10)
COLOR_MAP = {
    "normal function": "#2b6cb0",     # marine blue
    "decreased function": "#ed8936",  # orange
    "no function": "#e53e3e",         # red
    "unknown function": "#718096",    # gray50
}

CLASS_ORDER = [
    "no function",
    "decreased function",
    "normal function",
    "unknown function",
]

LABEL_PRETTY = {
    "normal function": "Normal function",
    "decreased function": "Decreased function",
    "no function": "No function",
    "unknown function": "Unknown / uncertain",
}


def place_labels(ax, xs, ys, texts):
    """Add point labels with simple deterministic vertical nudging to reduce overlap."""
    x_span = max(xs) - min(xs) if len(xs) > 1 else 1.0
    y_span = max(ys) - min(ys) if len(ys) > 1 else 1.0
    x_tol = 0.13 * x_span
    y_tol = 0.065 * y_span

    placed = []
    # Sort by x then y for predictable placement
    order = np.argsort(xs)
    for idx in order:
        x, y, txt = float(xs[idx]), float(ys[idx]), texts[idx]
        dy = 0.025 * y_span
        dx = 0.015 * x_span

        # Nudge vertically if too close to an already-placed label
        for _ in range(8):
            clash = any(
                abs((x + dx) - px) < x_tol and abs((y + dy) - py) < y_tol
                for px, py in placed
            )
            if not clash:
                break
            dy = -dy if dy > 0 else abs(dy) + 0.045 * y_span

        placed.append((x + dx, y + dy))
        ax.annotate(
            txt,
            xy=(x, y),
            xytext=(x + dx, y + dy),
            fontsize=8.5,
            fontweight="normal",
            color="#1a202c",
            arrowprops=dict(arrowstyle="-", color="#a0aec0", lw=0.6, shrinkA=0, shrinkB=4),
            ha="left",
            va="center",
        )


def plot_structural_metrics(mapped):
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6), dpi=300)
    fig.patch.set_facecolor("white")

    panels = [
        (
            axes[0],
            "rel_sasa",
            "Relative Solvent Accessibility (Shrake-Rupley)",
            "A. Heme Distance vs. Solvent Accessibility",
            0.20,
            "Buried threshold (rel_sasa = 0.20)",
        ),
        (
            axes[1],
            "entropy",
            "CYP2D-Family Shannon Entropy (bits; 0 = 100% conserved)",
            "B. Heme Distance vs. Evolutionary Conservation",
            0.50,
            "High conservation (H <= 0.50 bits)",
        ),
    ]

    for ax, y_col, y_label, title, y_ref, ref_text in panels:
        sub = mapped.dropna(subset=["min_dist_heme_Fe_A", y_col]).copy()

        # Subtle reference lines for structural interpretation
        ax.axvline(10.0, color="#cbd5e0", linestyle="--", linewidth=1.0, zorder=1)
        ax.axhline(y_ref, color="#cbd5e0", linestyle=":", linewidth=1.0, zorder=1)

        for fn in CLASS_ORDER:
            grp = sub[sub["function"] == fn]
            if grp.empty:
                continue
            ax.scatter(
                grp["min_dist_heme_Fe_A"],
                grp[y_col],
                c=COLOR_MAP[fn],
                label=f"{LABEL_PRETTY[fn]} (n={len(grp)})",
                s=95,
                edgecolors="white",
                linewidths=1.1,
                zorder=3,
            )

        place_labels(
            ax,
            sub["min_dist_heme_Fe_A"].values,
            sub[y_col].values,
            [f"{pc} ({al})" for pc, al in zip(sub["protein_change"], sub["allele"])],
        )

        ax.set_title(title, fontsize=11.5, fontweight="bold", pad=10, color="#1a202c")
        ax.set_xlabel("Min. Heavy-Atom Distance to Heme Fe (A)", fontsize=10, color="#2d3748")
        ax.set_ylabel(y_label, fontsize=10, color="#2d3748")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, linestyle="-", alpha=0.18, zorder=0)

        # Pad x/y limits slightly so annotations never clip on the right edge
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        ax.set_xlim(max(0, xmin - 1.5), xmax + 5.0)
        ax.set_ylim(max(-0.08, ymin - 0.08 * (ymax - ymin)), ymax + 0.14 * (ymax - ymin))

        ax.text(
            10.3, ax.get_ylim()[1] * 0.96, "10 A heme zone",
            fontsize=7.8, color="#718096", va="top", ha="left",
        )
        ax.text(
            ax.get_xlim()[1] * 0.98, y_ref + (-0.018 if y_col == "rel_sasa" else 0.02) * (ymax - ymin), ref_text,
            fontsize=7.8, color="#718096", va=("top" if y_col == "rel_sasa" else "bottom"), ha="right",
        )

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=4,
        frameon=False,
        fontsize=9.5,
        bbox_to_anchor=(0.5, -0.02),
    )

    fig.suptitle(
        f"CYP2D6 Mapped Missense Variants (PDB: 2F9Q Chain A, n = {len(mapped)})",
        fontsize=13,
        fontweight="bold",
        y=0.99,
        color="#1a202c",
    )
    plt.tight_layout(rect=[0, 0.05, 1, 0.94])

    out_path = FIG_DIR / "summary_structural_metrics.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_sas_frequencies(raw_df):
    df = raw_df.copy()
    df["freq_num"] = pd.to_numeric(df["freq_sas"], errors="coerce")
    df = df.dropna(subset=["freq_num"]).copy()

    # Keep only detected alleles (> 0%) and deduplicate multi-variant star alleles (e.g. *2)
    df = df[df["freq_num"] > 0].drop_duplicates(subset=["allele"]).copy()
    df["is_missense"] = df["variant_type"].str.strip().str.lower().eq("missense")
    df = df.sort_values("freq_num", ascending=True)

    fig, ax = plt.subplots(figsize=(9.5, 5.2), dpi=300)
    fig.patch.set_facecolor("white")

    y_labels = [
        f"{row['allele']}  ({row['protein_change']})"
        for _, row in df.iterrows()
    ]
    colors = [COLOR_MAP.get(fn, "#718096") for fn in df["function"]]
    hatches = ["" if m else "///" for m in df["is_missense"]]

    bars = ax.barh(
        y_labels,
        df["freq_num"] * 100.0,
        color=colors,
        edgecolor="#1a202c",
        linewidth=0.7,
        height=0.65,
    )
    for bar, h in zip(bars, hatches):
        if h:
            bar.set_hatch(h)

    for bar, freq in zip(bars, df["freq_num"]):
        pct = freq * 100.0
        ax.text(
            bar.get_width() + 0.35,
            bar.get_y() + bar.get_height() / 2,
            f"{pct:.2f}%",
            va="center",
            ha="left",
            fontsize=8.8,
            color="#1a202c",
        )

    ax.set_xlabel("South Asian / Indian Allele Frequency (%)", fontsize=10, color="#2d3748")
    ax.set_title(
        "South Asian Allele Frequencies: 3D-Mappable Missense vs. Non-Mappable Alleles",
        fontsize=11.5,
        fontweight="bold",
        pad=12,
        color="#1a202c",
    )
    ax.set_xlim(0, df["freq_num"].max() * 100.0 + 3.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", linestyle="-", alpha=0.2)

    legend_items = [
        Patch(facecolor=COLOR_MAP["normal function"], edgecolor="#1a202c", label="Normal function"),
        Patch(facecolor=COLOR_MAP["decreased function"], edgecolor="#1a202c", label="Decreased function"),
        Patch(facecolor=COLOR_MAP["no function"], edgecolor="#1a202c", label="No function"),
        Patch(facecolor=COLOR_MAP["unknown function"], edgecolor="#1a202c", label="Unknown / uncertain"),
        Patch(facecolor="white", edgecolor="#1a202c", hatch="///", label="Non-mappable / in-frame (splice, del)"),
    ]
    ax.legend(handles=legend_items, loc="lower right", frameon=True, facecolor="white", fontsize=8.5)

    fig.text(
        0.01, -0.02,
        "Sources: SG10K_Health Indian cohort (N=328; Maulana et al., 2024) and supplementary gnomAD SAS (*22). "
        "Alleles with 0% detected or unavailable SAS data omitted.",
        fontsize=7.8, color="#4a5568",
    )

    plt.tight_layout()
    out_path = FIG_DIR / "summary_sas_frequencies.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")


def print_summary_stats(final_df, mapped):
    print("\n=== SUMMARY STATISTICS BY CPIC FUNCTION CLASS (Mapped in 2F9Q Chain A) ===")
    for fn in CLASS_ORDER:
        grp = mapped[mapped["function"] == fn]
        if grp.empty:
            continue
        print(
            f"{fn:20s} | n={len(grp):2d} | "
            f"Heme dist median={grp['min_dist_heme_Fe_A'].median():5.2f} A "
            f"[{grp['min_dist_heme_Fe_A'].min():.2f}-{grp['min_dist_heme_Fe_A'].max():.2f}] | "
            f"rel_sasa median={grp['rel_sasa'].median():.2f} "
            f"[{grp['rel_sasa'].min():.2f}-{grp['rel_sasa'].max():.2f}] | "
            f"entropy median={grp['entropy'].median():.3f} "
            f"[{grp['entropy'].min():.3f}-{grp['entropy'].max():.3f}]"
        )

    unres = final_df[final_df["in_structure"].astype(str).str.lower() != "true"]
    if not unres.empty:
        print("\nUnresolved in 2F9Q Chain A (excluded from 3D structural scatter):")
        print(unres[["allele", "protein_change", "function", "entropy"]].to_string(index=False))


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    final_df = pd.read_csv(FINAL_CSV)
    raw_df = pd.read_csv(RAW_CSV)

    mapped = final_df[final_df["in_structure"].astype(str).str.lower() == "true"].copy()
    for col in ["min_dist_heme_Fe_A", "rel_sasa", "entropy"]:
        mapped[col] = pd.to_numeric(mapped[col], errors="coerce")
    mapped["entropy"] = mapped["entropy"].abs()

    plot_structural_metrics(mapped)
    plot_sas_frequencies(raw_df)
    print_summary_stats(final_df, mapped)


if __name__ == "__main__":
    main()