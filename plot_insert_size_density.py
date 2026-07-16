#!/usr/bin/env python3
"""Generate an insert size density plot PDF from ATP2 Guenther et al. PNAS 2015."""

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from collections import Counter

INPUT_FILE = "ATP2_guenther_count.txt"
OUTPUT_PDF = "ATP2_guenther_insert_size_density.pdf"


def load_insert_sizes(path):
    sizes = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                sizes.append(int(line))
    return np.array(sizes)


def main():
    sizes = load_insert_sizes(INPUT_FILE)

    counts = Counter(sizes)
    xs = np.array(sorted(counts.keys()))
    ys = np.array([counts[x] for x in xs], dtype=float)
    total = ys.sum()
    density = ys / total  # normalize to density (sums to ~1 over integer bins)

    fig, ax = plt.subplots(figsize=(8, 4.5))

    # Bar density
    ax.bar(xs, density, width=1.0, color="#4C72B0", alpha=0.55,
           linewidth=0, label="Observed frequency")

    ax.set_xlabel("Insert size (bp)", fontsize=13)
    ax.set_ylabel("Density", fontsize=13)
    ax.set_title(
        "ATP2 cfDNA insert size distribution\nGuenther et al., PNAS 2015",
        fontsize=13, fontweight="bold", pad=10
    )

    ax.set_xlim(xs.min() - 5, xs.max() + 5)
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
        lambda v, _: f"{v:.4f}" if v < 0.001 else f"{v:.3f}"
    ))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=11)

    n_reads = len(sizes)
    ax.text(0.97, 0.95, f"n = {n_reads:,} reads",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=10, color="#444444")

    ax.legend(frameon=False, fontsize=10, loc="upper right",
              bbox_to_anchor=(0.97, 0.88))

    fig.tight_layout()
    fig.savefig(OUTPUT_PDF, dpi=300, bbox_inches="tight")
    print(f"Saved: {OUTPUT_PDF}")


if __name__ == "__main__":
    main()
