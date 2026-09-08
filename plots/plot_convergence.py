"""
Plot the inner/outer integration-radius convergence test, one line per
sector (parity combination).

Two log-log panels: relative error of |Q_h| (vs. the most accurate point in
each sweep) against (rmin-2M)/M and against rmax/M.

Data expected from test/convergence_test.jl, one pair of files per sector:
    ../data/convergence_{tag}_inner.dat   columns: (rmin-2M)/M  rel_err
    ../data/convergence_{tag}_outer.dat   columns: rmax/M       rel_err
where tag = {sector}_{l1}{m1}_{l2}{m2}_{l}.

Usage:
    python3 plot_convergence.py
"""

import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Color palette — exact RGB from the project .tex files (see plot_q_ooo.py)
# ---------------------------------------------------------------------------

def rgb(r, g, b):
    return (r / 255, g / 255, b / 255)

COLORS = {
    "burntsienna":   rgb(220, 70, 40),
    "deepcerulean":  rgb(0, 123, 167),
    "sageleaf":      rgb(141, 182, 0),
    "lilacmist":     rgb(150, 110, 180),
    "ambergold":     rgb(230, 150, 40),
    "midnightblue":  rgb(0, 70, 140),
}

FIG_W_IN = 9.52 / 2.54
FIG_H_IN = 7.14 / 2.54
MAJOR_TICK_PT = 1.1 * 2.835
MINOR_TICK_PT = 0.75 * 2.835
LINE_WIDTH = 0.9
MARKER_SIZE = 4

mpl.rcParams.update({
    "font.family":        "serif",
    "font.size":          9,
    "axes.labelsize":     9,
    "legend.fontsize":    7,
    "xtick.labelsize":    8,
    "ytick.labelsize":    8,
    "text.usetex":        False,
    "mathtext.fontset":   "cm",
    "xtick.direction":    "in",
    "ytick.direction":    "in",
    "xtick.top":          True,
    "ytick.right":        True,
    "xtick.major.size":   MAJOR_TICK_PT,
    "xtick.minor.size":   MINOR_TICK_PT,
    "ytick.major.size":   MAJOR_TICK_PT,
    "ytick.minor.size":   MINOR_TICK_PT,
    "xtick.major.width":  0.4,
    "xtick.minor.width":  0.4,
    "ytick.major.width":  0.4,
    "ytick.minor.width":  0.4,
    "axes.linewidth":     0.5,
    "axes.grid":          False,
    "figure.dpi":         150,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
    "savefig.pad_inches": 0.01,
})

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(SCRIPT_DIR, "..", "data")

# sector -> (l1,m1,l2,m2,l), must match test/convergence_test.jl's CHANNELS
CHANNELS = {
    "eee": (2, 2, 2, -2, 2),
    "eeo": (2, 2, 3, -2, 2),
    "eoe": (2, 2, 2, -2, 3),
    "eoo": (2, 2, 2, -2, 2),
    "ooe": (2, 2, 2, -2, 2),
    "ooo": (2, 2, 3, -2, 2),
}
SECTOR_COLOR = {
    "eee": "burntsienna",
    "eeo": "deepcerulean",
    "eoe": "sageleaf",
    "eoo": "lilacmist",
    "ooe": "ambergold",
    "ooo": "midnightblue",
}
MARKERS = {
    "eee": "o", "eeo": "s", "eoe": "^", "eoo": "v", "ooe": "D", "ooo": "P",
}


def main():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(2 * FIG_W_IN, FIG_H_IN), sharey=True)

    for sector in CHANNELS:
        l1, m1, l2, m2, l = CHANNELS[sector]
        tag = f"{sector}_{l1}{m1}_{l2}{m2}_{l}"
        inner_path = os.path.join(DATA_DIR, f"convergence_{tag}_inner.dat")
        outer_path = os.path.join(DATA_DIR, f"convergence_{tag}_outer.dat")
        if not (os.path.isfile(inner_path) and os.path.isfile(outer_path)):
            print(f"skipping {sector}: data not found (run test/convergence_test.jl first)")
            continue

        inner = np.loadtxt(inner_path, comments="#")
        outer = np.loadtxt(outer_path, comments="#")
        c = COLORS[SECTOR_COLOR[sector]]
        mk = MARKERS[sector]

        ax1.plot(inner[:, 0], inner[:, 1], marker=mk, color=c,
                  linewidth=LINE_WIDTH, markersize=MARKER_SIZE, label=sector.upper())
        ax2.plot(outer[:, 0], outer[:, 1], marker=mk, color=c,
                  linewidth=LINE_WIDTH, markersize=MARKER_SIZE, label=sector.upper())

    for ax, xlabel in ((ax1, r"$(r_{\rm min}-2M)/M$"), (ax2, r"$r_{\rm max}/M$")):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(xlabel)

    ax1.set_ylabel(r"relative error in $|\mathcal{Q}_h|$")
    ax1.legend(loc="best", ncol=2)

    fig.tight_layout(pad=0.3)

    out = os.path.join(DATA_DIR, "convergence_all_sectors.pdf")
    fig.savefig(out)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
