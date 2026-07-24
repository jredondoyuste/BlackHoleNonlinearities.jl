"""
Plot Q vs ω for an OOO channel from pre-computed data files.

Edit the CHANNEL section below to choose which channel to plot, then run:
    python plot_q_ooo.py

Data file expected: ../data/q_ooo_{l1}{m1}_{l2}{m2}_{l}.dat
Columns: w  Re(Q_h)  Im(Q_h)  |Q_h|/w
"""

import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ---------------------------------------------------------------------------
# Color palette — exact RGB from the project .tex files
# ---------------------------------------------------------------------------

def rgb(r, g, b):
    return (r / 255, g / 255, b / 255)

COLORS = {
    # warm
    "lightpeach":    rgb(255, 229, 180),
    "apricot":       rgb(255, 205, 117),
    "amberglow":     rgb(255, 160,  58),
    "sunsetorange":  rgb(255, 120,  40),
    "coralrose":     rgb(255,  94,  77),
    "burntsienna":   rgb(220,  70,  40),
    "russet":        rgb(170,  50,  30),
    # cool
    "icemist":       rgb(220, 240, 255),
    "mistlavender":  rgb(180, 200, 255),
    "glacierblue":   rgb(130, 180, 255),
    "frostedteal":   rgb( 58, 183, 183),
    "deepcerulean":  rgb(  0, 123, 167),
    "midnightblue":  rgb(  0,  70, 140),
    "arcticnavy":    rgb( 20,  40,  80),
    # accent
    "sageleaf":      rgb(141, 182,   0),
    "crimsonred":    rgb(200,  50,  70),
    "tealstone":     rgb(  0, 105,  98),
    "ambergold":     rgb(230, 150,  40),
    "lilacmist":     rgb(150, 110, 180),
}

# ---------------------------------------------------------------------------
# Matplotlib style matching pgfplots defaults in the .tex files
# ---------------------------------------------------------------------------

# Figure size: 9.52cm × 7.14cm  (same as \axisdefaultwidth / \axisdefaultheight)
FIG_W_IN = 9.52 / 2.54
FIG_H_IN = 7.14 / 2.54

# Tick lengths in points (1 mm = 2.835 pt)
MAJOR_TICK_PT = 1.1 * 2.835   # 3.12 pt
MINOR_TICK_PT = 0.75 * 2.835  # 2.13 pt

# semithick in pgfplots ≈ 0.6 pt; use 0.7 for a comfortable screen rendering
LINE_WIDTH = 0.7

mpl.rcParams.update({
    # Font
    "font.family":        "serif",
    "font.size":          9,
    "axes.labelsize":     9,
    "legend.fontsize":    7,
    "xtick.labelsize":    8,
    "ytick.labelsize":    8,
    # Math (use built-in mathtext to avoid needing a LaTeX install)
    "text.usetex":        False,
    "mathtext.fontset":   "cm",
    # Ticks
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
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    # Axes
    "axes.linewidth":     0.5,
    "axes.xmargin":       0,    # enlarge x limits=false
    "axes.ymargin":       0.04,
    # Legend
    "legend.framealpha":  1.0,
    "legend.edgecolor":   "black",
    "legend.fancybox":    False,
    "legend.borderpad":   0.4,
    "legend.labelspacing": 0.3,
    # Lines
    "lines.linewidth":    LINE_WIDTH,
    # Grid
    "axes.grid":          False,
    # Figure
    "figure.dpi":         150,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
    "savefig.pad_inches": 0.01,
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(SCRIPT_DIR, "..", "data")

def data_path(l1, m1, l2, m2, l):
    fname = f"q_ooo_{l1}{m1}_{l2}{m2}_{l}.dat"
    return os.path.join(DATA_DIR, fname)

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

DASHES_SOLID  = (None, None)
DASHES_DENSE  = (4, 2)          # densely dashed
DASHES_DOT    = (1, 2)          # densely dotted

def format_m(m):
    """Format magnetic quantum number with sign for labels."""
    return f"+{m}" if m > 0 else str(m)

def channel_label(l1, m1, l2, m2, l):
    return (
        r"$|\mathcal{Q}_{(" + str(l1) + "," + format_m(m1) + r")\times("
        + str(l2) + "," + format_m(m2) + r")\to(" + str(l) + r",0)}|\,/\,\Omega_d$"
    )

def plot_channel(l1, m1, l2, m2, l, ax=None, color="crimsonred",
                 dashes=DASHES_SOLID, label=None):
    """
    Load data for channel (l1,m1)×(l2,m2)→(l,0) and plot |Q_h|/ω vs ω.
    Returns the axes object.
    """
    fpath = data_path(l1, m1, l2, m2, l)
    if not os.path.isfile(fpath):
        raise FileNotFoundError(
            f"Data file not found: {fpath}\n"
            f"Run compute_q_ooo.py first for channel "
            f"(l1={l1},m1={m1},l2={l2},m2={m2},l={l})."
        )

    data = np.loadtxt(fpath, comments="#")
    w      = data[:, 0]
    abs_qh = data[:, 3]   # |Q_h|/w  (column 4)

    if ax is None:
        fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))

    c = COLORS.get(color, color)
    lw = LINE_WIDTH
    line_label = label if label is not None else channel_label(l1, m1, l2, m2, l)

    if dashes == DASHES_SOLID:
        ax.plot(w, abs_qh, color=c, linewidth=lw, label=line_label)
    else:
        ax.plot(w, abs_qh, color=c, linewidth=lw, label=line_label,
                dashes=dashes)

    return ax


def make_figure(l1, m1, l2, m2, l):
    """
    Produce a single-channel Q vs ω figure matching the .tex style.
    """
    fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))

    plot_channel(l1, m1, l2, m2, l, ax=ax, color="crimsonred")

    # Axes labels  (match TeX: x = Ω_d M, y = |Q|/Ω_d)
    ax.set_xlabel(r"$\Omega_d M$")
    ax.set_ylabel(channel_label(l1, m1, l2, m2, l))

    ax.set_xlim(left=0)
    # ax.set_xlim(right=0.5)
    ax.set_ylim(bottom=0)
    # ax.set_ylim(top=0.0012)

    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

    # Legend top-left (matches \at={(0.02,0.98)}, anchor=north west)
    ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.98))

    fig.tight_layout(pad=0.3)
    return fig


# ---------------------------------------------------------------------------
# CHANNEL — edit these values and run the script
# ---------------------------------------------------------------------------

l1 = 5
m1 = 5
l2 = 5
m2 = -5
l  = 3

SAVE = True   # set True to also write a PDF to ../data/

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

if (l1 + l2 + l) % 2 == 0:
    print(f"WARNING: l1+l2+l = {l1+l2+l} is even — this channel is zero by the OOO parity rule.")
if not (abs(l1 - l2) <= l <= l1 + l2):
    print(f"WARNING: l={l} violates triangle rule |{l1}-{l2}|={abs(l1-l2)} ≤ l ≤ {l1+l2}.")

fig = make_figure(l1, m1, l2, m2, l)

if SAVE:
    out = os.path.join(DATA_DIR, f"plot_q_ooo_{l1}{m1}_{l2}{m2}_{l}.pdf")
    fig.savefig(out)
    print(f"Saved: {out}")

# plt.show()
