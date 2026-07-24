"""
Plot |Q_h|/(ω1+ω2) on a 2D (ω1, ω2) grid from pre-computed data files.

Edit the CHANNEL section below and run:
    python plot_q2d.py

Data file expected: ../data/q2d_{parity}_{l1}{m1}_{l2}{m2}_{l}.dat
Columns: w1  w2  |Q1|/w_out  |Q2|/w_out  err/w_out

Two panels are produced:
  left  — |Q2|/(ω1+ω2)  (better-converged estimate)
  right — |ΔQ|/(ω1+ω2) = convergence error
"""

import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm

# ---------------------------------------------------------------------------
# Color palette
# ---------------------------------------------------------------------------

def rgb(r, g, b):
    return (r / 255, g / 255, b / 255)

COLORS = {
    "lightpeach":    rgb(255, 229, 180),
    "apricot":       rgb(255, 205, 117),
    "amberglow":     rgb(255, 160,  58),
    "sunsetorange":  rgb(255, 120,  40),
    "coralrose":     rgb(255,  94,  77),
    "burntsienna":   rgb(220,  70,  40),
    "russet":        rgb(170,  50,  30),
    "icemist":       rgb(220, 240, 255),
    "mistlavender":  rgb(180, 200, 255),
    "glacierblue":   rgb(130, 180, 255),
    "frostedteal":   rgb( 58, 183, 183),
    "deepcerulean":  rgb(  0, 123, 167),
    "midnightblue":  rgb(  0,  70, 140),
    "arcticnavy":    rgb( 20,  40,  80),
    "sageleaf":      rgb(141, 182,   0),
    "crimsonred":    rgb(200,  50,  70),
    "tealstone":     rgb(  0, 105,  98),
    "ambergold":     rgb(230, 150,  40),
    "lilacmist":     rgb(150, 110, 180),
}

# ---------------------------------------------------------------------------
# Matplotlib style
# ---------------------------------------------------------------------------

FIG_W_IN = 2 * 9.52 / 2.54   # two panels side by side
FIG_H_IN =     7.14 / 2.54

MAJOR_TICK_PT = 1.1 * 2.835
MINOR_TICK_PT = 0.75 * 2.835
LINE_WIDTH = 0.7

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
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "axes.linewidth":     0.5,
    "axes.xmargin":       0,
    "axes.ymargin":       0,
    "lines.linewidth":    LINE_WIDTH,
    "axes.grid":          False,
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

def data_path(parity, l1, m1, l2, m2, l):
    fname = f"q2d_{parity}_{l1}{m1}_{l2}{m2}_{l}.dat"
    return os.path.join(DATA_DIR, fname)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

def load_q2d(parity, l1, m1, l2, m2, l):
    fpath = data_path(parity, l1, m1, l2, m2, l)
    if not os.path.isfile(fpath):
        raise FileNotFoundError(
            f"Data file not found: {fpath}\n"
            f"Run compute_q2d.py first."
        )

    raw = np.loadtxt(fpath, comments="#")

    # Infer grid size from metadata comment "reshape: data.reshape(N1, N2, -1)"
    N1 = N2 = None
    with open(fpath) as fh:
        for line in fh:
            if "reshape" in line and "reshape(" in line:
                import re
                m = re.search(r"reshape\((\d+),\s*(\d+)", line)
                if m:
                    N1, N2 = int(m.group(1)), int(m.group(2))
                break

    if N1 is None:
        # Fall back: assume square grid
        n = int(round(len(raw) ** 0.5))
        N1 = N2 = n

    grid = raw.reshape(N1, N2, -1)   # (N1, N2, ncols)
    # Columns: 0=w1, 1=w2, 2=|Q1|/w_out, 3=|Q2|/w_out, 4=err/w_out
    w1      = grid[:, 0, 0]          # shape (N1,)
    w2      = grid[0, :, 1]          # shape (N2,)
    abs_q2  = grid[:, :, 3]          # |Q2|/w_out — better converged
    abs_err = grid[:, :, 4]          # err/w_out  — convergence error

    return w1, w2, abs_q2, abs_err

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def format_m(m):
    return f"+{m}" if m > 0 else str(m)

def make_figure(parity, l1, m1, l2, m2, l, log_scale=False):
    w1, w2, abs_q2, abs_err = load_q2d(parity, l1, m1, l2, m2, l)

    fig, axes = plt.subplots(1, 2, figsize=(FIG_W_IN, FIG_H_IN),
                             sharey=True, sharex=False)

    title_base = (
        r"$(" + str(l1) + r"," + format_m(m1) + r")\times("
        + str(l2) + r"," + format_m(m2) + r")\to(" + str(l) + r", " + format_m(m1+m2) + r")$"
        + f"  [{parity.upper()}]"
    )

    cmap_main  = "inferno"
    cmap_err   = "viridis"

    for ax, data, cmap, cbar_label, title_suffix, use_log in [
        (axes[0], abs_q2,  cmap_main, "",       r"$|\mathcal{Q}|/\Omega_d$",      log_scale),
        (axes[1], abs_err, cmap_err,  "",  r"$|\Delta\mathcal{Q}|/\Omega_d$",  True),
    ]:
        pos = data[data > 0]
        norm = LogNorm(vmin=pos.min(), vmax=data.max()) if use_log else None

        pcm = ax.pcolormesh(w2, w1, data,
                            cmap=cmap, norm=norm,
                            shading="auto", rasterized=True)

        cbar = fig.colorbar(pcm, ax=ax, pad=0.02, aspect=20)
        cbar.set_label(cbar_label, fontsize=8)
        cbar.ax.tick_params(labelsize=7)

        ax.set_xlabel(r"$\omega_2 M$")
        ax.set_title(title_suffix, fontsize=8, pad=3)

        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

    axes[0].set_ylabel(r"$\omega_1 M$")
    fig.suptitle(title_base, fontsize=9, y=1.01)
    fig.tight_layout(pad=0.4)
    return fig

# ---------------------------------------------------------------------------
# CHANNEL — edit these values and run the script
# ---------------------------------------------------------------------------

parity = "ooe"
l1 = 2
m1 = 2
l2 = 4
m2 = -2
l  = 2

LOG_SCALE = False   # set True for log-colour scale
SAVE      = True    # write PDF to ../data/

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

fig = make_figure(parity, l1, m1, l2, m2, l, log_scale=LOG_SCALE)

if SAVE:
    out = os.path.join(DATA_DIR, f"plot_q2d_{parity}_{l1}{m1}_{l2}{m2}_{l}.pdf")
    fig.savefig(out)
    print(f"Saved: {out}")

# plt.show()
