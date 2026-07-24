"""
Low-frequency fit of |Q_h| for the channels stored in ../data_PM/.

For each data file the magnitude of the (better-converged) Q-factor is fitted
to two low-frequency power laws

    quadratic:  |Q(w)| = A w + B w^2          <=>  |Q|/w = A + B w
    cubic:      |Q(w)| = A w + B w^2 + C w^3  <=>  |Q|/w = A + B w + C w^2

(weighted least squares using the per-point error column).  A single multi-panel
figure shows, for every channel, the data with error bars and both fitted curves
overlaid.

Data files: ../data_PM/*.dat
Header (2 lines):
    # w  |Q1|/w  |Q2|/w  err/w
    # parity=ooe l1=2 l2=2 l=4 m1=2 m2=2 source=bruno lambda1=3.0 lambda2=3.5
Columns: w  |Q1|/w  |Q2|/w  err/w
  |Q2| is the better-converged estimate (central value);
  err = ||Q1| - |Q2||.

Run:
    python plot_q_lowfreq_fit.py
"""

import os
import glob
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
    "amberglow":     rgb(255, 160,  58),
    "sunsetorange":  rgb(255, 120,  40),
    "coralrose":     rgb(255,  94,  77),
    "burntsienna":   rgb(220,  70,  40),
    "crimsonred":    rgb(200,  50,  70),
    "deepcerulean":  rgb(  0, 123, 167),
    "midnightblue":  rgb(  0,  70, 140),
    "frostedteal":   rgb( 58, 183, 183),
    "tealstone":     rgb(  0, 105,  98),
    "sageleaf":      rgb(141, 182,   0),
    "lilacmist":     rgb(150, 110, 180),
    "arcticnavy":    rgb( 20,  40,  80),
}

DATA_COLOR  = COLORS["deepcerulean"]
FIT2_COLOR  = COLORS["crimsonred"]    # quadratic |Q| = A w + B w^2
FIT3_COLOR  = COLORS["sageleaf"]      # cubic     |Q| = A w + B w^2 + C w^3

# ---------------------------------------------------------------------------
# Matplotlib style matching the project's pgfplots look
# ---------------------------------------------------------------------------

mpl.rcParams.update({
    "font.family":         "serif",
    "font.size":           9,
    "axes.labelsize":      9,
    "axes.titlesize":      8,
    "legend.fontsize":     6.5,
    "xtick.labelsize":     7,
    "ytick.labelsize":     7,
    "text.usetex":         False,
    "mathtext.fontset":    "cm",
    "xtick.direction":     "in",
    "ytick.direction":     "in",
    "xtick.top":           True,
    "ytick.right":         True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.width":   0.4,
    "ytick.major.width":   0.4,
    "xtick.minor.width":   0.4,
    "ytick.minor.width":   0.4,
    "axes.linewidth":      0.5,
    "axes.xmargin":        0,
    "axes.ymargin":        0.05,
    "legend.framealpha":   1.0,
    "legend.edgecolor":    "black",
    "legend.fancybox":     False,
    "figure.dpi":          150,
    "savefig.dpi":         300,
    "savefig.bbox":        "tight",
    "savefig.pad_inches":  0.02,
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(SCRIPT_DIR, "..", "data_PM")
OUT_PATH   = os.path.join(SCRIPT_DIR, "q_lowfreq_fit.pdf")

# Restrict the fit to the low-frequency window (None = use all points).
FIT_WMAX = None

# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_header(fpath):
    """Read the 2nd comment line and return its key=value pairs as a dict."""
    meta = {}
    with open(fpath) as f:
        for line in f:
            if not line.startswith("#"):
                break
            for tok in line[1:].split():
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    meta[k] = v
    return meta


def format_m(m):
    return f"+{m}" if m > 0 else str(m)


def channel_label(meta):
    l1, m1 = meta["l1"], meta["m1"]
    l2, m2 = meta["l2"], meta["m2"]
    l = meta["l"]
    m = int(m1) + int(m2)
    return (
        r"$(" + l1 + "," + format_m(int(m1)) + r")\times("
        + l2 + "," + format_m(int(m2)) + r")\to("
        + l + "," + format_m(m) + r")$"
    )

# ---------------------------------------------------------------------------
# Fits on y = |Q|/w, polynomial in w of given degree:
#   degree 1 (quadratic |Q|):  y = A + B w
#   degree 2 (cubic     |Q|):  y = A + B w + C w^2
# Linear in the coefficients -> weighted linear least squares.
# ---------------------------------------------------------------------------

def fit_lowfreq(w, q_over_w, err_over_w, degree, wmax=None):
    """
    Weighted least-squares fit of y = |Q|/w to a polynomial of `degree` in w
    (so |Q| is a polynomial of degree+1 with no constant term).
    Returns (coeffs, sigmas), each of length degree+1, ordered [A, B, (C)].
    """
    mask = np.ones_like(w, dtype=bool)
    if wmax is not None:
        mask = w <= wmax
    ww   = w[mask]
    y    = q_over_w[mask]
    sig  = err_over_w[mask]
    sig  = np.where(sig > 0, sig, np.nanmedian(sig[sig > 0]))

    # design matrix columns: w^0, w^1, ..., w^degree
    Xdes = np.vstack([ww**k for k in range(degree + 1)]).T
    W    = 1.0 / sig**2
    XtW  = Xdes.T * W
    cov  = np.linalg.inv(XtW @ Xdes)
    beta = cov @ (XtW @ y)
    sig_beta = np.sqrt(np.diag(cov))
    return beta, sig_beta


def poly_q_over_w(coeffs, w):
    """Evaluate y = sum_k coeffs[k] w^k."""
    return sum(c * w**k for k, c in enumerate(coeffs))

# ---------------------------------------------------------------------------
# Build the figure
# ---------------------------------------------------------------------------

def main():
    files = sorted(glob.glob(os.path.join(DATA_DIR, "*.dat")))
    if not files:
        raise FileNotFoundError(f"No .dat files found in {DATA_DIR}")

    n = len(files)
    ncols = 3 if n > 4 else 2
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(3.1 * ncols, 2.5 * nrows),
        squeeze=False,
    )
    axes_flat = axes.flatten()

    for ax, fpath in zip(axes_flat, files):
        meta = parse_header(fpath)
        data = np.loadtxt(fpath, comments="#")
        w   = data[:, 0]
        y   = data[:, 2]   # |Q2|/w  (central value)
        err = data[:, 3]   # err/w

        c2, s2 = fit_lowfreq(w, y, err, degree=1, wmax=FIT_WMAX)  # [A, B]
        c3, s3 = fit_lowfreq(w, y, err, degree=2, wmax=FIT_WMAX)  # [A, B, C]

        # data with error bars
        ax.errorbar(
            w, y, yerr=err,
            fmt="o", ms=2.0, mfc=DATA_COLOR, mec=DATA_COLOR,
            ecolor=DATA_COLOR, elinewidth=0.5, capsize=0,
            alpha=0.6, zorder=2, label="data",
        )

        # fit overlays:  |Q|/w = polynomial in w
        wf = np.linspace(w.min(), w.max(), 400)
        ax.plot(wf, poly_q_over_w(c2, wf), color=FIT2_COLOR, lw=1.2, zorder=3,
                label=r"$A\,w + B\,w^2$")
        ax.plot(wf, poly_q_over_w(c3, wf), color=FIT3_COLOR, lw=1.2,
                dashes=(4, 2), zorder=4, label=r"$A\,w + B\,w^2 + C\,w^3$")

        ax.set_title(channel_label(meta), pad=3)
        # ax.set_xlim(left=0)
        ax.set_ylim(bottom=-1e-3)
        ax.set_xscale('log')
        ax.set_xlim(1e-3,1e-1)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

        # annotate fitted coefficients (cubic fit)
        txt = (f"$A = {c3[0]:.3e}$\n$B = {c3[1]:.3e}$\n$C = {c3[2]:.3e}$")
        ax.text(0.04, 0.96, txt, transform=ax.transAxes,
                va="top", ha="left", fontsize=5.5,
                bbox=dict(boxstyle="round,pad=0.25", fc="white",
                          ec="0.7", lw=0.4, alpha=0.9))

        print(f"{os.path.basename(fpath):28s}")
        print(f"    quad:  A = {c2[0]:.6e} ± {s2[0]:.2e}   "
              f"B = {c2[1]:.6e} ± {s2[1]:.2e}")
        print(f"    cubic: A = {c3[0]:.6e} ± {s3[0]:.2e}   "
              f"B = {c3[1]:.6e} ± {s3[1]:.2e}   "
              f"C = {c3[2]:.6e} ± {s3[2]:.2e}")

    # turn off any unused panels
    for ax in axes_flat[n:]:
        ax.axis("off")

    # shared axis labels
    fig.supxlabel(r"$\Omega_d M$", y=0.01, fontsize=10)
    fig.supylabel(r"$|\mathcal{Q}_h|\,/\,\Omega_d$", x=0.005, fontsize=10)

    # single shared legend
    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right",
               bbox_to_anchor=(0.995, 0.998), ncol=2)

    fig.tight_layout(rect=(0.02, 0.02, 1, 0.98))
    fig.savefig(OUT_PATH)
    print(f"\nSaved: {OUT_PATH}")


if __name__ == "__main__":
    main()
