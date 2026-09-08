"""
Plot |Q| and arg(Q) vs ω for all channels of a given sector (3×3 channel grid,
each cell split into a modulus panel and a phase panel).

Usage:
    python3 plot_q_allsectors.py [sector] [--norm in|out]

sector: ooo, ooe, eoo, eoe, eeo, eee  (default: ooo)
--norm: out (default, Q normalised by outgoing linear amplitudes) or
        in (normalised by ingoing amplitudes; arg depends on the r* origin,
        r* = r + 2log(r/2-1), M=1)

Reads ../data/q_{sector}_*.dat files (13-column ReQout/ImQout/ReQin/ImQin at
two lambda domains plus ReQcf/ImQcf/ReQrm/ImQrm error-budget columns at
lambda1; falls back to the 9-column format with domain-only error bars, and
to the old 4-column |Q|-only format, plotting modulus only).
To generate them, run:
    python3 compute_q_allsectors.py <sector> [npts]
"""

import argparse
import os, sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec

def rgb(r, g, b):
    return (r / 255, g / 255, b / 255)

COLORS = {
    "deepcerulean":  rgb(  0, 123, 167),
    "crimsonred":    rgb(200,  50,  70),
    "frostedteal":   rgb( 58, 183, 183),
    "sageleaf":      rgb(141, 182,   0),
    "midnightblue":  rgb(  0,  70, 140),
    "amberglow":     rgb(255, 160,  58),
    "arcticnavy":    rgb( 20,  40,  80),
}

DATA_COLOR = COLORS["deepcerulean"]
LIGHT_COLOR = COLORS["frostedteal"]
QNM_COLOR  = COLORS["crimsonred"]

mpl.rcParams.update({
    "font.family":         "serif",
    "font.size":           9,
    "axes.labelsize":      9,
    "axes.titlesize":      9,
    "legend.fontsize":     7,
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
    "lines.linewidth":     0.7,
    "figure.dpi":          150,
    "savefig.dpi":         300,
    "savefig.bbox":        "tight",
    "savefig.pad_inches":  0.02,
})

QNM_RE = {2: 0.37367, 3: 0.59944, 4: 0.80918, 5: 1.01233, 6: 1.21210, 7: 1.41118, 8: 1.60981}

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(SCRIPT_DIR, "..", "data")

VALID_SECTORS = ("ooo", "ooe", "eoo", "eoe", "eeo", "eee")

BASE_CHANNELS = [
    (2,  2, 2, -2),
    (2,  1, 2, -1),
    (2,  2, 2,  1),
    (2,  2, 3, -2),
    (2,  2, 3, -1),
    (2,  1, 3, -1),
    (3,  3, 3, -3),
    (3,  2, 3, -2),
    (2,  2, 4, -2),
]


def find_l(l1, l2, m1, m2, sector):
    need_odd = sector in ("ooo", "eoe", "eeo")
    m = abs(m1 + m2)
    l_min = max(abs(l1 - l2), m, 2)
    for l in range(l_min, l1 + l2 + 1):
        if ((l + l1 + l2) % 2 == 1) == need_odd:
            return l
    return None


def exchange_vanishes(l1, l2, l, sector):
    """When ω1=ω2, l1=l2, same input parity, and l is odd, the result
    vanishes by exchange symmetry of the source and (-1)^l from the 3j."""
    return l1 == l2 and sector[0] == sector[1] and l % 2 == 1


def omega_max_for_l(l):
    """2ω ≤ 1.5 × ω_QNM(l) → ω ≤ 0.75 × ω_QNM(l)."""
    return 0.75 * QNM_RE[l]


def format_m(m):
    return f"+{m}" if m > 0 else str(m)


def channel_label(l1, m1, l2, m2, l):
    m = m1 + m2
    return (r"$(" + str(l1) + "," + format_m(m1) + r")\times("
            + str(l2) + "," + format_m(m2) + r")\to("
            + str(l) + "," + format_m(m) + r")$")


def empty_cell(fig, gs_cell, title, text):
    ax = fig.add_subplot(gs_cell[:])
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center", va="center",
            fontsize=9, color="0.5")
    ax.set_title(title, fontsize=10, pad=3)
    ax.set_xticks([])
    ax.set_yticks([])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sector", nargs="?", default="ooo")
    parser.add_argument("--norm", choices=["in", "out"], default="out")
    args = parser.parse_args()

    sector = args.sector
    assert sector in VALID_SECTORS, f"Invalid sector: {sector}"
    norm = args.norm

    fig = plt.figure(figsize=(11, 9.5))
    outer = gridspec.GridSpec(3, 3, figure=fig, wspace=0.32, hspace=0.42)

    for i, (l1, m1, l2, m2) in enumerate(BASE_CHANNELS):
        row, col = divmod(i, 3)
        cell = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[row, col], height_ratios=[3, 2], hspace=0.08)

        l = find_l(l1, l2, m1, m2, sector)
        if l is None:
            empty_cell(fig, cell, f"({l1},{format_m(m1)})$\\times$({l2},{format_m(m2)})",
                       "forbidden\n(parity selection rule)")
            continue

        if exchange_vanishes(l1, l2, l, sector):
            empty_cell(fig, cell, channel_label(l1, m1, l2, m2, l),
                       "$Q = 0$\n(exchange symmetry)")
            continue

        fname = f"q_{sector}_{l1}{m1 if m1 <= 0 else m1}_{l2}{m2 if m2 <= 0 else m2}_{l}.dat"
        fpath = os.path.join(DATA_DIR, fname)

        if not os.path.isfile(fpath):
            empty_cell(fig, cell, channel_label(l1, m1, l2, m2, l),
                       "no data\n(run compute_q_allsectors.py)")
            continue

        data = np.loadtxt(fpath, comments="#")
        if data.ndim < 2 or len(data) == 0:
            empty_cell(fig, cell, channel_label(l1, m1, l2, m2, l), "empty data")
            continue

        with open(fpath) as fh:
            header = fh.readline()

        w = data[:, 0]
        old_format = data.shape[1] < 9

        ax_mod = fig.add_subplot(cell[0])
        ax_ph  = fig.add_subplot(cell[1], sharex=ax_mod)

        if old_format:
            # Old format: # w |Q1| |Q2| err
            y1 = data[:, 1]
            y2 = data[:, 2]
            mask1 = (y1 > 0) & np.isfinite(y1)
            mask2 = (y2 > 0) & np.isfinite(y2)
            if mask2.sum() == 0 and mask1.sum() == 0:
                ax_mod.text(0.5, 0.5, "$Q = 0$", transform=ax_mod.transAxes,
                           ha="center", va="center", fontsize=9, color="0.5")
            else:
                ax_mod.plot(w[mask1], y1[mask1], "-", color=LIGHT_COLOR, lw=0.5, alpha=0.8, zorder=1)
                ax_mod.plot(w[mask2], y2[mask2], "o-", ms=2.5, color=DATA_COLOR, lw=0.6, zorder=2)
            ax_mod.set_ylim(bottom=0)
            ax_ph.text(0.5, 0.5, "phase n/a (old data)", transform=ax_ph.transAxes,
                      ha="center", va="center", fontsize=7, color="0.5")
            ax_ph.set_yticks([])
        else:
            # Qout-type quantities: always present at lambda1/lambda2 (cols 1-4);
            # error-budget columns (Qcf, Qrm) are Qout-type too, evaluated at
            # lambda1 only (cols 9-12), present only in the new format.
            Qout1 = data[:, 1] + 1j * data[:, 2]
            Qout2 = data[:, 3] + 1j * data[:, 4]
            Qin1  = data[:, 5] + 1j * data[:, 6]
            Qin2  = data[:, 7] + 1j * data[:, 8]
            has_err = data.shape[1] >= 13
            if has_err:
                Qcf = data[:, 9]  + 1j * data[:, 10]
                Qrm = data[:, 11] + 1j * data[:, 12]
            else:
                Qcf = np.full_like(Qout1, np.nan + 1j * np.nan)
                Qrm = np.full_like(Qout1, np.nan + 1j * np.nan)

            Qcentral = Qout2 if norm == "out" else Qin2

            # Error terms are always built from Qout-based ratios (domain,
            # cfree perturbation, quad_rmin shift): for --norm in this makes
            # them normalisation-independent up to the |Qin2|/|Qout2|
            # amplitude ratio, which is why the fractional error is applied
            # to |Qcentral| rather than reusing an absolute Qout difference.
            magout2 = np.abs(Qout2)
            np.seterr(divide="ignore", invalid="ignore")  # Q=0 channels give 0/0, masked below
            frac_mod = np.nanmax(np.stack([
                np.abs(np.abs(Qout2) - np.abs(Qout1)),
                np.abs(np.abs(Qout2) - np.abs(Qcf)),
                np.abs(np.abs(Qout2) - np.abs(Qrm)),
            ]), axis=0) / magout2
            yerr_mod = frac_mod * np.abs(Qcentral)

            yerr_arg = np.nanmax(np.stack([
                np.abs(np.angle(Qout2 / Qout1)),
                np.abs(np.angle(Qout2 / Qcf)),
                np.abs(np.angle(Qout2 / Qrm)),
            ]), axis=0) / np.pi

            mag2 = np.abs(Qcentral)
            mask2 = (mag2 > 0) & np.isfinite(mag2)

            if mask2.sum() == 0:
                ax_mod.text(0.5, 0.5, "$Q = 0$", transform=ax_mod.transAxes,
                           ha="center", va="center", fontsize=9, color="0.5")
                ax_ph.set_yticks([])
            else:
                ax_mod.errorbar(
                    w[mask2], mag2[mask2], yerr=yerr_mod[mask2],
                    fmt="o-", ms=2.5, mfc=DATA_COLOR, mec=DATA_COLOR,
                    color=DATA_COLOR, ecolor=DATA_COLOR, elinewidth=0.5,
                    capsize=1.5, lw=0.6, zorder=2,
                )
                ax_mod.set_ylim(bottom=0)

                ph2 = np.unwrap(np.angle(Qcentral[mask2])) / np.pi
                ax_ph.errorbar(
                    w[mask2], ph2, yerr=yerr_arg[mask2],
                    fmt="o-", ms=2.5, mfc=DATA_COLOR, mec=DATA_COLOR,
                    color=DATA_COLOR, ecolor=DATA_COLOR, elinewidth=0.5,
                    capsize=1.5, lw=0.6, zorder=2,
                )

        omega_qnm_half = QNM_RE.get(l, 0.5) / 2
        for ax in (ax_mod, ax_ph):
            ax.axvline(omega_qnm_half, color=QNM_COLOR, ls="--", lw=0.8,
                       alpha=0.7, zorder=1, label=r"$\omega_{\rm QNM}/2$")
            ax.set_xlim(0, omega_max_for_l(l))
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

        plt.setp(ax_mod.get_xticklabels(), visible=False)
        ax_mod.set_title(channel_label(l1, m1, l2, m2, l), fontsize=10, pad=3)
        ax_ph.set_xlabel(r"$\omega M$", fontsize=8)
        ax_mod.set_ylabel(r"$|\mathcal{Q}|$", fontsize=8)
        ax_ph.set_ylabel(r"$\arg\mathcal{Q}/\pi$", fontsize=8)

    handles, labels = [], []
    for ax in fig.axes:
        h, lab = ax.get_legend_handles_labels()
        if h:
            handles, labels = h[:1], lab[:1]
            break
    if handles:
        fig.legend(handles, labels, loc="upper right",
                   bbox_to_anchor=(0.98, 0.99), fontsize=8)

    norm_label = r"$\mathcal{Q}_{\rm out}$" if norm == "out" else r"$\mathcal{Q}_{\rm in}$"
    fig.suptitle(f"{sector.upper()} sector    " + r"($\omega_1 = \omega_2 = \omega$)   " + norm_label,
                 fontsize=12, y=1.0)

    out_path = os.path.join(SCRIPT_DIR, f"q_allfreq_{sector}.pdf")
    fig.savefig(out_path)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
