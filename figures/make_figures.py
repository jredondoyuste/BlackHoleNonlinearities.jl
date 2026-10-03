#!/usr/bin/env python3
"""Redraw the figures of the paper from figures/data (numpy + matplotlib).

    python figures/make_figures.py [--outdir DIR]

Every Q in figures/data is final: strain amplitudes with the phases of Eq. (19),
the (-1)^m' of Eq. (39), and the factor 1/2 for identical parents (l, m, omega),
exactly as returned by jl/QFactor.jl. The scripts only select and draw.
"""
from pathlib import Path
import argparse
import os
import re
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "bhn-mpl"))
import matplotlib
matplotlib.use("pdf")
import matplotlib.colors as mc
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import FixedFormatter, FixedLocator, MaxNLocator, NullFormatter

DATA = Path(__file__).resolve().parent / "data"
COLW = 4.6
TOL = 0.01  # pointwise accuracy criterion of the appendix

QNM_RE_OMEGA = {2: 0.37367, 3: 0.59943, 4: 0.80918, 5: 1.01230, 6: 1.21200, 7: 1.40970,
                8: 1.6061937282727146, 9: 1.8017947810838226, 10: 1.996787794118052,
                11: 2.191333780445272, 12: 2.3855413807894985}
SECTORS = ["eee", "eeo", "eoe", "eoo", "ooe", "ooo"]
SECTOR_COLOR = {"eee": "#C8102E", "eeo": "#2E5EAA", "eoe": "#0E8A6E",
                "eoo": "#D97706", "ooe": "#7B4EA8", "ooo": "#A8562A"}
SECTOR_MARKER = {"eee": "o", "eeo": "s", "eoe": "^", "eoo": "D", "ooe": "v", "ooo": "P"}
MASTER = {"e": "+", "o": "-"}
QLABEL = r"$|\mathcal{Q}|$"


def style():
    plt.rcParams.update({
        "savefig.bbox": "tight", "savefig.pad_inches": 0.08, "pdf.fonttype": 42,
        "font.family": "serif", "font.serif": ["DejaVu Serif"], "mathtext.fontset": "cm",
        "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 8.5, "legend.fontsize": 6.2,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.linewidth": 0.6, "axes.labelpad": 2.0,
        "lines.linewidth": 1.2, "lines.markersize": 3.0, "lines.markeredgewidth": 0.0,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.45,
        "ytick.minor.width": 0.45, "xtick.major.size": 3.0, "ytick.major.size": 3.0,
        "xtick.minor.size": 1.7, "ytick.minor.size": 1.7, "xtick.minor.visible": True,
        "ytick.minor.visible": True, "legend.frameon": False, "legend.handlelength": 1.6,
        "legend.handletextpad": 0.5, "legend.labelspacing": 0.28, "legend.borderaxespad": 0.4,
        "legend.columnspacing": 1.0, "axes.grid": False})


def header(path):
    meta = {}
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            meta.update(re.findall(r"(\w+)=(-?\w+)", line))
    return meta


class Channel:
    """One equal-frequency channel (omega1 = omega2) read from figures/data."""

    def __init__(self, path):
        meta = header(path)
        self.sector = meta["sector"]
        self.l1, self.m1, self.l2, self.m2, self.l = (int(meta[k]) for k in ("l1", "m1", "l2", "m2", "l"))
        self.representative = meta.get("representative") == "1"
        self.order = int(meta.get("order", 0))
        a = np.loadtxt(path, ndmin=2)
        self.w, self.Q = a[:, 0], a[:, 1] + 1j * a[:, 2]
        self.rel, self.gap = a[:, 5], a[:, 6].astype(bool)
        self.trusted = a[:, 7].astype(bool) if a.shape[1] > 7 else np.ones(len(self.w), bool)

    m = property(lambda self: self.m1 + self.m2)
    wp = property(lambda self: 2.0 * self.w)
    x = property(lambda self: self.wp / QNM_RE_OMEGA[self.l])
    absQ = property(lambda self: np.abs(self.Q))
    key = property(lambda self: (self.sector, self.l1, self.m1, self.l2, self.m2))

    @property
    def phase(self):
        ph = np.unwrap(np.angle(self.Q)) / np.pi
        return ph - 2.0 * np.round(ph[0] / 2.0) if len(ph) else ph

    def broken(self, values):
        values = np.array(values, dtype=float)
        values[self.gap] = np.nan
        return values

    def label(self):
        s = self.sector
        return (rf"$({self.l1},{self.m1},{MASTER[s[0]]})\!\times\!({self.l2},{self.m2},{MASTER[s[1]]})"
                rf"\!\to\!({self.l},{self.m},{MASTER[s[2]]})$")


def channels(folder):
    return [Channel(p) for p in sorted((DATA / folder).glob("*.dat"))]


def shade(hexcolor, f):
    r, g, b = mc.to_rgb(hexcolor)
    if f >= 0:
        return (r * (1 - f), g * (1 - f), b * (1 - f))
    return (r - (1 - r) * f, g - (1 - g) * f, b - (1 - b) * f)


def fig_lowfreq(out):
    x_subres = 0.85
    sel = sorted(channels("lowfreq"), key=lambda r: SECTORS.index(r.sector))
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(COLW, 3.8), sharex=True,
                                 gridspec_kw=dict(height_ratios=[1.8, 1.0], hspace=0.08))
    lo, hi = np.inf, 0.0
    for r in sel:
        k = r.x <= x_subres
        lo, hi = min(lo, r.wp[k].min()), max(hi, r.wp[k].max())
        me = max(1, int(k.sum()) // 14)
        kw = dict(marker=SECTOR_MARKER[r.sector], color=SECTOR_COLOR[r.sector], zorder=3, ms=2.6, markevery=me)
        ax.plot(r.wp[k], r.broken(r.absQ)[k], label=r.label(), **kw)
        bx.plot(r.wp[k], r.broken(r.phase)[k], **kw)
    ymin = min(r.absQ[r.x <= x_subres].min() for r in sel)
    ymax = max(r.absQ[r.x <= x_subres].max() for r in sel)
    smallest = min(sel, key=lambda r: r.absQ[r.x <= x_subres].max())
    k = smallest.x <= x_subres
    x0, y0 = float(smallest.wp[k][0]), float(smallest.absQ[k][0]) * 0.55
    xr = np.array([x0, min(2.2 * x0, 0.5 * hi)])
    ax.plot(xr, y0 * (xr / x0) ** 2, color="0.3", ls="--", lw=0.9, zorder=1)
    ax.annotate(r"$\propto(M\omega')^{2}$", xy=(xr[-1], y0 * (xr[-1] / x0) ** 2), xytext=(3.0, -1.0),
                textcoords="offset points", fontsize=6.8, color="0.3", ha="left", va="top")
    ax.set(xscale="log", yscale="log", ylabel=QLABEL, xlim=(0.8 * lo, 1.08 * hi), ylim=(0.3 * ymin, 130.0 * ymax))
    ticks = [t for t in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0) if 0.8 * lo <= t <= 1.55 * hi]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in ticks]))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.legend(loc="upper left", ncol=1, fontsize=5.4)
    bx.set_xlabel(r"$M\omega'$")
    bx.set_ylabel(r"$\arg\mathcal{Q}/\pi$")
    bx.set_ylim(-1.08, 1.08)
    bx.set_yticks([-1, -0.5, 0, 0.5, 1])
    bx.set_yticklabels([r"$-1$", r"$-\frac{1}{2}$", r"$0$", r"$\frac{1}{2}$", r"$1$"])
    fig.savefig(out / "results_lowfreq.pdf")
    plt.close(fig)


RESONANCE_COLORS = {2: "#C8102E", 3: "#B8860B", 4: "#2E5EAA", 5: "#0E8A6E", 6: "#7B4EA8", 7: "#D97706"}


def fig_resonance(out):
    groups = {}
    for r in channels("resonance"):
        groups.setdefault(r.l, []).append(r)
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(COLW, 3.7), sharex=True,
                                 gridspec_kw=dict(height_ratios=[1.7, 1.0], hspace=0.08))
    all_r = [r for chans in groups.values() for r in chans]
    xlo, xhi = 0.15, 1.02 * max(r.wp.max() for r in all_r)
    for l in sorted(groups):
        chans = sorted(groups[l], key=lambda r: r.order)
        shades = np.linspace(-0.25, 0.35, len(chans))
        for j, r in enumerate(chans):
            rep = r.representative
            kw = dict(color=shade(RESONANCE_COLORS[l], shades[j]), lw=1.75 if rep else 0.55,
                      alpha=1.0 if rep else 0.42, zorder=3 if rep else 2)
            ax.plot(r.wp, r.broken(r.absQ), label=rf"$\ell'={l}$" if j == 0 else None, **kw)
            bx.plot(r.wp, r.broken(r.phase), **kw)
    vis = [r.absQ[(r.wp >= xlo) & (r.wp <= xhi)] for r in all_r]
    ax.set_yscale("log")
    ax.set_ylim(0.85 * min(v.min() for v in vis if len(v)), 1.5 * max(v.max() for v in vis if len(v)))
    for l in sorted(groups):
        for a in (ax, bx):
            a.axvline(QNM_RE_OMEGA[l], color="0.45", ls=":", lw=0.7, zorder=0)
        ax.annotate(rf"$\ell\!=\!{l}$", xy=(QNM_RE_OMEGA[l], 0.0), xycoords=("data", "axes fraction"),
                    xytext=(1.8, 4.0), textcoords="offset points", fontsize=6.0, color="0.45",
                    rotation=90, ha="left", va="bottom")
    ax.set_ylabel(QLABEL)
    ax.set_xlim(xlo, xhi)
    ax.legend(loc="upper left", ncol=3, fontsize=6.0, handlelength=1.3, handletextpad=0.4,
              columnspacing=1.0, labelspacing=0.25, borderaxespad=0.1)
    bx.set_ylabel(r"$\arg\mathcal{Q}/\pi$")
    bx.set_xlabel(r"$M\omega'$")
    fig.subplots_adjust(top=0.96)
    fig.savefig(out / "results_resonance.pdf")
    plt.close(fig)


def fig_selfcoup(out, td_data=None):
    recs = sorted(channels("selfcoup"), key=lambda r: r.l1)
    shades = dict(zip([r.l1 for r in recs], np.linspace(-0.30, 0.45, len(recs))))
    fig, ax = plt.subplots(figsize=(COLW, 2.5))
    xmax = 0.0
    for r in recs:
        k = r.trusted
        xmax = max(xmax, r.x[k].max())
        ax.plot(r.x, np.where(k, r.absQ, np.nan), color=shade("#1F77B4", shades[r.l1]), lw=1.3,
                marker="s", ms=2.2, label=rf"$\ell={r.l1}\!\to\!{r.l}$", zorder=3)
    if td_data:
        w, re_q, im_q = np.loadtxt(td_data, unpack=True, usecols=(0, 1, 2))
        ax.plot(2 * w / QNM_RE_OMEGA[4], np.hypot(re_q, im_q), ls="none", marker="o", ms=3.0, mfc="none",
                mec="#C8102E", mew=0.7, label=r"$\ell=2\!\to\!4$, time domain", zorder=4)
    ax.axvline(1.0, color="0.5", ls=":", lw=0.7, zorder=0)
    ax.set_yscale("log")
    ax.set_ylim(bottom=3e-3)
    ax.set_xlim(right=xmax + 0.01)
    ax.set_xlabel(r"$\omega'/\mathrm{Re}[\omega_{\ell'0}]$")
    ax.set_ylabel(QLABEL)
    ax.legend(loc="upper left", fontsize=5.6, ncol=2, handlelength=1.3, columnspacing=0.9, labelspacing=0.2)
    fig.savefig(out / "results_selfcoup.pdf")
    plt.close(fig)


def fig_mdep(out):
    maps = {}
    for p in sorted((DATA / "mdep").glob("*.dat")):
        meta = header(p)
        a = np.loadtxt(p, ndmin=2)
        maps[(int(meta["l1"]), int(meta["l"]), meta["sector"])] = dict(
            m1=a[:, 0].astype(int), m2=a[:, 1].astype(int), Q=a[:, 2] + 1j * a[:, 3],
            ok=(a[:, 4] < TOL) | (a[:, 5] == 1))  # source 1: Zhen Zhong's values where this code fails 1%
    label = {"eee": r"(+,+)\to+", "eoe": r"(+,-)\to+"}
    rows = []
    for ell in (2, 3, 4, 5):
        row = []
        for ell_out in range(2, 2 * ell + 1):
            if ell == 5 and ell_out in (5, 7):
                continue
            sector = "eee" if ell_out % 2 == 0 else "eoe"
            rec = maps.get((ell, ell_out, sector))
            status = ("data\npending" if rec is None else
                      "ok" if rec["ok"].any() else "fails 1%\nvalidation")
            row.append(dict(rec=rec if status == "ok" else None, l1=ell, l=ell_out, sector=sector, status=status))
        rows.append(row)
    allq = np.concatenate([np.abs(r["Q"])[r["ok"]] for r in maps.values()])
    norm = LogNorm(vmin=max(allq.min(), allq.max() / 1e3), vmax=allq.max())
    fig = plt.figure(figsize=(7.0, 5.6))
    outer = fig.add_gridspec(len(rows), 1, left=0.085, right=0.89, bottom=0.045, top=0.975, hspace=0.015,
                             height_ratios=[1.0 / len(row) for row in rows])
    im = None
    for i, row in enumerate(rows):
        grid = outer[i].subgridspec(1, len(row), wspace=0.18)
        for j, slot in enumerate(row):
            ax = fig.add_subplot(grid[0, j])
            ax.set_title(rf"$\ell'={slot['l']},\ {label[slot['sector']]}$", fontsize=6.8, pad=1.2)
            r = slot["rec"]
            if j == 0:
                p = ax.get_position()
                fig.text(0.035, 0.5 * (p.y0 + p.y1), rf"$\ell_1=\ell_2={slot['l1']}$", ha="right",
                         va="center", rotation=90, fontsize=8.0, color="0.2")
            if r is None:
                ax.set_facecolor("0.975")
                for spine in ax.spines.values():
                    spine.set_color("0.72")
                    spine.set_linestyle(":")
                ax.text(0.5, 0.5, slot["status"], transform=ax.transAxes, ha="center", va="center",
                        fontsize=7.2, color="0.45")
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            ms = np.arange(-slot["l1"], slot["l1"] + 1)
            A = np.full((len(ms), len(ms)), np.nan)
            for a_, b_, q, ok in zip(r["m1"], r["m2"], np.abs(r["Q"]), r["ok"]):
                if ok:
                    A[b_ + slot["l1"], a_ + slot["l1"]] = q
            im = ax.pcolormesh(ms, ms, A, cmap="viridis", norm=norm, shading="nearest", rasterized=True)
            ax.set_aspect("equal")
            ax.set_xticks(ms[::max(1, len(ms) // 4)])
            ax.set_yticks(ms[::max(1, len(ms) // 4)])
            ax.tick_params(labelsize=6.6, pad=0.8)
            if j == 0:
                ax.set_ylabel(r"$m_2$", fontsize=8, labelpad=1)
            else:
                ax.tick_params(labelleft=False)
            if i == len(rows) - 1:
                ax.set_xlabel(r"$m_1$", fontsize=8, labelpad=0.5)
    cb = fig.colorbar(im, cax=fig.add_axes([0.925, 0.15, 0.022, 0.72]))
    cb.set_label(QLABEL, fontsize=8.0, labelpad=3.0)
    cb.ax.tick_params(labelsize=6.2, length=1.5)
    cb.outline.set_linewidth(0.4)
    fig.savefig(out / "results_mdep.pdf", bbox_inches=None)
    plt.close(fig)


def fig_parity_plane(out):
    coef = {}
    for line in (DATA / "parity_plane" / "sectors_22x22_to44.dat").read_text().splitlines():
        if not line.startswith("#"):
            name, _, re_q, im_q = line.split()
            coef[name] = float(re_q) + 1j * float(im_q)
    a, b, c = coef["a"], coef["b"], coef["c"]
    axis = np.linspace(-2.0, 2.0, 801)
    pre, pim = np.meshgrid(axis, axis)
    p = pre + 1j * pim
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        amplitude = np.abs((a + 2j * p * b - p**2 * c) / (1 + 1j * p) ** 2)
    amplitude[np.abs(1 + 1j * p) < 0.055] = np.nan  # pole at p = +i: no incident +m strain
    vmin, vmax = np.nanpercentile(amplitude[np.isfinite(amplitude) & (amplitude > 0)], [1, 99])
    vmin = max(vmin, 1e-4)
    with plt.style.context("default"), plt.rc_context({"font.size": 10, "axes.labelsize": 11,
                                                        "xtick.labelsize": 9.5, "ytick.labelsize": 9.5}):
        fig, ax = plt.subplots(figsize=(4.15, 3.45), constrained_layout=True)
        image = ax.pcolormesh(pre, pim, np.ma.masked_invalid(amplitude), shading="nearest", cmap="viridis",
                              norm=LogNorm(vmin=vmin, vmax=vmax), rasterized=True)
        ax.contour(pre, pim, amplitude, levels=np.geomspace(vmin, vmax, 7), colors="white",
                   linewidths=0.45, alpha=0.65)
        ax.plot(0, 0, "o", ms=5.5, color="white", mec="black", mew=0.8, label=r"$p=0$")
        ax.plot(0, 1, "*", ms=8, color="#d95f02", mec="black", mew=0.5, label=r"$p=\mathrm{i}$",
                clip_on=False, zorder=5)
        zeros = np.roots([c, -2j * b, -a])
        ax.plot(zeros.real, zeros.imag, "x", ms=5, mew=1.0, color="white", label=r"$\mathcal{Q}=0$", zorder=5)
        ax.set(xlabel=r"$\operatorname{Re}(p)$", ylabel=r"$\operatorname{Im}(p)$",
               xlim=(axis[0], axis[-1]), ylim=(axis[0], axis[-1]))
        ax.set_aspect("equal")
        ax.legend(loc="upper left", frameon=True, framealpha=0.82, fontsize=8)
        ticks = [v for v in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0) if vmin <= v <= vmax]
        cb = fig.colorbar(image, ax=ax, pad=0.025, ticks=ticks)
        cb.set_label(QLABEL)
        cb.ax.set_yticklabels([f"{v:g}" for v in ticks])
        fig.savefig(out / "results_parity_plane_22x22_to44.pdf")
        plt.close(fig)


def load_map(name, symmetric=False):
    a = np.loadtxt(DATA / "unequal" / f"eee_{name}.dat")
    xs, ys = np.unique(a[:, 0]), np.unique(a[:, 1])
    ix, iy = np.searchsorted(xs, a[:, 0]), np.searchsorted(ys, a[:, 1])
    q = a[:, 2] + 1j * a[:, 3]
    mag = np.full((len(ys), len(xs)), np.nan)
    phase = np.full_like(mag, np.nan)
    ok = a[:, 4] < TOL
    mag[iy[ok], ix[ok]], phase[iy[ok], ix[ok]] = np.abs(q[ok]), np.angle(q[ok]) / np.pi
    if symmetric:  # the 22x22 kernel is exchange symmetric; only one triangle is stored
        fill = ~np.isfinite(mag) & np.isfinite(mag.T)
        mag[fill], phase[fill] = mag.T[fill], phase.T[fill]
    return xs, ys, mag, phase


def fig_unequal(out):
    datasets = (load_map("22x22_to44", symmetric=True), load_map("22x33_to55"), load_map("22x44_to66"))
    headers = (r"$(2,2)\times(2,2)\to(4,4)$", r"$(2,2)\times(3,3)\to(5,5)$", r"$(2,2)\times(4,4)\to(6,6)$")
    xlabels = (r"$M\omega_{22}$", r"$M\omega_{33}$", r"$M\omega_{44}$")
    resonances = (0.8091783775322389, 1.012295312135351, 1.212009820652131)
    mags = np.concatenate([d[2][np.isfinite(d[2]) & (d[2] > 0)] for d in datasets])
    norm = LogNorm(float(mags.min()), float(mags.max()))
    ymin, ymax = max(datasets[0][1].min(), datasets[1][0].min(), datasets[2][0].min()), 0.51
    with plt.rc_context({"savefig.pad_inches": 0.02}):
        fig, axes = plt.subplots(2, 3, figsize=(7.05, 4.55), sharey=True,
                                 gridspec_kw={"hspace": 0.055, "wspace": 0.12})
        for col, ((xs, ys, mag, phase), head, res) in enumerate(zip(datasets, headers, resonances)):
            if col:  # vertical axis is the 22-parent frequency in every column
                xs, ys, mag, phase = ys, xs, mag.T, phase.T
            mag_im = axes[0, col].pcolormesh(xs, ys, np.ma.masked_invalid(mag), shading="nearest",
                                             cmap="viridis", norm=norm, rasterized=True)
            ph_im = axes[1, col].pcolormesh(xs, ys, np.ma.masked_invalid(phase), shading="nearest",
                                            cmap="twilight", vmin=-1, vmax=1, rasterized=True)
            axes[0, col].set_title(head, fontsize=10.5, pad=3)
            axes[0, col].tick_params(axis="x", bottom=False, labelbottom=False)
            axes[1, col].set_xlabel(xlabels[col], fontsize=10.5, labelpad=2)
            for row in range(2):
                ax = axes[row, col]
                ax.tick_params(axis="both", labelsize=9, pad=2)
                if col:
                    ax.tick_params(axis="y", labelleft=False)
                ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
                ax.set_yticks((0.1, 0.2, 0.3, 0.4, 0.5))
                ax.set_box_aspect(1)
                ax.set_ylim(ymin, ymax)
                xlim, ylim = ax.get_xlim(), ax.get_ylim()
                lo, hi = max(xlim[0], res - ylim[1]), min(xlim[1], res - ylim[0])
                if lo < hi:
                    line, = ax.plot([lo, hi], [res - lo, res - hi], color="#d73027", lw=0.75,
                                    ls=(0, (4, 2.5)), dash_capstyle="round", zorder=4)
                    line.set_path_effects([path_effects.Stroke(linewidth=1.35, foreground="white"),
                                           path_effects.Normal()])
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
        fig.supylabel(r"$M\omega_{22}$", fontsize=10.5, x=0.058)
        cb = fig.colorbar(mag_im, ax=list(axes[0, :]), pad=0.018, fraction=0.026, aspect=22)
        cb.set_label(QLABEL, fontsize=10.5)
        cb.ax.tick_params(labelsize=8.5)
        cb = fig.colorbar(ph_im, ax=list(axes[1, :]), pad=0.018, fraction=0.026, aspect=22, ticks=(-1, 0, 1))
        cb.set_label(r"$\arg\mathcal{Q}/\pi$", fontsize=10.5)
        cb.ax.tick_params(labelsize=8.5)
        fig.savefig(out / "results_unequal_combined.pdf")
        plt.close(fig)


def fig_accuracy(out):
    labels = {s: "$" + "".join(MASTER[c] for c in s) + "$" for s in SECTORS}
    fig, axes = plt.subplots(2, 1, figsize=(3.45, 4.25), gridspec_kw=dict(hspace=0.34))
    for i, (ax, which, xlabel) in enumerate(((axes[0], "inner", r"$(r_{\min}-2M)/M$"),
                                             (axes[1], "outer", r"$r_{\max}/M$"))):
        rows = []
        for path in sorted((DATA / "accuracy").glob(f"convergence_*_{which}.dat")):
            a = np.atleast_2d(np.loadtxt(path))
            rows.append((path.name.split("_")[1], a[:, 0], a[:, 1]))
        for sector, x, err in rows:
            ax.plot(x, err, color=SECTOR_COLOR[sector], marker=SECTOR_MARKER[sector], ms=2.7, lw=0.85,
                    label=labels[sector])
        ax.set(xscale="log", yscale="log", xlabel=xlabel)
        top = max(max(float(np.nanmax(e)) for _, _, e in rows) * 1.7, 3e-2)
        ax.set_ylim(min(float(np.nanmin(e[e > 0])) for _, _, e in rows) / 1.7, top)
        ax.axhspan(1e-2, top, color="#D55E00", alpha=0.09, linewidth=0, zorder=0)
        ax.axhline(1e-2, color="0.38", ls=":", lw=0.8, zorder=1)
        if i == 0:
            ax.set_ylabel("relative error")
        ax.tick_params(axis="both", labelsize=7.5)
        ax.set_box_aspect(0.62)
    axes[0].legend(loc="upper left", ncol=2, fontsize=5.5, handlelength=1.0, columnspacing=0.55,
                   handletextpad=0.25, borderpad=0.22, labelspacing=0.22, frameon=True,
                   facecolor="white", framealpha=0.88, edgecolor="none")
    fig.savefig(out / "results_accuracy.pdf")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parent / "output")
    parser.add_argument("--td-data", type=Path, default=DATA / "selfcoup_time_domain" / "eee_22_22_4_TD.dat",
                        help="time-domain table (w ReQ ImQ ...) overlaid on the self-coupling figure")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    style()
    fig_lowfreq(args.outdir)
    fig_resonance(args.outdir)
    fig_selfcoup(args.outdir, args.td_data)
    fig_mdep(args.outdir)
    fig_parity_plane(args.outdir)
    fig_unequal(args.outdir)
    fig_accuracy(args.outdir)
    print("wrote 7 figures to", args.outdir)


if __name__ == "__main__":
    main()
