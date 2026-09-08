"""
Plot the convergence of the complex Q (modulus and phase) from
test/convergence_phase.jl: data/convphase_{sector}_*.dat with columns
    kind value w |Qout| argQout |Qin| argQin
Four columns of panels (rmax, quad_rmin, tolerance, cfree), two rows
(relative |Q| error, phase error in rad). One line per (sector, w); the
reference is the last point of each sweep (largest rmax, smallest rmin,
tightest tol, cfree=0).

Usage: python3 plot_convphase.py [out|in]
"""
import glob, os, sys
import numpy as np
import matplotlib.pyplot as plt

norm = sys.argv[1] if len(sys.argv) > 1 else "out"
col = 1 if norm == "out" else 3  # index into num (after kind,value)
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
KINDS = ["rmax", "rmin", "tol", "cfree"]
XL = {"rmax": r"$r_{\max}/M$", "rmin": r"$(r_{\min}-2M)/M$", "tol": "ODE/quad tol", "cfree": "cfree set"}
REF = {"rmax": -1, "rmin": -1, "tol": -1, "cfree": 0}

fig, ax = plt.subplots(2, 4, figsize=(14, 6))
for fn in sorted(glob.glob(os.path.join(D, "convphase_*.dat"))):
    sector = os.path.basename(fn).split("_")[1]
    rows = [l.split() for l in open(fn) if not l.startswith("#")]
    kinds = np.array([r[0] for r in rows]); vals = np.array([float(x) for r in rows for x in [r[1]]])
    num = np.array([[float(x) for x in r[2:]] for r in rows])
    w = num[:, 0]; Q = num[:, col] * np.exp(1j * num[:, col + 1])
    for j, k in enumerate(KINDS):
        for wi in np.unique(w):
            m = (kinds == k) & (w == wi)
            if m.sum() < 2: continue
            v, q = vals[m], Q[m]
            qref = q[REF[k]]
            dmod = np.abs(np.abs(q) - np.abs(qref)) / np.abs(qref)
            darg = np.abs(np.angle(q / qref))
            lab = f"{sector} w={wi:.3f}"
            if k == "cfree":
                ax[0, j].plot(v, dmod, "o-", ms=3, label=lab); ax[1, j].plot(v, darg, "o-", ms=3)
            else:
                ax[0, j].loglog(v, dmod + 1e-16, "o-", ms=3, label=lab); ax[1, j].loglog(v, darg + 1e-16, "o-", ms=3)
for j, k in enumerate(KINDS):
    ax[1, j].set_xlabel(XL[k]); ax[0, j].set_title(k)
ax[0, 0].set_ylabel(r"$||Q|-|Q_{\rm ref}||/|Q_{\rm ref}|$"); ax[1, 0].set_ylabel(r"$|\arg(Q/Q_{\rm ref})|$ [rad]")
ax[0, 0].legend(fontsize=5, ncol=2)
fig.suptitle(f"convergence of $Q_{{\\rm {norm}}}$"); fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"convphase_{norm}.pdf"); fig.savefig(out); print("Saved", out)
