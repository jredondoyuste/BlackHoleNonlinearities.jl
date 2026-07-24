"""
Compute Q_h on a 2D (ω1, ω2) grid for nonlinear QNM channels.

ω1 and ω2 are independent — unlike compute_q.py which fixes ω1=ω2.
For each grid point the Q-factor is computed at two integration-domain
scales (lambda1, lambda2); err = ||Q1| - |Q2|| is a convergence estimate.

Output files: ../data/q2d_{parity}_{l1}{m1}_{l2}{m2}_{l}.dat
Columns:
  w1  w2  |Q1|/w_out  |Q2|/w_out  err/w_out
where w_out = w1 + w2.

The file can be loaded and reshaped as:
  data = np.loadtxt(file, comments="#")
  grid = data.reshape(N1, N2, -1)  # axis 0: ω1, axis 1: ω2
"""

import subprocess
import os
import sys
import numpy as np

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

PARITY    = "ooe"      # "ooo", "ooe", ...
SOURCE    = "bruno"

# ω1 axis
W1_MIN, W1_MAX, N1 = 0.05, 0.50, 20
# ω2 axis
W2_MIN, W2_MAX, N2 = 0.05, 0.50, 20

LAMBDA1   = 3.7        # first  integration-domain scale
LAMBDA2   = 4.0        # second integration-domain scale

# Channels: (l1, l2, l, m1, m2)
CHANNELS = [
    (2, 3, 3,  2, -1),
    (2, 4, 2,  2, -2),
]

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.join(SCRIPT_DIR, "..")
RUN_SCRIPT   = os.path.join(SCRIPT_DIR, "runs.jl")
DATA_DIR     = os.path.join(PROJECT_ROOT, "data")

# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------

def run_channel(l1, l2, l, m1, m2):
    outfile = os.path.join(DATA_DIR, f"q2d_{PARITY}_{l1}{m1}_{l2}{m2}_{l}.dat")

    cmd = [
        "julia", f"--project={PROJECT_ROOT}",
        RUN_SCRIPT,
        "qscan2d", PARITY,
        str(l1), str(l2), str(l), str(m1), str(m2),
        str(W1_MIN), str(W1_MAX), str(N1),
        str(W2_MIN), str(W2_MAX), str(N2),
        str(LAMBDA1), str(LAMBDA2),
        SOURCE,
    ]

    print(f"  Running {PARITY} 2D: l1={l1} l2={l2} l={l} m1={m1} m2={m2} "
          f"grid={N1}×{N2}  λ=({LAMBDA1},{LAMBDA2}) ...", flush=True)
    # stderr is not captured so ProgressMeter output shows in the terminal
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=None, text=True)

    if result.returncode != 0:
        print(f"  ERROR for channel (l1={l1},l2={l2},l={l})", file=sys.stderr)
        return

    rows = []
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        w1, w2   = float(parts[0]), float(parts[1])
        abs_q1   = float(parts[2])
        abs_q2   = float(parts[3])
        err      = float(parts[4])
        w_out    = w1 + w2
        rows.append((w1, w2, abs_q1 / w_out, abs_q2 / w_out, err / w_out))

    os.makedirs(DATA_DIR, exist_ok=True)
    with open(outfile, "w") as f:
        f.write("# w1  w2  |Q1|/w_out  |Q2|/w_out  err/w_out\n")
        f.write(f"# parity={PARITY} l1={l1} l2={l2} l={l} m1={m1} m2={m2} "
                f"source={SOURCE} lambda1={LAMBDA1} lambda2={LAMBDA2}\n")
        f.write(f"# grid: N1={N1} w1=[{W1_MIN},{W1_MAX}]  "
                f"N2={N2} w2=[{W2_MIN},{W2_MAX}]\n")
        f.write(f"# reshape: data.reshape({N1}, {N2}, -1)\n")
        for row in rows:
            f.write(("  ".join(f"{v:.6f}" if i < 2 else f"{v:.10e}"
                               for i, v in enumerate(row))) + "\n")


    print(f"  -> {outfile}  ({len(rows)} points,  {N1}×{N2} grid)")


def main():
    print(f"Q-factor 2D scan [{PARITY}]:  {N1}×{N2} grid")
    print(f"  ω1 ∈ [{W1_MIN}, {W1_MAX}],  ω2 ∈ [{W2_MIN}, {W2_MAX}]")
    print(f"  Source: {SOURCE},  λ = ({LAMBDA1}, {LAMBDA2})\n")

    for (l1, l2, l, m1, m2) in CHANNELS:
        run_channel(l1, l2, l, m1, m2)

    print("\nDone.")


if __name__ == "__main__":
    main()
