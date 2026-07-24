"""
Compute Q_h vs ω for multiple parity channels via runs.jl.

Output files: ../data/q_{parity}_{l1}{m1}_{l2}{m2}_{l}.dat

Columns:
  w  |Q1|/w  |Q2|/w  err/w
where Q1/Q2 are computed at lambda1/lambda2 and err = ||Q1| - |Q2||.

Normalization (atoh conversion ψ → h) is applied inside Julia's qfactor function.
"""

import subprocess
import os
import sys

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

PARITY    = "ooe"          # "ooo", "ooe", ...
OMEGA_MIN = 0.0001
OMEGA_MAX = 0.1
N_PTS     = 100            # ~0.01 step
LAMBDA1   = 3.0            # first integration-domain scale
LAMBDA2   = 3.5            # second integration-domain scale
SOURCE    = "bruno"        # "adrien" or "bruno"

# Channels: (l1, l2, l, m1, m2)
# OOO requires l+l1+l2=odd, l≥2, triangle rule
# OOE requires l+l1+l2=even, same geometric rules
CHANNELS = [
    # (2, 2, 2,  2, -2),
    # (2, 2, 4,  2, -2),
    # (2, 2, 4,  2, 2),
    # (2, 3, 3,  2, -3),
    # (2, 3, 5,  2, -3),
    # (3, 3, 3,  3, -3),
    # (3, 3, 6,  3, -3),
    (2, 4, 2, 2, -2),
    (2, 4, 4, 2, -2),
    # (2, 4, 6, 2, -2)
]

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.join(SCRIPT_DIR, "..")
RUN_SCRIPT   = os.path.join(SCRIPT_DIR, "runs.jl")
DATA_DIR     = os.path.join(PROJECT_ROOT, "data_PM")

# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------

def run_channel(l1, l2, l, m1, m2):
    outfile = os.path.join(DATA_DIR, f"q_{PARITY}_{l1}{m1}_{l2}{m2}_{l}.dat")

    cmd = [
        "julia", f"--project={PROJECT_ROOT}",
        RUN_SCRIPT,
        "qscan", PARITY,
        str(l1), str(l2), str(l), str(m1), str(m2),
        str(OMEGA_MIN), str(OMEGA_MAX), str(N_PTS),
        str(LAMBDA1), str(LAMBDA2),
        SOURCE,
    ]

    print(f"  Running {PARITY}: l1={l1} l2={l2} l={l} m1={m1} m2={m2} "
          f"λ=({LAMBDA1},{LAMBDA2}) ...", flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ERROR for channel (l1={l1},l2={l2},l={l}):", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        return

    rows = []
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 4:
            continue
        w      = float(parts[0])
        abs_q1 = float(parts[1])
        abs_q2 = float(parts[2])
        err    = float(parts[3])
        rows.append((w, abs_q1 / w, abs_q2 / w, err / w))

    os.makedirs(DATA_DIR, exist_ok=True)
    with open(outfile, "w") as f:
        f.write("# w  |Q1|/w  |Q2|/w  err/w\n")
        f.write(f"# parity={PARITY} l1={l1} l2={l2} l={l} m1={m1} m2={m2} "
                f"source={SOURCE} lambda1={LAMBDA1} lambda2={LAMBDA2}\n")
        for row in rows:
            f.write(("  ".join(f"{v:.10e}" if i > 0 else f"{v:.6f}"
                               for i, v in enumerate(row))) + "\n")

    print(f"  -> {outfile}  ({len(rows)} points)")


def main():
    print(f"Q-factor scan [{PARITY}]: {N_PTS} points in ω ∈ [{OMEGA_MIN}, {OMEGA_MAX}]")
    print(f"Source: {SOURCE},  λ = ({LAMBDA1}, {LAMBDA2})\n")

    for (l1, l2, l, m1, m2) in CHANNELS:
        run_channel(l1, l2, l, m1, m2)

    print("\nDone.")


if __name__ == "__main__":
    main()
