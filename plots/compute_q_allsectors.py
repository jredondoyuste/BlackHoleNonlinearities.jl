"""
Compute Q_h vs ω for a 3×3 grid of channels, one sector at a time.
Runs all channels in a SINGLE Julia process to avoid repeated JIT.

Usage:
    python3 compute_q_allsectors.py [sector] [npts] [lambda2]

sector: ooo, ooe, eoo, eoe, eeo, eee  (default: ooo)
npts:   number of frequency points (default: 10)

Output: ../data/q_{sector}_{l1}{m1}_{l2}{m2}_{l}.dat per channel
"""

import subprocess, os, sys, tempfile

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.join(SCRIPT_DIR, "..")
DATA_DIR     = os.path.join(PROJECT_ROOT, "data")
JULIA        = os.path.expanduser("~/.juliaup/bin/julia")

LAMBDA1   = 4.0
LAMBDA2   = 5.0
OMEGA_MIN = 0.1

# QNM fundamental frequencies (real part, M=1)
QNM_RE = {2: 0.37367, 3: 0.59944, 4: 0.80918, 5: 1.01233, 6: 1.21210, 7: 1.41118, 8: 1.60981}

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
    return str(m) if m <= 0 else str(m)


def main():
    sector = sys.argv[1] if len(sys.argv) > 1 else "ooo"
    npts   = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    global LAMBDA2
    LAMBDA2 = float(sys.argv[3]) if len(sys.argv) > 3 else LAMBDA2   # optional: second-domain scale
    assert sector in ("eee", "ooo", "ooe", "eoo", "eoe", "eeo"), f"Invalid sector: {sector}"

    channels = []
    for l1, m1, l2, m2 in BASE_CHANNELS:
        l = find_l(l1, l2, m1, m2, sector)
        if l is None:
            print(f"  ({l1},{m1:+d})×({l2},{m2:+d}): forbidden by parity selection rule — skipping", flush=True)
            continue
        if exchange_vanishes(l1, l2, l, sector):
            print(f"  ({l1},{m1:+d})×({l2},{m2:+d})→l={l}: vanishes by exchange symmetry (l1=l2, ω1=ω2, l odd) — skipping", flush=True)
            continue
        channels.append((l1, m1, l2, m2, l))

    if not channels:
        print("No valid channels.", flush=True)
        return

    # Build a Julia script that computes everything in one process
    n_channels = len(channels)
    julia_lines = [
        'include("jl/QFactor.jl")',
        'using ProgressMeter',
        '',
        'function lambda_domain(lam, w)',
        '    rmin = lam + 1.0',
        '    rmax = 10^(lam - 1) / (2 * w)',
        '    return rmin, rmax',
        'end',
        '',
        f'sector = "{sector}"',
        f'npts = {npts}',
        f'w_min = {OMEGA_MIN}',
        f'lam1 = {LAMBDA1}',
        f'lam2 = {LAMBDA2}',
        '',
        '# First qfactor call in this process pays a one-time JIT-compile',
        '# cost (roughly constant, ~20-25s, independent of sector) — this is',
        '# unavoidable startup overhead, not a sign anything is stuck.',
        'println(stderr, "[warming up: first call compiles the ODE/quadrature stack, ~20-25s]")',
        'flush(stderr)',
        't_warmup = @elapsed qfactor(sector, (2, 2, w_min), (2, -2, w_min), 2; rmin=5.0, rmax=100.0)',
        'println(stderr, "[warmup done in $(round(t_warmup,digits=1))s — remaining points are fast]\\n")',
        'flush(stderr)',
        '',
    ]

    for ci, (l1, m1, l2, m2, l) in enumerate(channels, start=1):
        m_out = m1 + m2
        w_max = omega_max_for_l(l)
        fname = f"q_{sector}_{l1}{format_m(m1)}_{l2}{format_m(m2)}_{l}.dat"
        fpath = os.path.join(DATA_DIR, fname)
        julia_lines.append(f'println(stderr, "[channel {ci}/{n_channels}] ({l1},{m1:+d})×({l2},{m2:+d})→({l},{m_out:+d})  ω∈[{OMEGA_MIN},{w_max:.4f}]")')
        julia_lines.append(f'flush(stderr)')
        julia_lines.append(f'open("{fpath}", "w") do f')
        julia_lines.append(f'    println(f, "# w  ReQout1 ImQout1 ReQout2 ImQout2 ReQin1 ImQin1 ReQin2 ImQin2 ReQcf ImQcf ReQrm ImQrm")')
        julia_lines.append(f'    println(f, "# parity={sector} l1={l1} l2={l2} l={l} m1={m1} m2={m2} lambda1={LAMBDA1} lambda2={LAMBDA2} rstar_origin: r*=r+2log(r/2-1) errcols: Qcf=cfree(c1m1=0.5+0.2im,c2m1=-0.3)@lambda1 Qrm=quad_rmin shifted by ±1 @lambda1")')
        julia_lines.append(f'    prog = Progress(npts; desc="  [channel {ci}/{n_channels}] w sweep: ", output=stderr, showspeed=true)')
        julia_lines.append(f'    for w in range(w_min, {w_max}, length=npts)')
        julia_lines.append(f'        rmin1, rmax1 = lambda_domain(lam1, w)')
        julia_lines.append(f'        rmin2, rmax2 = lambda_domain(lam2, w)')
        julia_lines.append(f'        try')
        julia_lines.append(f'            R1 = qfactor_full(sector, ({l1}, {m1}, w), ({l2}, {m2}, w), {l}; rmin=rmin1, rmax=rmax1)')
        julia_lines.append(f'            R2 = qfactor_full(sector, ({l1}, {m1}, w), ({l2}, {m2}, w), {l}; rmin=rmin2, rmax=rmax2)')
        julia_lines.append(f'            Qcf = NaN + NaN*im')
        julia_lines.append(f'            try')
        julia_lines.append(f'                Qcf = qfactor(sector, ({l1}, {m1}, w), ({l2}, {m2}, w), {l}; rmin=rmin1, rmax=rmax1, cfree=(c1m1=0.5+0.2im, c2m1=-0.3))')
        julia_lines.append(f'            catch e')
        julia_lines.append(f'                println(stderr, "    w=", w, " Qcf failed: ", e)')
        julia_lines.append(f'            end')
        julia_lines.append(f'            Qrm = NaN + NaN*im')
        julia_lines.append(f'            try')
        julia_lines.append(f'                Qrm = qfactor(sector, ({l1}, {m1}, w), ({l2}, {m2}, w), {l}; rmin=rmin1, rmax=rmax1, quad_rmin=QUAD_RMIN_DEFAULT[sector] + (sector == "eoe" ? -1 : 1))')
        julia_lines.append(f'            catch e')
        julia_lines.append(f'                println(stderr, "    w=", w, " Qrm failed: ", e)')
        julia_lines.append(f'            end')
        julia_lines.append(f'            println(f, w, "  ", real(R1.Qout), "  ", imag(R1.Qout), "  ", real(R2.Qout), "  ", imag(R2.Qout), "  ",')
        julia_lines.append(f'                     real(R1.Qin), "  ", imag(R1.Qin), "  ", real(R2.Qin), "  ", imag(R2.Qin), "  ",')
        julia_lines.append(f'                     real(Qcf), "  ", imag(Qcf), "  ", real(Qrm), "  ", imag(Qrm))')
        julia_lines.append(f'        catch e')
        julia_lines.append(f'            println(stderr, "    w=", w, " failed: ", e)')
        julia_lines.append(f'        end')
        julia_lines.append(f'        next!(prog)')
        julia_lines.append(f'    end')
        julia_lines.append(f'end')
        julia_lines.append(f'println(stderr, "  → {fname}\\n")')
        julia_lines.append(f'flush(stderr)')
        julia_lines.append('')

    julia_lines.append('println(stderr, "Done.")')

    os.makedirs(DATA_DIR, exist_ok=True)

    # Write temp Julia script and run it in one process
    with tempfile.NamedTemporaryFile(mode='w', suffix='.jl', dir=PROJECT_ROOT,
                                     delete=False, prefix='_compute_') as tmp:
        tmp.write('\n'.join(julia_lines))
        tmp_path = tmp.name

    print(f"Computing Q-factor for sector {sector.upper()}", flush=True)
    print(f"Λ=({LAMBDA1},{LAMBDA2}), {npts} points", flush=True)
    print(f"ω_max per channel: 0.75 × ω_QNM(l)  [cutoff: 2ω ≤ 1.5 × ω_QNM]", flush=True)
    print(f"{len(channels)} channels in ONE Julia process\n", flush=True)

    try:
        # Stream stderr live line-by-line instead of buffering it until the
        # whole process exits — with subprocess.run(..., stderr=PIPE) you'd
        # see nothing at all until everything is already done.
        proc = subprocess.Popen(
            [JULIA, f"--project={PROJECT_ROOT}", tmp_path],
            text=True, stderr=subprocess.PIPE, bufsize=1,
        )
        for line in proc.stderr:
            sys.stderr.write(line)
            sys.stderr.flush()
        returncode = proc.wait()
        if returncode != 0:
            print(f"Julia exited with code {returncode}", file=sys.stderr)
    finally:
        os.unlink(tmp_path)


if __name__ == "__main__":
    main()
