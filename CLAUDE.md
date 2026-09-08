# BlackHoleNonlinearities.jl — cluster branch

Second-order perturbation theory of a Schwarzschild black hole: given two linear
gravitational-wave modes, compute the quadratic coupling amplitude Q of the mode
they source. Extends Cardoso, Redondo-Yuste, Sperhake, Tuncer (2603.04501);
builds on Bucciotti, Juliano, Kuntz, Trincherini (2406.14611).

This branch is the compute-only subset for Della. The Mathematica sources
(`mma/`), the mma→Julia translator (`tools/`), analysis scripts (`plots/`) and
the convergence studies live on `main`.

## Conventions

- M = 1, so the horizon is at r = 2. `r` is areal, not tortoise.
- Parity: odd = axial = Regge-Wheeler; even = polar = Zerilli. l ≥ 2.
- Sector `XYZ` = parity of mode 1, mode 2, output. Six of them: eee eeo eoe eoo ooe ooo.
- Modes are `(l, m, ω)`; the output multipole `l` is passed separately, ω_out = ω1 + ω2.
- `Q_out = A2/(aout1*aout2)` — phase independent of the r* origin. `Q_in` uses the
  ingoing amplitudes and is not, with r* = r + 2log(r/2 − 1).
- Coupling selection rules come from Wigner 3j symbols, so forbidden combinations
  return zero rather than erroring.

## Architecture

    jl/scan_modes.jl        batch driver: one sector, threaded over a mode list (SLURM entry point)
    jl/runs.jl              interactive CLI: source / qfactor / qscan / qscan2d
      └ jl/QFactor.jl       qfactor(), qfactor_full(), qfactor_row_full(); Green's-function quadrature
        ├ jl/HomogeneousSolutions.jl   linear RW/Zerilli ODE solves, extract_amps()
        └ jl/Source.jl      dispatches coefficients()/source()/make_source() by sector
          └ jl/generated/SourceXYZ.jl  GENERATED — never hand-edit; regenerate on main

`slurm/run_sector.slurm` is a job array, one task per sector, reading mode
combinations from `slurm/modes.txt`.

## Accuracy caveats

- Q still depends on the free regularization constants at the 1–2.5% level
  (a horizon boundary term missing from the regulator). `test/test_qfactor.jl`
  tolerance is set to 2.5e-2 for this reason.
- The overall sign of Q (a π offset in its phase) is unsettled — the code uses
  +i/(2ω a_in) where the ±S convention would give −i.
- Per-sector `QUAD_RMIN_DEFAULT` in QFactor.jl is tuned and should not be changed
  casually: eoe needs 5, eee/eeo must stay at 3.

## Running

    julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'   # once, before submitting
    julia --project=. test/test_qfactor.jl                                  # smoke test, expect FINAL: PASS
    sbatch slurm/run_sector.slurm                                           # all six sectors
