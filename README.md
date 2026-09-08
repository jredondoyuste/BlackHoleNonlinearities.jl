## BlackHoleNonlinearities.jl 

Calculation of the excitation of higher harmonics when scattering gravitational waves off a Schwarzschild Black Hole. 
Code on constant development. 

> **`cluster` branch** — compute-only subset for HPC. `mma/`, `tools/`, `plots/`
> and the convergence studies are on `main`. See `CLAUDE.md`.

### Layout

| | |
|---|---|
| `mma/` | regularized source coefficients (Mathematica), input to the translator |
| `tools/` | `mma2julia.jl` translator and its Mathematica oracle |
| `jl/` | solver: `QFactor.jl`, `HomogeneousSolutions.jl`, `Source.jl` |
| `jl/generated/` | translated sources, one per sector — do not hand-edit |
| `jl/runs.jl` | CLI: `source` / `qfactor` / `qscan` / `qscan2d` |
| `jl/scan_modes.jl`, `slurm/` | batch driver and job array, one task per sector |
| `test/`, `plots/` | test suite and analysis scripts |

Sectors are named `XYZ` for the parities of mode 1, mode 2 and the output
(`e` = even/Zerilli, `o` = odd/Regge-Wheeler): eee, eeo, eoe, eoo, ooe, ooo.

    julia --project=. jl/runs.jl qfactor ooo 2 2 3 2 -2 0.3   # l1 l2 l m1 m2 omega
    julia --project=. test/test_qfactor.jl

[Code in development]

### Authors
Jaime Redondo-Yuste

### Citation

If you find this code useful, please cite

```
@article{Cardoso:2026llh,
    author = "Cardoso, Vitor and Redondo-Yuste, Jaime and Sperhake, Ulrich and Tuncer, Furkan",
    title = "{Nonlinear Dynamics in General Relativity}",
    eprint = "2603.04501",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "3",
    year = "2026"
}
```

This code builds upon the [QuadraticQNM](https://github.com/akuntz00/QuadraticQNM) package by Adrien Kuntz and Bruno Bucciotti. 
Please cite also their work:
```
@article{Bucciotti:2024jrv,
    author = "Bucciotti, Bruno and Juliano, Leonardo and Kuntz, Adrien and Trincherini, Enrico",
    title = "{Amplitudes and polarizations of quadratic quasi-normal modes for a Schwarzschild black hole}",
    eprint = "2406.14611",
    archivePrefix = "arXiv",
    primaryClass = "hep-th",
    doi = "10.1007/JHEP09(2024)119",
    journal = "JHEP",
    volume = "09",
    pages = "119",
    year = "2024"
}
```
