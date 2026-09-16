# BlackHoleNonlinearities.jl

Numerical calculation of quadratic gravitational-wave scattering by a Schwarzschild black hole. The code solves the Regge-Wheeler and Zerilli equations, evaluates the regularized second-order sources, and computes the quadratic coupling coefficient \(\mathcal Q\).

The conventions and results correspond to:

> B. Bucciotti, J. Redondo-Yuste, A. Kuntz, and V. Cardoso, “Second-order scattering response of a Schwarzschild BH,” [arXiv:2609.xxxxx](https://arxiv.org/abs/2609.xxxxx).
> V. Cardoso, J. Redondo-Yuste, U. Sperhake, and F. Tuncer, “Nonlinear Dynamics in General Relativity,” [arXiv:2603.04501](https://arxiv.org/abs/2603.04501).

## Contents

- `mma/`: Mathematica notebooks and source terms.
- `jl/`: Solver for RW-Z homogeneous equations, evaluation of regularized source terms, and the Green-function integral for \(\mathcal Q\).
- `examples/`: Julia examples reproducing the figures in the paper and a tutorial notebook.

## Installation

Install Julia 1.12, clone the repository, and instantiate the environment:

```sh
git clone https://github.com/jredondoyuste/BlackHoleNonlinearities.jl.git
cd BlackHoleNonlinearities.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Computing a coupling coefficient

```julia
include("jl/QFactor.jl")

mode1 = (2, 2, 0.15) #l1, m1, w1
mode2 = (2, -2, 0.15) #l2, m2, w2
result = qfactor_full("eee", mode1, mode2, 2; rmin=5.0, rmax=2000.0) #the lower integration limit is 2M(1+10^{-rmin})

result.Qin
result.Qout
```

`Qin` is normalized by the incoming amplitudes of the parent waves; `Qout` uses their outgoing amplitudes. 

## Examples

Run the compact figure examples individually or in one Julia process:

```sh
julia --project=. examples/low_frequency.jl
julia --project=. examples/resonance.jl
julia --project=. examples/run_all.jl
```

The scripts write PDFs to `examples/output/`. 


## Citation

```bibtex
@article{Bucciotti:2026xyz,
  author = {Bucciotti, Bruno and Redondo-Yuste, Jaime and Kuntz, Adrien and Cardoso, Vitor},
  title = {Second-order scattering response of a Schwarzschild black hole},
  eprint = {2609.XXXXX},
  archivePrefix = {arXiv},
  primaryClass = {gr-qc},
  year = {2026}
}

@article{Cardoso:2026llh,
  author = {Cardoso, Vitor and Redondo-Yuste, Jaime and Sperhake, Ulrich and Tuncer, Furkan},
  title = {Nonlinear Dynamics in General Relativity},
  eprint = {2603.04501},
  archivePrefix = {arXiv},
  primaryClass = {gr-qc},
  year = {2026}
}
```

The source derivation builds on [QuadraticQNM](https://github.com/akuntz00/QuadraticQNM) by Adrien Kuntz and Bruno Bucciotti. Please also cite:

```bibtex
@article{Bucciotti:2024jrv,
  author = {Bucciotti, Bruno and Juliano, Leonardo and Kuntz, Adrien and Trincherini, Enrico},
  title = {Amplitudes and polarizations of quadratic quasi-normal modes for a Schwarzschild black hole},
  eprint = {2406.14611},
  archivePrefix = {arXiv},
  primaryClass = {hep-th},
  doi = {10.1007/JHEP09(2024)119},
  journal = {JHEP},
  volume = {09},
  pages = {119},
  year = {2024}
}
```

## License

See [LICENSE](LICENSE).
