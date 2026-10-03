# BlackHoleNonlinearities.jl

Numerical calculation of quadratic gravitational-wave scattering by a Schwarzschild black hole. The code solves the Regge-Wheeler and Zerilli equations, evaluates the regularized second-order sources, and computes the quadratic coupling coefficient \(\mathcal Q\).

The conventions and results correspond to:

> B. Bucciotti, J. Redondo-Yuste, Z. Zhong, A. Kuntz, and V. Cardoso, “Second-order scattering response of a Schwarzschild black hole,” [arXiv:2609.20934](https://arxiv.org/abs/2609.20934).
> V. Cardoso, J. Redondo-Yuste, U. Sperhake, and F. Tuncer, “Nonlinear Dynamics in General Relativity,” [arXiv:2603.04501](https://arxiv.org/abs/2603.04501).

## Contents

- `mma/`: Mathematica notebooks and source terms.
- `jl/`: Solver for RW-Z homogeneous equations, evaluation of regularized source terms, and the Green-function integral for \(\mathcal Q\).
- `examples/`: small Julia examples that recompute coarse versions of the figures, and a tutorial notebook.
- `figures/`: the data behind every figure of the paper and the script that draws them.

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

`Qin` is the coupling coefficient of the paper, \(\mathcal Q = G_N M\,\mathcal A^{(2)}_{\rm out}/(\mathcal A^{(1)}_{\rm in,1}\mathcal A^{(1)}_{\rm in,2})\), normalized by the incoming strain amplitudes of the parent waves; `Qout` uses their outgoing amplitudes instead.

Conventions (\(M=1\), \(r_* = r + 2\log(r/2-1)\)):

- Strain amplitudes, \(\mathcal A = \tfrac{\mu\lambda}{2}(A_+ \mp i A_-)\): an odd (axial) master amplitude enters the outgoing strain with \(-i\) and the incoming strain with \(+i\).
- The angular coupling is \((-1)^{m'}\) times the 3j symbol \(\begin{pmatrix}\ell_1&\ell_2&\ell'\\ m_1&m_2&-m'\end{pmatrix}\), with \(m'=m_1+m_2\).
- The source is the exchange-summed bilinear kernel. When both parents have the same \((\ell, m, \omega)\) the result includes the factor \(1/2\) of a single driving mode; pass `symmetry_factor=false` to obtain the kernel itself (e.g. for maps that are continuous across \(\omega_1=\omega_2\)).

## Examples

Run the compact figure examples individually or in one Julia process:

```sh
julia --project=. examples/low_frequency.jl
julia --project=. examples/resonance.jl
julia --project=. examples/run_all.jl
```

The scripts write PDFs to `examples/output/`. 

## Figures of the paper

`figures/data/` holds the values plotted in every figure of the paper, in the conventions above (one table per curve or map, with headers describing the columns). Redraw all figures with

```sh
python figures/make_figures.py   # numpy + matplotlib; writes figures/output/
```

The time-domain points in the self-coupling figure come from an independent time-domain evolution by Zhen Zhong. In the azimuthal-dependence figure, the few cells where this code does not reach 1% accuracy are filled with Zhen Zhong's independent frequency-domain values (flagged `source=1`; the two agree to better than 1% on every common cell).


## Citation

```bibtex
@article{Bucciotti:2026xyz,
  author = {Bucciotti, Bruno and Redondo-Yuste, Jaime and Zhong, Zhen and Kuntz, Adrien and Cardoso, Vitor},
  title = {Second-order scattering response of a Schwarzschild black hole},
  eprint = {2609.20934},
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
