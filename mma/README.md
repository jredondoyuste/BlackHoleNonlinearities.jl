# Mathematica sources

These are Bruno Bucciotti's final derivation notebooks and exported expressions used for the calculation in arXiv:2603.04501.

- `EEE.nb` through `OOO.nb`: derivations for the six input/output parity sectors.
- `H.wl`, `Linear.wl`, `Source.wl`, `diffeoQ.wl`: shared definitions.
- `SourceXYZ.wl`: regularized sources.
- `CoefficientsXYZ.wl`: coefficient form consumed by the Julia implementation.

The generated Julia files in `../jl/generated/` are the numerical form of the six `CoefficientsXYZ.wl` files. The derivation builds on [QuadraticQNM](https://github.com/akuntz00/QuadraticQNM); cite Bucciotti et al., JHEP 09 (2024) 119 in addition to the paper associated with this repository.
