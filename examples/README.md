# Examples

Each script recomputes a small version of a numerical figure:

- `low_frequency.jl`: low-frequency scaling in all six parity sectors.
- `resonance.jl`: response around the first three output QNM frequencies.
- `self_coupling.jl`: equal-mode self-couplings toward the eikonal regime.
- `angular_dependence.jl`: dependence on the parent azimuthal numbers.
- `parity_mixing.jl`: response over the complex parity-mixing plane.
- `unequal_frequencies.jl`: a coarse two-frequency response map.
- `convergence.jl`: sensitivity to the inner and outer integration endpoints.
- `tutorial.ipynb`: homogeneous solutions, source convergence, and the `Q` interface.

Run one example with `julia --project=. examples/NAME.jl`, or run all scripts in one process with `julia --project=. examples/run_all.jl`. Output is written to `examples/output/`.

These examples favor runtime over publication accuracy. Increase the grids, use the default `1e-10` tolerances, and enlarge the radial domain for production calculations.
