# Quantum Harmonic Oscillator — Path Integral MCMC

[![CI](https://github.com/TommasoSerraCC/harmonic-oscillator-mcmc/actions/workflows/ci.yml/badge.svg)](https://github.com/TommasoSerraCC/harmonic-oscillator-mcmc/actions/workflows/ci.yml)
[![Documentation](https://github.com/TommasoSerraCC/harmonic-oscillator-mcmc/actions/workflows/docs.yml/badge.svg)](https://tommasoserracc.github.io/harmonic-oscillator-mcmc/)

Markov Chain Monte Carlo simulation of the one-dimensional quantum harmonic
oscillator in the Euclidean path integral formalism.

## Physics

After a Wick rotation to imaginary time, the quantum partition function becomes
a sum over periodic paths weighted by `e^{−S}`, with `S` the Euclidean action.
Discretising imaginary time into `N_t` slices of spacing `η = aω` turns the
oscillator into a one-dimensional chain of variables coupled to their nearest
neighbours — an ordinary classical statistical system that can be sampled by
Monte Carlo. Two limits connect the lattice back to the physics: `η → 0`
removes the discretisation error, and `βħω → ∞` projects onto the ground
state.

## Simulation

The sampling core is written in Fortran. Because the discretised action is
Gaussian, the conditional distribution of each site given its neighbours is
known in closed form, so the chain is updated with a **heat-bath** sweep: every
site is redrawn exactly from its conditional Gaussian (Box–Muller), with
acceptance 1 and no step-size tuning. Each heat-bath sweep is followed by five
**over-relaxation** sweeps — a deterministic reflection of each site about its
conditional mean that leaves the action invariant — which decorrelate
successive paths at negligible cost.

## Analysis

The Python layer estimates errors on the correlated Monte Carlo series by
**blocking**, scanning the block size until the estimated variance reaches a
plateau, and by **blocked jackknife** for nonlinear quantities such as
connected correlators and effective energy gaps. Autocorrelation functions and
integrated autocorrelation times are computed for every observable, and all
fits are weighted least squares with `χ²/ndf` and normalised residuals as
quality checks.

From the sampled ensemble the code extracts:

- the thermal moments `⟨y⟩`, `⟨y²⟩`, `⟨y³⟩`, with the odd ones serving as a
  symmetry check;
- the energy via the lattice virial estimator, whose temperature dependence
  reproduces `E(T) = 1/2 + 1/(e^{βħω} − 1)` and whose continuum extrapolation
  in `η²` recovers the exact ground-state energy `ħω/2`;
- the low-lying spectrum: effective gaps
  `ΔE(n) = (1/η) · log[C(n)/C(n+1)]` from the connected correlators of `y`,
  `y²` and `y³` show clean plateaux at `ħω`, `2ħω` and `3ħω`;
- the ground-state probability density `|ψ₀(y)|²`, reconstructed from the
  histogram of sampled positions and matching the exact Gaussian.

## License

[GPL-3.0](LICENSE)
