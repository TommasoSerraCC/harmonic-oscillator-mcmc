# Harmonic Oscillator MCMC

Thermodynamics and spectrum of the one-dimensional quantum harmonic
oscillator, obtained from a Markov Chain Monte Carlo simulation of the
Euclidean path integral.

The project has three layers:

* a **Fortran core** that generates the Markov chain and measures the
  observables along each sampled path;
* a **NumPy analysis layer** that turns the raw chains into physical
  results, with blocking and jackknife error estimates;
* a **Tkinter GUI** that drives both and provides interactive fitting.

The physics results, including the extraction of the first three energy
gaps and the reconstruction of the ground-state wave function, are
collected in the accompanying report.

## What the simulation measures

| Quantity | Method | Exact value |
|---|---|---|
| $\langle E\rangle$, $\langle y^2\rangle$ | direct estimators extrapolated to $\eta\to0$ | $\tfrac12\coth(\beta\hbar\omega/2)$ |
| $E_1-E_0$ | plateau of the effective gap from $C_y(n)$ | 1 |
| $E_2-E_0$ | plateau of the effective gap from $C_{y^2}(n)$ | 2 |
| $E_3-E_0$ | plateau of the effective gap from $C_A(n)$ | 3 |
| $\lvert\psi_0(y)\rvert^2$ | histogram of sampled positions | $e^{-y^2}/\sqrt{\pi}$ |

```{toctree}
:maxdepth: 2
:caption: Contents

installation
physics
algorithm
usage
dataformats
fortran
api
testing
```

## Indices

* {ref}`genindex`
* {ref}`modindex`
