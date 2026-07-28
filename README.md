# Quantum Harmonic Oscillator — Path Integral MCMC

[![CI](https://github.com/TommasoSerraCC/harmonic-oscillator-mcmc/actions/workflows/ci.yml/badge.svg)](https://github.com/TommasoSerraCC/harmonic-oscillator-mcmc/actions/workflows/ci.yml)
[![Documentation](https://readthedocs.org/projects/harmonic-oscillator-mcmc/badge/?version=latest)](https://harmonic-oscillator-mcmc.readthedocs.io/en/latest/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

Thermodynamics and spectrum of the 1D quantum harmonic oscillator, studied with a
Markov Chain Monte Carlo simulation of the Euclidean path integral.

A Fortran core generates the Markov chain, a NumPy layer does the statistical
analysis (blocking, jackknife, autocorrelations), and a Tkinter GUI ties the two
together and provides interactive fitting. The full write-up is in
[`relazione/relazione.pdf`](relazione/relazione.pdf); `relazione/oscillator.pdf`
keeps the original Italian version.

Full documentation is built with Sphinx from [`docs/`](docs/) — see
[Documentation](#documentation) below.

---

## What it computes

After a Wick rotation the partition function `Z = Tr(exp(-βH))` becomes a functional
integral over **periodic** paths weighted by `exp(-S_E)`. Discretising Euclidean time
into `N_t` slices maps the quantum problem onto a 1D statistical system of `N_t` real
variables `y_i`, with action (in dimensionless units `η = aω`, `y = x·sqrt(mω)`):

```
S_E = Σ_i [ y_i² (η/2 + 1/η) − (1/η) y_i y_{i+1} ]
```

From the sampled paths the code extracts:

| Quantity | Method | Exact value |
|---|---|---|
| `⟨E⟩`, `⟨y²⟩` | direct estimators, extrapolated to `η → 0` | `½·coth(βℏω/2)` |
| `E₁ − E₀` | plateau of the effective gap from `C_y(n)` | 1 |
| `E₂ − E₀` | plateau of the effective gap from `C_{y²}(n)` | 2 |
| `E₃ − E₀` | plateau of the effective gap from `C_A(n)`, `A = y³ − 3y/2` | 3 |
| `\|ψ₀(y)\|²` | histogram of sampled positions | `exp(−y²)/√π` |

The Markov chain is **heat bath + over-relaxation**, both with acceptance exactly 1.
The heat bath samples each site directly from its Gaussian conditional
`P(y_i | y_{j≠i})`; the microcanonical (over-relaxation) sweep reflects each site
about that conditional mean, which leaves the action invariant and decorrelates
cheaply. One recorded measurement corresponds to `10 × (1 heat bath + 5 over-relaxation)`
sweeps, which keeps autocorrelation times of order 1.

---

## Repository layout

```
Makefile                  build, test, docs and report targets
main.f                    simulation driver, reads its parameters as a Fortran namelist
mcmc/
  mcmc_oscillator.f       heat bath, over-relaxation, Box-Muller, Euclidean action
  phys_sub.f              observables (y, y², y³, A), correlators, energy estimator
  rand.f                  ran2 RNG (Numerical Recipes) with seed persistence
test/test.f               unit tests for the Fortran core
analysis/
  statistics.py           blocking, jackknife and autocorrelation estimators
  pipeline.py             raw chain in, result tables out
  formatting.py           value ± uncertainty rounding
python_scripts/
  analyze_and_save.py     CLI wrapper around analysis.pipeline
  fit_correlation_length.py   exponential fits of the connected correlators
  plot_results.py         CLI plots: blocking, jackknife, gaps, autocorrelation
  plot_summary.py         CLI plots: observables vs η
  plot_correlators.py     CLI plots: the four connected correlators
ui/
  app.py                  Tkinter notebook with the seven analysis tabs
  core/data_manager.py    dataset discovery, file I/O, saved-fit archive
  core/plotting.py        embedded matplotlib canvas, shared widgets and styling
  core/ground_state_sim.py    vectorised NumPy re-implementation of the chain
  tabs/                   one module per tab
tests/                    pytest suite for the Python layer
docs/                     Sphinx documentation sources
run_ui.py                 GUI entry point
relazione/                LaTeX report and its PDF
results/                  post-processed analysis tables and saved fits (tracked)
plots/                    figures used by the report (tracked)
data/                     raw Markov chains (NOT tracked — regenerate locally)
```

---

## Requirements

- **gfortran** (any recent version; the code is fixed-form Fortran 77)
- **GNU make**
- **Python 3.9+** with the packages in `requirements.txt`:
  `pip install -r requirements.txt`
- **pdflatex** — only needed to rebuild the report

For the test suite and the documentation install `requirements-dev.txt` instead.

Tkinter ships with the standard CPython installer on Windows and macOS; on Debian
or Ubuntu install `python3-tk`. The analysis layer itself does not import Tkinter,
so the command-line workflow and the tests run without it.

---

## Build

```sh
make            # builds ./main (main.exe on Windows)
make test       # builds and runs the Fortran unit tests
make pytest     # runs the Python test suite
make docs       # builds the Sphinx documentation
make clean      # removes executables and LaTeX leftovers
make help       # lists all targets
```

---

## Quick start

The GUI is the intended way to drive the whole pipeline:

```sh
make ui         # or: python run_ui.py
```

The window has seven tabs (labels are in Italian):

| Tab | What it does |
|---|---|
| **Analisi** | Runs the Fortran simulation and the Python analysis; browses `data/` and `results/` |
| **Osservabili** | Means vs `η` with the continuum-limit fit; blocking/jackknife plateaus; autocorrelation and `τ_exp` |
| **Correlatori** | Connected correlators in log scale with exponential fits; multi-`N_t` comparison |
| **Gap Energetici** | Effective gap `ΔE(n)` with a weighted constant fit over the plateau |
| **E vs T** | Energy vs temperature at fixed `η`, fitted with `a + b/(exp(βℏω) − 1)` |
| **Istogramma** | Standalone ground-state simulation and position histogram vs `\|ψ₀\|²` |
| **Risultati** | File browser with text and image preview |

A typical session:

1. In **Analisi**, set `βℏω`, the number of steps and the list of `N_t`, then
   *Avvia acquisizione*. This writes `data/bhw{N}_nstep{M}/raw_data_nt{T}.dat`.
   Tick *Solo energia* to record only the energy — much smaller files, enough for
   the energy-vs-temperature scan.
2. Still in **Analisi**, set the thermalisation `skip` and press
   *Analizza tutti mancanti*. This runs `analyze_and_save.py` on every raw file
   that has no matching result yet.
3. Use the remaining tabs to plot and fit. Every fit reports its parameters,
   `χ²/ndof` and the parameter correlation, can be shown with normalised
   residuals, and can be archived with *Salva fit* (JSON + PNG under
   `results/.../fits/`, browsable from the side panel).

---

## Command-line workflow

The GUI is a front end; each stage can also be run on its own.

**1. Generate a chain.** `main` reads a Fortran namelist from stdin:

```sh
echo "&params bhw=10, nsteps=1000000, n_nt=3, nt_vals(1)=50, nt_vals(2)=100, nt_vals(3)=200, energy_only=0 /" | ./main
```

| Field | Meaning |
|---|---|
| `bhw` | `βℏω`, the inverse temperature in units of `ℏω` |
| `nsteps` | number of recorded measurements (each = 10 heat-bath sweeps) |
| `n_nt` | how many entries of `nt_vals` to use |
| `nt_vals(i)` | the lattice sizes `N_t` to simulate, one run each |
| `energy_only` | `1` records only the energy and skips the correlators |

**2. Analyse.**

```sh
python python_scripts/analyze_and_save.py --bhw 10 --nt 200 --nstep 1000000 --skip 50000
```

`--skip` is the number of leading measurements discarded as thermalisation; it
becomes part of the output directory name, so several choices can coexist.

**3. Plot.**

```sh
python python_scripts/plot_results.py     --bhw 10 --nt 200 --skip 50000 --save
python python_scripts/plot_correlators.py --bhw 10 --nt 200 --skip 50000 --save
python python_scripts/plot_summary.py     --bhw 10 --save
python python_scripts/fit_correlation_length.py --bhw 10 --nt 200 --skip 50000 --save
```

All scripts must be run from the repository root, since they use paths relative
to it.

---

## Data formats

### Raw chain — `data/bhw{N}_nstep{M}/raw_data_nt{T}.dat`

One line per recorded measurement, `5 + 4·(N_t/2)` columns:

```
y  y²  y³  A  E   then, for n = 1 … N_t/2:   y(0)y(n)  y²(0)y²(n)  y³(0)y³(n)  A(0)A(n)
```

Each entry is already averaged over the `N_t` sites of the path. In
`energy_only` mode the file is named `raw_energy_nt{T}.dat` and has a single
column. These files are large (gigabytes for the bigger runs) and are **not**
tracked by git.

### Analysis output — `results/bhw{N}_nstep{M}/nt{T}_therm{S}/`

| File | Contents |
|---|---|
| `observables.dat` | mean and blocking error of `y`, `y²`, `y³`, `A`, `E` |
| `blocking_observables.dat` | `σ(mean)` vs block size `k = 2⁰ … 2^kmax` — used to locate the plateau |
| `connected_correlators.dat` | `C_O(n) = ⟨O(0)O(n)⟩ − ⟨O⟩²` with jackknife errors |
| `blocking_jackknife_correlators.dat` | jackknife error vs block size at three representative `n` |
| `energy_gaps.dat` | `ΔE(n) = log[C(n)/C(n+1)]/η` with jackknife errors |
| `tau_int.dat`, `tau_exp_fit.dat` | integrated and exponential autocorrelation times |
| `correlation_length.dat` | `ξ` and `ΔE` from the exponential fits (written by `fit_correlation_length.py`) |
| `fits/` | archived interactive fits (JSON + PNG) |

The jackknife is applied directly to the ratio `log[C(n)/C(n+1)]`, not propagated
from the individual correlators, because numerator and denominator are strongly
correlated.

---

## Reproducing the report

```sh
make report     # runs pdflatex twice to resolve cross-references
```

The report includes figures from `plots/` and `results/`, both of which are
tracked, so it builds from a fresh clone without re-running any simulation.

---

## Tests

```sh
make test      # Fortran unit tests
make pytest    # Python test suite
```

Both run in CI on every push (`.github/workflows/ci.yml`), together with a
documentation build that treats warnings as errors.

The **Fortran suite** checks the uniform RNG (mean and variance), the Box-Muller
transform (mean, variance and kurtosis), periodic neighbour indexing, **exact
conservation of the Euclidean action by the over-relaxation sweep**, the
equilibrium values of `⟨y⟩` and `⟨y²⟩` for the heat bath and for the production
update, the energy estimator, and the consistency of the correlator machinery at
zero separation. It exits with status 1 if any check fails.

The **Python suite** validates the estimators against independently known
answers rather than against previous output:

- AR(1) series, whose autocorrelation `ρ^k` and integrated time `ρ/(1-ρ)` are
  known analytically, pin the blocking plateau and the autocorrelation estimator;
- the vectorised jackknife is cross-checked against an explicit leave-one-block-out
  loop, for both the connected correlator and the energy gap;
- exact special cases (unit blocks reduce to the standard error; a pure
  exponential correlator gives a flat gap `1/(ξη)`);
- the NumPy sampler must reproduce the ground state — mean 0, variance ½,
  vanishing odd moments, Gaussian kurtosis;
- a real analysis is written and read back through every `data_manager` loader,
  pinning the file-format contract between the pipeline and the GUI.

## Documentation

```sh
make docs      # output in docs/_build/html
```

The sources live in `docs/` and cover the physics, the algorithm, the usage of
both the GUI and the CLI, the data formats, a reference for the Fortran
routines, and the autodoc-generated Python API.

`.readthedocs.yaml` builds the same pages on Read the Docs. To publish them,
import the repository from the Read the Docs dashboard; the badge above assumes
the project keeps the default slug `harmonic-oscillator-mcmc`, so adjust the two
URLs if you pick a different one.

---

## Results

From `βℏω = 10` with `N_t` from 4 to 200, `10⁶` measurements per lattice size:

| Quantity | Measured | Exact |
|---|---|---|
| `E₀` (continuum extrapolation, linear in `η²`) | 0.4996(3) | 0.5 |
| `E₁ − E₀` | 0.9999(3) | 1 |
| `E₂ − E₀` | 2.0005(8) | 2 |
| `E₃ − E₀` | 2.999(2) | 3 |

The gap figures come from a naive constant fit over the plateau that ignores the
correlation between neighbouring `n`, so the quoted uncertainties are indicative
only; this is discussed in the report.

---

## Notes

- `main.f` creates its output directory with a Windows `mkdir` call. On Linux or
  macOS create `data/` beforehand, or adapt that one line.
- `randomseed` stores the RNG state between runs, so consecutive invocations
  continue the stream instead of repeating it. Delete it to start from the
  default seed.
- `plots/old_plots/` holds superseded figures from earlier versions of the
  analysis and is deliberately left untracked.

---

## License

See [LICENSE](LICENSE).
