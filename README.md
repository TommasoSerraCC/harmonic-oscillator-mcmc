# Quantum Harmonic Oscillator — Path Integral MCMC

[![CI](https://github.com/TommasoSerraCC/harmonic-oscillator-mcmc/actions/workflows/ci.yml/badge.svg)](https://github.com/TommasoSerraCC/harmonic-oscillator-mcmc/actions/workflows/ci.yml)
[![Documentation](https://github.com/TommasoSerraCC/harmonic-oscillator-mcmc/actions/workflows/docs.yml/badge.svg)](https://tommasoserracc.github.io/harmonic-oscillator-mcmc/)

Markov Chain Monte Carlo simulation of the one-dimensional quantum harmonic
oscillator in the Euclidean path integral formalism.

After a Wick rotation, discretising imaginary time into `N_t` slices turns the
quantum problem into a chain of coupled variables that can be sampled like an
ordinary classical system. From the sampled paths the code measures thermal
averages, energy gaps and the ground-state position distribution.

## How it works

Three layers:

- **Fortran core** — generates the chain with a heat-bath plus over-relaxation
  update, both with acceptance 1, and measures observables and correlators along
  each path.
- **Python analysis** — blocking and jackknife error estimates, autocorrelation
  times, connected correlators and effective energy gaps.
- **Tkinter GUI** — drives both and provides interactive plotting and fitting.

## Requirements

`gfortran`, GNU `make`, Python 3.9 or newer.

```sh
pip install -r requirements.txt
```

The GUI needs Tkinter, which ships with CPython on Windows and macOS; on Debian
or Ubuntu install `python3-tk`. The analysis layer does not need it.

## Build

```sh
make            # build the simulation executable
make test       # Fortran unit tests
make pytest     # Python test suite
make docs       # build the documentation
make ui         # launch the GUI
make help       # list all targets
```

## Usage

The GUI is the easiest way in:

```sh
make ui
```

It runs the simulation, analyses the output, and lets you plot and fit the
results interactively, saving each fit alongside its parameters.

Each stage also works on its own. The simulation reads its parameters as a
Fortran namelist on standard input:

```sh
echo "&params bhw=10, nsteps=1000000, n_nt=1, nt_vals(1)=100, energy_only=0 /" | ./main
```

`bhw` is the inverse temperature in units of `ħω`, `nsteps` the number of
measurements and `nt_vals` the lattice sizes to simulate. Raw chains are written
to `data/`; the analysis reads them and writes result tables to `results/`:

```sh
python python_scripts/analyze_and_save.py --bhw 10 --nt 100 --skip 50000
```

The estimators are importable on their own if you only need the numerics:

```python
from analysis import statistics as st

sigma = st.blocking_sigma(series, block_size=1000)
```

## Layout

```
main.f              simulation driver
mcmc/               algorithm, observables, random number generator
analysis/           blocking, jackknife, autocorrelation, analysis pipeline
python_scripts/     command line entry points
ui/                 Tkinter interface
tests/              Python test suite
test/               Fortran unit tests
docs/               documentation sources
```

## Tests

```sh
make test      # Fortran: RNG, action invariance under over-relaxation, equilibrium
make pytest    # Python: estimators checked against analytically known answers
```

Both run in CI on every push.

## Documentation

<https://tommasoserracc.github.io/harmonic-oscillator-mcmc/>

## License

[GPL-3.0](LICENSE)
