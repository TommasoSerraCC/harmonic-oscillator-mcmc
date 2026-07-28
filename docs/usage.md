# Usage

The GUI is the intended way to drive the pipeline, but every stage can
be run on its own from the command line.

## The graphical interface

```sh
make ui        # or: python run_ui.py
```

The window is a notebook of seven tabs. The labels are in Italian, since
that is the language of the accompanying report.

| Tab | Purpose |
|---|---|
| **Analisi** | Launches the Fortran simulation and the Python analysis; browses `data/` and `results/` |
| **Osservabili** | Means against $\eta$ with the continuum-limit fit; blocking and jackknife plateaus; autocorrelation and $\tau_{\text{exp}}$ |
| **Correlatori** | Connected correlators on a logarithmic scale with exponential fits, and a multi-$N_t$ comparison |
| **Gap Energetici** | Effective gap $\Delta E(n)$ with a weighted constant fit over the plateau |
| **E vs T** | Energy against temperature at fixed $\eta$, fitted with $a + b/(e^{\beta\hbar\omega}-1)$ |
| **Istogramma** | Standalone ground-state simulation and position histogram against $\lvert\psi_0\rvert^2$ |
| **Risultati** | File browser with text and image preview |

### A typical session

1. In **Analisi**, set $\beta\hbar\omega$, the number of measurements and
   the list of $N_t$, then press *Avvia acquisizione*. Tick *Solo
   energia* to record only the energy, which produces far smaller files
   and is enough for the energy-versus-temperature scan.
2. Set the thermalisation `skip` and press *Analizza tutti mancanti*.
   Every raw file without a matching result is analysed.
3. Use the remaining tabs to plot and fit. Each fit reports its
   parameters, $\chi^2/\text{ndof}$ and the parameter correlation, can be
   displayed with normalised residuals, and can be archived with *Salva
   fit* as a JSON file plus a PNG under `results/.../fits/`.

The simulation and the analysis run on worker threads, so the interface
stays responsive; the ground-state simulation can also be cancelled
while it runs.

## Command line

### 1. Generate a chain

The executable reads a Fortran namelist from standard input:

```sh
echo "&params bhw=10, nsteps=1000000, n_nt=3, nt_vals(1)=50, nt_vals(2)=100, nt_vals(3)=200, energy_only=0 /" | ./main
```

| Field | Meaning |
|---|---|
| `bhw` | $\beta\hbar\omega$, the inverse temperature in units of $\hbar\omega$ |
| `nsteps` | number of recorded measurements |
| `n_nt` | how many entries of `nt_vals` to use |
| `nt_vals(i)` | the lattice sizes $N_t$ to simulate, one run each |
| `energy_only` | `1` records only the energy and skips the correlators |

Each run writes `data/bhw{N}_nstep{M}/raw_data_nt{T}.dat`.

The random seed is persisted in a file named `randomseed`, so
consecutive invocations continue the stream instead of repeating it.
Delete that file to start again from the default seed.

### 2. Analyse

```sh
python python_scripts/analyze_and_save.py \
    --bhw 10 --nt 200 --nstep 1000000 --skip 50000
```

`--skip` is the number of leading measurements discarded as
thermalisation. It becomes part of the output directory name, so several
choices can coexist and be compared.

### 3. Plot

```sh
python python_scripts/plot_results.py     --bhw 10 --nt 200 --skip 50000 --save
python python_scripts/plot_correlators.py --bhw 10 --nt 200 --skip 50000 --save
python python_scripts/plot_summary.py     --bhw 10 --save
python python_scripts/fit_correlation_length.py --bhw 10 --nt 200 --skip 50000 --save
```

All scripts must be run from the repository root, because they resolve
`data/`, `results/` and `plots/` relative to the working directory.

## Using the analysis as a library

The numerics are importable and free of file-system assumptions, so they
can be reused directly:

```python
from analysis import statistics as st
from analysis.pipeline import analyze_dataset

# error on a correlated series
sigma = st.blocking_sigma(series, block_size=1000)

# full analysis of one raw chain
summary = analyze_dataset('data/bhw10_nstep1000000/raw_data_nt200.dat',
                          'results/bhw10_nstep1000000/nt200_therm50000',
                          eta=0.05, skip=50000)
print(summary['means']['E'], summary['errors']['E'])
```
