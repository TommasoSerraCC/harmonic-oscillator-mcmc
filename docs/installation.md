# Installation

## Requirements

| Component | Needed for |
|---|---|
| `gfortran` | building the simulation and the Fortran tests |
| GNU `make` | the build targets |
| Python 3.9 or newer | the analysis layer and the GUI |
| `pdflatex` | rebuilding the report (optional) |

Tkinter is required by the GUI. It ships with the standard CPython
installer on Windows and macOS; on Debian or Ubuntu install
`python3-tk`. The analysis layer itself does not import Tkinter, so the
command line workflow and the test suite work without it.

## Setting up

```sh
git clone <repository-url>
cd harmonic-oscillator-mcmc

python -m venv .venv
# Windows:      .venv\Scripts\activate
# Linux/macOS:  source .venv/bin/activate

pip install -r requirements.txt
make
```

`requirements.txt` covers running the code. The test suite additionally needs
`pytest`, and building this documentation needs the toolchain listed in
`docs/requirements.txt`:

```sh
pip install pytest
pip install -r docs/requirements.txt
```

## Build targets

```sh
make            # build the simulation executable
make test       # build and run the Fortran unit tests
make pytest     # run the Python test suite
make docs       # build this documentation into docs/_build/html
make ui         # launch the GUI
make report     # compile the LaTeX report
make clean      # remove build artifacts and LaTeX leftovers
make help       # list the targets
```

## Platform note

`main.f` creates its output directory through a Windows `mkdir` call. On
Linux or macOS, create the directory before the first run:

```sh
mkdir -p data/bhw10_nstep1000000
```

Everything else, including the Fortran unit tests, is portable.
