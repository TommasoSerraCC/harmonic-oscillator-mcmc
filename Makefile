# ===========================================================
#  Quantum harmonic oscillator - path integral MCMC
# ===========================================================
#  make            build the simulation executable (main)
#  make test       build and run the Fortran unit tests
#  make ui         launch the Tkinter analysis GUI
#  make report     compile relazione/relazione.tex (twice)
#  make clean      remove build artifacts and LaTeX leftovers
# ===========================================================

FC     := gfortran
FFLAGS := -O2
PY     := python

ifeq ($(OS),Windows_NT)
  EXE := .exe
else
  EXE :=
endif

MCMC_SRC := mcmc/mcmc_oscillator.f mcmc/phys_sub.f mcmc/rand.f
MAIN_SRC := main.f $(MCMC_SRC)
TEST_SRC := test/test.f $(MCMC_SRC)

CLEANFILES := main$(EXE) test$(EXE) texput.log \
              relazione/*.aux relazione/*.log \
              relazione/*.toc relazione/*.out

.PHONY: all test ui report clean help

all: main$(EXE)

main$(EXE): $(MAIN_SRC)
	$(FC) $(FFLAGS) -o $@ $(MAIN_SRC)

test$(EXE): $(TEST_SRC)
	$(FC) $(FFLAGS) -o $@ $(TEST_SRC)

test: test$(EXE)
	./test$(EXE)

ui:
	$(PY) run_ui.py

report:
	cd relazione && pdflatex -interaction=nonstopmode relazione.tex
	cd relazione && pdflatex -interaction=nonstopmode relazione.tex

clean:
	@$(PY) -c "import glob,os;[os.remove(f) for p in '$(CLEANFILES)'.split() for f in glob.glob(p)]"
	@echo Cleaned.

help:
	@echo "make          - build the simulation executable"
	@echo "make test     - build and run the Fortran unit tests"
	@echo "make ui       - launch the analysis GUI"
	@echo "make report   - compile the LaTeX report"
	@echo "make clean    - remove build artifacts"
