# Fortran core

The simulation is fixed-form Fortran 77, built by `make` from four
files.

## `main.f`

The driver. It reads its parameters from a namelist on standard input,
creates the output directory, and loops over the requested lattice
sizes. For each $N_t$ it starts from a cold path, applies
`total_update` `nsteps` times, and writes one line of measurements per
step.

`cold_start(y, nt)`
: Initialises every site of the path to zero.

## `mcmc/mcmc_oscillator.f`

The algorithm.

`heat_bath_sweep(y, nt, sigma, alpha, eta)`
: One sweep of the lattice, replacing each site with a draw from its
  Gaussian conditional distribution. Periodic neighbours are handled
  explicitly at both ends.

`microcanonical_sweep(y, nt, alpha, eta)`
: One over-relaxation sweep, reflecting each site about the conditional
  mean. Leaves the Euclidean action invariant.

`total_update(y, nt, sigma, alpha, eta)`
: The production update: ten repetitions of one heat-bath sweep followed
  by five over-relaxation sweeps. One call corresponds to one recorded
  measurement.

`box_muller(x, mu, sigma)`
: Returns a normal deviate of mean `mu` and standard deviation `sigma`,
  built from two uniform variates.

`get_indexes(idx, nt, il, ir)`
: Left and right neighbour indices with periodic wrap-around.

`euclidean_action(s, y, nt, eta, alpha)`
: Evaluates the discrete action of a path. Used by the test suite to
  verify that over-relaxation conserves it.

## `mcmc/phys_sub.f`

Observables and correlators.

`y1`, `y2`, `y3`, `A`
: Local observables at a site: $y$, $y^2$, $y^3$ and
  $A = y^3 - \tfrac32 y$.

`path_observable(y, nt, obs_func, result)`
: Averages any of the above over the sites of a path. The observable is
  passed as an external function.

`y1_corr`, `y2_corr`, `y3_corr`, `A_corr`
: Correlator wrappers evaluating $O(i)\,O(i+n)$ with wrap-around. The
  separation `n` is set beforehand through `set_corr_param` and shared
  via the `/corr_params/` common block.

`path_ene(y, nt, eta, energy)`
: The lattice energy estimator.

## `mcmc/rand.f`

`ran2()`
: The `ran2` combined generator from *Numerical Recipes*, returning a
  uniform variate in $[0,1)$.

`ranstart` / `ranfinish`
: Read and write the generator state to a file named `randomseed`, so
  consecutive runs continue the same stream. If the file is missing or
  incomplete, a default seed is used.

## `test/test.f`

The unit tests; see {doc}`testing`.

## A note on fixed-form line length

The sources are fixed-form Fortran, where anything past column 72 is
ignored by the compiler. Continuation lines carry a marker in column 6.
Keep statements within the limit when editing: exceeding it truncates
code silently rather than raising an error.
