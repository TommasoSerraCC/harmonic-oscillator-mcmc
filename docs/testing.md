# Testing

There are two suites, one per language, and both run in continuous
integration on every push.

```sh
make test      # Fortran unit tests
make pytest    # Python test suite
```

## Fortran suite

Built from `test/test.f` and run as `./run_tests`. It exits with a
non-zero status if any check fails, so it can gate a build.

| Test | Property checked |
|---|---|
| `test_ran2` | the uniform generator has mean $1/2$ and variance $1/12$ |
| `test_box_muller` | the normal deviates have mean 0, variance 1 and kurtosis 3 |
| `test_get_indexes` | neighbour indices wrap around correctly |
| `test_action_invariance` | **the over-relaxation sweep leaves the Euclidean action invariant**, to within round-off |
| `test_heat_bath` | heat-bath sampling alone reaches equilibrium: $\langle y\rangle = 0$, $\langle y^2\rangle \approx 1/2$ |
| `test_total_update` | the production update reaches the same equilibrium |
| `test_path_ene` | the energy estimator returns $\approx 1/2$ at low temperature |
| `test_correlator` | at zero separation the correlator machinery reduces exactly to $y^2$ |

The action-invariance test is the sharpest of these: it is not a
statistical statement but an exact algebraic property, and a relative
drift above $10^{-8}$ means the update is wrong.

## Python suite

Run with `pytest` from the repository root. No test touches the real
`data/` directory; inputs are generated into temporary directories.

The suite is built around checking estimators against answers that are
known independently, rather than against previous output:

**Analytically known series.** An AR(1) process has autocorrelation
$\rho^k$ and integrated time $\rho/(1-\rho)$, so the blocking plateau
must sit at $\sqrt{1+2\tau_{\text{int}}}$ times the naive standard
error. The autocorrelation estimator is checked against $\rho^k$
directly.

**Brute-force cross-checks.** The vectorised jackknife computes
leave-one-block-out means by subtracting block sums from the total,
which is fast but easy to get subtly wrong. The tests compare it against
an explicit loop that rebuilds each reduced sample, for both the
connected correlator and the energy gap.

**Exact special cases.** A constant correlator with a zero-mean
observable must give a zero uncertainty; a pure exponential correlator
must give a flat gap equal to $1/(\xi\eta)$; blocking with unit blocks
must reduce exactly to the standard error of the mean.

**Physics from the sampler.** The NumPy chain is run at low temperature
and small lattice spacing, and the sampled positions are required to
reproduce the ground state: mean 0, variance $1/2$, vanishing odd
moments and Gaussian kurtosis. A separate test verifies that a coarse
lattice undershoots the continuum variance, matching the sign of the
$\eta^2$ correction fitted in the report.

**Round-trip through the loaders.** A real analysis is run into a
temporary tree and read back through every function in
{mod}`ui.core.data_manager`. This pins the contract between the file
formats written by the pipeline and the parsers used by the GUI, so a
change to either side that breaks the other fails the suite.

**End-to-end pipeline.** A synthetic chain with a built-in correlation
length is analysed, and the recovered connected correlators and energy
gaps are compared with the values that went in.

## Continuous integration

Two workflows run on GitHub Actions:

`ci.yml`
: Builds the Fortran on Ubuntu and runs its unit tests, then runs the
  Python suite against Python 3.9, 3.11 and 3.12.

`docs.yml`
: Builds this documentation and, on pushes to `main`, deploys it to
  GitHub Pages. Warnings are treated as errors, so a broken
  cross-reference fails the build instead of shipping.
