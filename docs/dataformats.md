# Data formats

Three directories carry data, with a naming convention that encodes the
run parameters so that several simulations can coexist:

```
data/bhw{N}_nstep{M}/            raw chains          (not tracked by git)
results/bhw{N}_nstep{M}/nt{T}_therm{S}/    analysis output  (tracked)
plots/bhw{N}_nstep{M}/nt{T}_therm{S}/      figures          (tracked)
```

`N` is $\beta\hbar\omega$, `M` the number of recorded measurements, `T`
the lattice size $N_t$ and `S` the number of thermalisation
measurements discarded during the analysis.

## Raw chains

### `raw_data_nt{T}.dat`

One line per recorded measurement, with $5 + 4\,(N_t/2)$ columns:

| Columns | Contents |
|---|---|
| 1–5 | $y$, $y^2$, $y^3$, $A$, $E$ |
| $6+4k$ … $9+4k$ | $y(0)y(n)$, $y^2(0)y^2(n)$, $y^3(0)y^3(n)$, $A(0)A(n)$ for $n = k+1$ |

Every entry is already averaged over the $N_t$ sites of the path. The
correlators are stored **raw**, not connected: the subtraction of
$\langle O\rangle^2$ happens during the analysis, inside the jackknife,
so that the correlation between the two terms is handled correctly.

These files are large — gigabytes for the longer runs — and are excluded
from version control.

### `raw_energy_nt{T}.dat`

Produced when `energy_only=1`. A single column holding the energy
estimator, nothing else.

## Analysis output

Written by {func}`analysis.pipeline.analyze_dataset`.

| File | Contents |
|---|---|
| `observables.dat` | one line per observable: name, mean, blocking error |
| `blocking_observables.dat` | `k` followed by $\sigma$ for each of the five observables, one row per block size $k = 2^0 \dots 2^{k_{\max}}$ |
| `connected_correlators.dat` | `n` followed by mean and jackknife error for each of the four correlators |
| `blocking_jackknife_correlators.dat` | `k` followed by the jackknife error of each correlator at three representative separations; the separations are recorded in the header |
| `energy_gaps.dat` | `n` followed by $\Delta E(n)$ and its error for each correlator |
| `tau_int.dat` | integrated autocorrelation time per observable |
| `tau_exp_fit.dat` | lag and $C_{y^2}$, with the fitted `A` and `tau` in the header comment |
| `correlation_length.dat` | $\xi$ and $\Delta E$ from the exponential fits, written by `fit_correlation_length.py` |
| `energy_only.marker` | present only for energy-only datasets |
| `fits/` | archived interactive fits, one JSON and one PNG per entry |

All tabular files carry a `#` header line naming the columns, and are
loaded by the functions in {mod}`ui.core.data_manager`. The test suite
runs a real analysis and reads it back through those loaders, so the two
sides cannot drift apart unnoticed.

## Error estimates

Two different procedures are used, because the quantities they apply to
behave differently.

**Blocking**, for the primary observables. The chain is cut into blocks
of increasing size, and the error is the standard error of the block
means. Once the block length exceeds the autocorrelation time the block
means become independent and the estimate stops growing; that plateau is
the reported error. The full scan is saved so the plateau can be
inspected rather than assumed.

**Jackknife over blocks**, for the connected correlators and the energy
gaps. These are non-linear functions of ensemble averages, so naive
propagation is invalid. Each block is removed in turn, the quantity is
recomputed on the remainder, and the spread of those estimates gives the
error. For the gaps the jackknife is applied directly to
$\log[C(n)/C(n+1)]$ rather than to the two correlators separately,
because they are strongly correlated.
