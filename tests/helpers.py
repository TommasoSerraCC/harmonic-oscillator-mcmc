"""Synthetic data generators shared by the tests.

Both generators produce series whose statistical properties are known in
closed form, so the estimators can be checked against exact answers
rather than against a previous run of themselves.
"""

import numpy as np


def ar1_series(rng, n, rho, scale=1.0):
    """Generate an AR(1) series ``x[t] = rho x[t-1] + noise``.

    The stationary process has autocorrelation ``C(k) = rho**k`` and
    integrated autocorrelation time ``rho / (1 - rho)``, which makes it a
    convenient reference for the blocking and autocorrelation estimators.

    Parameters
    ----------
    rng : numpy.random.Generator
        Source of randomness.
    n : int
        Length of the series.
    rho : float
        Lag-one autocorrelation, ``0 <= rho < 1``.
    scale : float, optional
        Standard deviation of the stationary distribution.

    Returns
    -------
    ndarray
        The generated series.
    """
    noise = rng.normal(0.0, scale * np.sqrt(1.0 - rho ** 2), size=n)
    x = np.empty(n)
    x[0] = rng.normal(0.0, scale)
    for t in range(1, n):
        x[t] = rho * x[t - 1] + noise[t]
    return x


def make_raw_array(rng, nsteps, ncorr, xi, noise=0.02):
    """Build a synthetic raw chain with a known correlation length.

    The four observables that feed the correlators are centred on zero,
    so the connected correlator reduces to the raw product and equals
    ``exp(-n / xi)`` up to noise. The effective gap is therefore
    ``1 / (xi * eta)`` at every separation, which the tests can check.

    Parameters
    ----------
    rng : numpy.random.Generator
        Source of randomness.
    nsteps : int
        Number of measurements.
    ncorr : int
        Number of separations, i.e. ``N_t / 2``.
    xi : float
        Correlation length in lattice units.
    noise : float, optional
        Relative noise added to the correlators.

    Returns
    -------
    ndarray
        Array of shape ``(nsteps, 5 + 4 * ncorr)``, matching the column
        layout written by the Fortran simulation.
    """
    data = np.empty((nsteps, 5 + 4 * ncorr))

    # y, y2, y3, A centred on zero; E offset so it is distinguishable.
    for column in range(4):
        data[:, column] = rng.normal(0.0, 0.1, size=nsteps)
    data[:, 4] = 0.5 + rng.normal(0.0, 0.01, size=nsteps)

    separations = np.arange(1, ncorr + 1)
    decay = np.exp(-separations / xi)
    for j in range(4):
        jitter = rng.normal(0.0, noise, size=(nsteps, ncorr))
        data[:, 5 + j::4] = decay[None, :] * (1.0 + jitter)

    return data
