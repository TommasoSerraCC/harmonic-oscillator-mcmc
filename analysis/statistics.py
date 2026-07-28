"""Estimators for correlated Markov chain data.

Successive MCMC measurements are not independent, so the naive standard
error of the mean underestimates the true uncertainty. Two standard
remedies are implemented here:

* **blocking** for primary observables — average over blocks of
  increasing size until the estimated error stops growing (the plateau);
* **jackknife over blocks** for quantities that are non-linear functions
  of ensemble averages, such as the connected correlator
  ``C(n) = <O(0)O(n)> - <O>^2`` and the effective gap
  ``log[C(n)/C(n+1)]/eta``, where naive error propagation fails because
  the terms are strongly correlated.

All functions take raw per-measurement arrays and are free of file I/O,
so they can be exercised directly by the test suite.
"""

import warnings

import numpy as np

__all__ = [
    "autocorrelation",
    "blocking_scan",
    "blocking_sigma",
    "fit_exponential_autocorrelation",
    "integrated_time",
    "jackknife_connected",
    "jackknife_connected_all_n",
    "jackknife_energy_gap",
    "max_block_exponent",
    "select_check_points",
]


# --------------------------------------------------------------------
#  Blocking
# --------------------------------------------------------------------

def blocking_sigma(x, block_size):
    """Standard error of the mean of ``x`` using blocks of ``block_size``.

    The series is cut into non-overlapping blocks, the mean of each block
    is taken, and the error is the standard deviation of the block means
    divided by the square root of their number. As ``block_size`` grows
    past the autocorrelation time the block means become independent and
    the estimate reaches a plateau.

    Parameters
    ----------
    x : array_like
        One-dimensional series of measurements.
    block_size : int
        Number of consecutive measurements per block. A trailing partial
        block is discarded.

    Returns
    -------
    float
        Estimated standard error of ``mean(x)``.

    Raises
    ------
    ValueError
        If fewer than two complete blocks fit in the series.
    """
    x = np.asarray(x, dtype=float).ravel()
    block_size = int(block_size)
    if block_size < 1:
        raise ValueError("block_size must be positive")
    nblocks = x.size // block_size
    if nblocks < 2:
        raise ValueError(
            f"need at least 2 blocks, got {nblocks} "
            f"({x.size} points, block_size={block_size})"
        )
    block_means = x[:nblocks * block_size].reshape(nblocks, block_size).mean(axis=1)
    return float(block_means.std(ddof=1) / np.sqrt(nblocks))


def max_block_exponent(nsteps, min_blocks=4):
    """Largest ``p`` such that blocks of size ``2**p`` leave ``min_blocks``.

    Parameters
    ----------
    nsteps : int
        Length of the series.
    min_blocks : int, optional
        Minimum number of blocks to retain. Default is 4.

    Returns
    -------
    int
        The exponent ``p``.
    """
    return int(np.log2(nsteps // min_blocks))


def blocking_scan(x, kmax=None):
    """Blocking error as a function of the block size, for plateau hunting.

    Parameters
    ----------
    x : array_like
        One-dimensional series of measurements.
    kmax : int, optional
        Largest exponent to probe; block sizes are ``2**0 … 2**kmax``.
        Defaults to :func:`max_block_exponent`.

    Returns
    -------
    k_values : ndarray of int
        The block sizes probed.
    sigmas : ndarray of float
        The corresponding estimates of the error on the mean.
    """
    x = np.asarray(x, dtype=float).ravel()
    if kmax is None:
        kmax = max_block_exponent(x.size)
    k_values = 2 ** np.arange(kmax + 1)
    sigmas = np.array([blocking_sigma(x, int(k)) for k in k_values])
    return k_values, sigmas


# --------------------------------------------------------------------
#  Jackknife
# --------------------------------------------------------------------

def _leave_one_block_out(values, block_size):
    """Leave-one-block-out means of ``values`` along the first axis.

    Instead of rebuilding each reduced sample explicitly, the per-block
    sums are computed once and subtracted from the total. This is what
    makes the jackknife affordable on chains of millions of steps.

    Parameters
    ----------
    values : ndarray
        Array of shape ``(n,)`` or ``(n, m)``.
    block_size : int
        Block length.

    Returns
    -------
    ndarray
        Array of shape ``(nblocks,)`` or ``(nblocks, m)`` holding the mean
        of the series with one block removed.
    """
    n = values.shape[0]
    nblocks = n // block_size
    if nblocks < 2:
        raise ValueError(
            f"need at least 2 blocks, got {nblocks} "
            f"({n} points, block_size={block_size})"
        )
    used = nblocks * block_size
    trimmed = values[:used]
    shape = (nblocks, block_size) + values.shape[1:]
    block_sums = trimmed.reshape(shape).sum(axis=1)
    return (block_sums.sum(axis=0) - block_sums) / (used - block_size)


def _jackknife_error(samples, mean):
    """Jackknife uncertainty from leave-one-out estimates."""
    nblocks = samples.shape[0]
    scale = (nblocks - 1) / nblocks
    return np.sqrt(scale * np.nansum((samples - mean) ** 2, axis=0))


def jackknife_connected(corr, obs, block_size):
    """Connected correlator and its jackknife error at one separation.

    Parameters
    ----------
    corr : array_like
        Series of the raw product ``O(i) O(i+n)``, one entry per
        measurement.
    obs : array_like
        Series of the observable ``O`` itself, same length.
    block_size : int
        Jackknife block length.

    Returns
    -------
    mean : float
        Jackknife estimate of ``<O(0)O(n)> - <O>^2``.
    sigma : float
        Its uncertainty.
    """
    corr = np.asarray(corr, dtype=float).ravel()
    obs = np.asarray(obs, dtype=float).ravel()
    corr_jack = _leave_one_block_out(corr, block_size)
    obs_jack = _leave_one_block_out(obs, block_size)
    conn = corr_jack - obs_jack ** 2
    mean = float(conn.mean())
    return mean, float(_jackknife_error(conn, mean))


def jackknife_connected_all_n(corr_matrix, obs, block_size):
    """Vectorised :func:`jackknife_connected` over every separation at once.

    Parameters
    ----------
    corr_matrix : array_like
        Array of shape ``(nsteps, ncorr)`` with the raw products for each
        separation ``n``.
    obs : array_like
        Series of the observable, shape ``(nsteps,)``.
    block_size : int
        Jackknife block length.

    Returns
    -------
    means : ndarray
        Connected correlators, shape ``(ncorr,)``.
    sigmas : ndarray
        Their uncertainties, shape ``(ncorr,)``.
    """
    corr_matrix = np.asarray(corr_matrix, dtype=float)
    obs = np.asarray(obs, dtype=float).ravel()
    corr_jack = _leave_one_block_out(corr_matrix, block_size)
    obs_jack = _leave_one_block_out(obs, block_size)
    conn = corr_jack - obs_jack[:, None] ** 2
    means = conn.mean(axis=0)
    return means, _jackknife_error(conn, means[None, :])


def jackknife_energy_gap(corr_matrix, obs, block_size, eta):
    """Effective energy gaps ``log[C(n)/C(n+1)]/eta`` with jackknife errors.

    The jackknife is applied to the ratio itself rather than propagated
    from the two correlators, because ``C(n)`` and ``C(n+1)`` are strongly
    correlated and independent propagation would overestimate the error.

    Ratios that come out non-positive, which happens once the correlator
    is buried in noise, are discarded from the average instead of
    producing a NaN gap.

    Parameters
    ----------
    corr_matrix : array_like
        Raw products, shape ``(nsteps, ncorr)``.
    obs : array_like
        Series of the observable, shape ``(nsteps,)``.
    block_size : int
        Jackknife block length.
    eta : float
        Lattice spacing in units of ``1/omega``.

    Returns
    -------
    gaps : ndarray
        Effective gaps, shape ``(ncorr - 1,)``.
    sigmas : ndarray
        Their uncertainties, same shape.
    """
    corr_matrix = np.asarray(corr_matrix, dtype=float)
    obs = np.asarray(obs, dtype=float).ravel()
    corr_jack = _leave_one_block_out(corr_matrix, block_size)
    obs_jack = _leave_one_block_out(obs, block_size)
    conn = corr_jack - obs_jack[:, None] ** 2

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = conn[:, :-1] / conn[:, 1:]
        ratio = np.where(ratio > 0, ratio, np.nan)
        gap = np.log(ratio) / eta

    # A separation whose ratio is negative in every jackknife sample
    # yields an all-NaN column; that is an expected outcome deep in the
    # noise, so the resulting "empty slice" warning is suppressed.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mean = np.nanmean(gap, axis=0)
    return mean, _jackknife_error(gap, mean[None, :])


# --------------------------------------------------------------------
#  Autocorrelation
# --------------------------------------------------------------------

def autocorrelation(x, max_lag):
    """Normalised autocorrelation function of a series.

    The estimator divides by the number of terms actually summed at each
    lag, so it stays unbiased as the lag approaches the series length.

    Parameters
    ----------
    x : array_like
        One-dimensional series.
    max_lag : int
        Number of lags to return, starting at zero.

    Returns
    -------
    ndarray
        Autocorrelation for lags ``0 … max_lag - 1``; the first entry is 1.
    """
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    if max_lag > n:
        raise ValueError("max_lag exceeds the length of the series")
    centred = x - x.mean()
    var = float(np.dot(centred, centred))
    if var == 0.0:
        raise ValueError("series has zero variance")
    return np.array([
        np.dot(centred[:n - t], centred[t:]) / var * n / (n - t)
        for t in range(max_lag)
    ])


def integrated_time(acf):
    """Integrated autocorrelation time, ``sum_{k>=1} C(k)``.

    Parameters
    ----------
    acf : array_like
        Autocorrelation function as returned by :func:`autocorrelation`.

    Returns
    -------
    float
        The summed tail, excluding the trivial lag-zero term.
    """
    return float(np.sum(np.asarray(acf, dtype=float).ravel()[1:]))


def fit_exponential_autocorrelation(acf, fallback_tau=5.0):
    """Fit ``A exp(-t/tau)`` to the leading part of an autocorrelation function.

    The fit range stops where the autocorrelation has decayed below one
    per cent, with a floor of ten points so that a very fast decay still
    leaves something to fit.

    Parameters
    ----------
    acf : array_like
        Autocorrelation function, lag zero first.
    fallback_tau : float, optional
        Value returned for ``tau`` if the fit does not converge.

    Returns
    -------
    amplitude : float
        Fitted ``A``.
    tau : float
        Fitted exponential autocorrelation time.
    """
    from scipy.optimize import curve_fit

    acf = np.asarray(acf, dtype=float).ravel()

    def exp_decay(t, amplitude, tau):
        return amplitude * np.exp(-t / tau)

    idx_zero = max(int(np.argmax(acf[1:] < 0.01)) + 1, 10)
    idx_zero = min(idx_zero, acf.size)
    t_fit = np.arange(1, idx_zero)
    if t_fit.size < 2:
        return 1.0, fallback_tau
    try:
        popt, _ = curve_fit(
            exp_decay, t_fit, acf[1:idx_zero], p0=[1.0, 5.0],
            bounds=([0.0, 0.1], [np.inf, 100.0]),
        )
    except Exception:
        return 1.0, fallback_tau
    return float(popt[0]), float(popt[1])


# --------------------------------------------------------------------
#  Misc
# --------------------------------------------------------------------

def select_check_points(ncorr):
    """Pick up to three representative separations for the blocking scan.

    Running the full blocking scan for every separation would be wasteful,
    so the plateau is inspected at a short, an intermediate and a long
    separation.

    Parameters
    ----------
    ncorr : int
        Number of available separations.

    Returns
    -------
    list of int
        Sorted, de-duplicated indices, at most three, all below ``ncorr``.
    """
    candidates = [1, ncorr // 4, ncorr // 2]
    points = sorted({min(v, ncorr - 1) for v in candidates})
    if len(points) < 3 and ncorr > 3:
        candidates = [1, ncorr // 3, 2 * ncorr // 3]
        points = sorted({min(v, ncorr - 1) for v in candidates})
    return points[:3]
