"""Tests for the estimators in :mod:`analysis.statistics`.

Wherever possible the vectorised implementations are checked against an
independent brute-force version written inside the test, or against a
series whose exact statistical properties are known analytically.
"""

import numpy as np
import pytest

from analysis import statistics as st
from helpers import ar1_series


# --------------------------------------------------------------------
#  Blocking
# --------------------------------------------------------------------

def test_blocking_with_unit_blocks_is_the_naive_standard_error(rng):
    """With one measurement per block, blocking must reduce to the SEM."""
    x = rng.normal(size=500)
    expected = x.std(ddof=1) / np.sqrt(x.size)
    assert st.blocking_sigma(x, 1) == pytest.approx(expected)


def test_blocking_ignores_a_trailing_partial_block(rng):
    """A partial trailing block is dropped, not padded."""
    x = rng.normal(size=101)
    assert st.blocking_sigma(x, 10) == pytest.approx(
        st.blocking_sigma(x[:100], 10))


def test_blocking_grows_towards_the_correlated_plateau(rng):
    """For AR(1) the plateau sits at sqrt(1 + 2 tau_int) times the SEM."""
    rho = 0.8
    x = ar1_series(rng, 200_000, rho)

    naive = st.blocking_sigma(x, 1)
    plateau = st.blocking_sigma(x, 2048)

    tau_int = rho / (1.0 - rho)
    expected_ratio = np.sqrt(1.0 + 2.0 * tau_int)
    assert plateau / naive == pytest.approx(expected_ratio, rel=0.1)


def test_blocking_is_flat_for_independent_data(rng):
    """Uncorrelated data show no plateau growth."""
    x = rng.normal(size=100_000)
    assert st.blocking_sigma(x, 512) / st.blocking_sigma(x, 1) == pytest.approx(
        1.0, rel=0.15)


def test_blocking_rejects_too_few_blocks(rng):
    x = rng.normal(size=10)
    with pytest.raises(ValueError, match='at least 2 blocks'):
        st.blocking_sigma(x, 8)


def test_blocking_rejects_non_positive_block_size(rng):
    with pytest.raises(ValueError, match='must be positive'):
        st.blocking_sigma(rng.normal(size=10), 0)


def test_blocking_scan_shape_and_sizes(rng):
    x = rng.normal(size=4096)
    k_values, sigmas = st.blocking_scan(x, kmax=5)
    assert list(k_values) == [1, 2, 4, 8, 16, 32]
    assert sigmas.shape == (6,)
    assert np.all(np.isfinite(sigmas))


def test_max_block_exponent_leaves_the_requested_blocks():
    assert 2 ** st.max_block_exponent(4000, min_blocks=4) <= 1000
    assert st.max_block_exponent(1024, min_blocks=4) == 8


# --------------------------------------------------------------------
#  Jackknife
# --------------------------------------------------------------------

def _brute_force_connected(corr, obs, block_size):
    """Straightforward leave-one-block-out reference implementation."""
    n = len(obs)
    nblocks = n // block_size
    used = nblocks * block_size
    samples = []
    for i in range(nblocks):
        keep = np.ones(used, dtype=bool)
        keep[i * block_size:(i + 1) * block_size] = False
        samples.append(corr[:used][keep].mean() - obs[:used][keep].mean() ** 2)
    samples = np.array(samples)
    mean = samples.mean()
    sigma = np.sqrt((nblocks - 1) / nblocks * np.sum((samples - mean) ** 2))
    return mean, sigma


def test_leave_one_block_out_average_recovers_the_global_mean(rng):
    """The leave-one-out means must average back to the overall mean."""
    x = rng.normal(size=600)
    reduced = st._leave_one_block_out(x, 50)
    assert reduced.mean() == pytest.approx(x[:600].mean())


def test_jackknife_connected_matches_brute_force(rng):
    """The vectorised jackknife must equal the explicit loop."""
    obs = rng.normal(0.3, 1.0, size=1000)
    corr = obs * np.roll(obs, 3) + rng.normal(0.0, 0.1, size=1000)

    mean, sigma = st.jackknife_connected(corr, obs, 50)
    ref_mean, ref_sigma = _brute_force_connected(corr, obs, 50)

    assert mean == pytest.approx(ref_mean)
    assert sigma == pytest.approx(ref_sigma)


def test_jackknife_all_n_matches_the_single_separation_version(rng):
    """Handling every separation at once must not change the answer."""
    nsteps, ncorr = 1200, 5
    obs = rng.normal(0.2, 1.0, size=nsteps)
    corr = rng.normal(1.0, 0.2, size=(nsteps, ncorr))

    means, sigmas = st.jackknife_connected_all_n(corr, obs, 60)
    for n in range(ncorr):
        mean, sigma = st.jackknife_connected(corr[:, n], obs, 60)
        assert means[n] == pytest.approx(mean)
        assert sigmas[n] == pytest.approx(sigma)


def test_jackknife_connected_recovers_a_known_value(rng):
    """A constant correlator with a zero-mean observable is exact."""
    obs = np.zeros(1000)
    corr = np.full(1000, 0.75)
    mean, sigma = st.jackknife_connected(corr, obs, 100)
    assert mean == pytest.approx(0.75)
    assert sigma == pytest.approx(0.0, abs=1e-12)


def test_jackknife_energy_gap_recovers_an_exponential_decay():
    """A pure exponential correlator gives a flat, exact gap."""
    nsteps, ncorr, xi, eta = 400, 6, 2.5, 0.4
    separations = np.arange(1, ncorr + 1)
    corr = np.tile(np.exp(-separations / xi), (nsteps, 1))
    obs = np.zeros(nsteps)

    gaps, sigmas = st.jackknife_energy_gap(corr, obs, 50, eta)

    assert gaps == pytest.approx(np.full(ncorr - 1, 1.0 / (xi * eta)))
    assert sigmas == pytest.approx(np.zeros(ncorr - 1), abs=1e-10)


def test_jackknife_energy_gap_matches_brute_force(rng):
    """Cross-check the ratio jackknife against an explicit loop."""
    nsteps, ncorr, block, eta = 800, 4, 40, 0.5
    obs = rng.normal(0.0, 0.05, size=nsteps)
    decay = np.exp(-np.arange(1, ncorr + 1) / 2.0)
    corr = decay[None, :] * (1.0 + rng.normal(0.0, 0.01, (nsteps, ncorr)))

    gaps, _ = st.jackknife_energy_gap(corr, obs, block, eta)

    nblocks = nsteps // block
    samples = []
    for i in range(nblocks):
        keep = np.ones(nsteps, dtype=bool)
        keep[i * block:(i + 1) * block] = False
        conn = corr[keep].mean(axis=0) - obs[keep].mean() ** 2
        samples.append(np.log(conn[:-1] / conn[1:]) / eta)
    assert gaps == pytest.approx(np.mean(samples, axis=0))


def test_jackknife_energy_gap_survives_a_negative_correlator():
    """Noise-dominated tails must not poison the whole result."""
    nsteps, eta = 400, 0.5
    corr = np.tile(np.array([1.0, 0.5, 0.25, -0.01]), (nsteps, 1))
    obs = np.zeros(nsteps)

    gaps, _ = st.jackknife_energy_gap(corr, obs, 50, eta)

    assert np.isfinite(gaps[:2]).all()
    assert np.isnan(gaps[2])


# --------------------------------------------------------------------
#  Autocorrelation
# --------------------------------------------------------------------

def test_autocorrelation_starts_at_one(rng):
    acf = st.autocorrelation(rng.normal(size=2000), 20)
    assert acf[0] == pytest.approx(1.0)


def test_autocorrelation_follows_the_ar1_law(rng):
    """For AR(1) the autocorrelation must decay as rho**k."""
    rho = 0.7
    acf = st.autocorrelation(ar1_series(rng, 200_000, rho), 8)
    expected = rho ** np.arange(8)
    assert acf == pytest.approx(expected, abs=0.03)


def test_integrated_time_matches_the_ar1_prediction(rng):
    rho = 0.7
    acf = st.autocorrelation(ar1_series(rng, 200_000, rho), 60)
    assert st.integrated_time(acf) == pytest.approx(rho / (1.0 - rho), rel=0.15)


def test_autocorrelation_rejects_a_constant_series():
    with pytest.raises(ValueError, match='zero variance'):
        st.autocorrelation(np.ones(100), 5)


def test_autocorrelation_rejects_an_excessive_lag(rng):
    with pytest.raises(ValueError, match='exceeds'):
        st.autocorrelation(rng.normal(size=10), 50)


def test_exponential_fit_recovers_the_decay_constant():
    tau = 7.0
    acf = np.exp(-np.arange(60) / tau)
    amplitude, fitted = st.fit_exponential_autocorrelation(acf)
    assert fitted == pytest.approx(tau, rel=0.02)
    assert amplitude == pytest.approx(1.0, rel=0.05)


def test_exponential_fit_falls_back_when_it_cannot_converge():
    """A degenerate autocorrelation returns the supplied fallback."""
    acf = np.array([1.0, 0.0])
    _, tau = st.fit_exponential_autocorrelation(acf, fallback_tau=3.5)
    assert tau == pytest.approx(3.5)


# --------------------------------------------------------------------
#  Helpers
# --------------------------------------------------------------------

@pytest.mark.parametrize('ncorr, expected', [
    (100, [1, 25, 50]),
    (16, [1, 4, 8]),
    (4, [1, 2]),
])
def test_select_check_points(ncorr, expected):
    assert st.select_check_points(ncorr) == expected


def test_select_check_points_stays_inside_the_range():
    for ncorr in range(2, 40):
        points = st.select_check_points(ncorr)
        assert len(points) <= 3
        assert points == sorted(set(points))
        assert all(0 <= p < ncorr for p in points)
