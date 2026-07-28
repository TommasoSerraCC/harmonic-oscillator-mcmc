"""Tests for the NumPy re-implementation of the Markov chain.

:mod:`ui.core.ground_state_sim` mirrors the Fortran algorithm with a
checkerboard decomposition. These tests check that it reproduces the
physics it is supposed to: at low temperature and small lattice spacing
the sampled positions must follow the ground-state distribution
``|psi_0(y)|^2 = exp(-y^2)/sqrt(pi)``, which has mean 0 and variance 1/2.
"""

import numpy as np
import pytest

from ui.core.ground_state_sim import run_ground_state_simulation

# beta*hbar*omega = 10 freezes the system into the ground state,
# eta = 0.1 keeps the discretisation error small.
BHW = 10.0
NT = 100
NSTEPS = 400
THERM = 100


@pytest.fixture(scope='module')
def positions():
    """One short simulation reused by the physics assertions."""
    np.random.seed(20260211)
    return run_ground_state_simulation(BHW, NT, NSTEPS, THERM)


def test_returns_one_sample_per_site_and_configuration(positions):
    assert positions.shape == (NSTEPS * NT,)
    assert np.all(np.isfinite(positions))


def test_distribution_is_symmetric_about_the_origin(positions):
    assert positions.mean() == pytest.approx(0.0, abs=0.05)


def test_variance_matches_the_ground_state(positions):
    """<y^2> = 1/2 for the continuum ground state."""
    assert positions.var() == pytest.approx(0.5, abs=0.05)


def test_odd_moments_vanish(positions):
    """The potential is even, so <y> and <y^3> must vanish."""
    assert np.mean(positions ** 3) == pytest.approx(0.0, abs=0.05)


def test_kurtosis_matches_a_gaussian(positions):
    """|psi_0|^2 is Gaussian, so <y^4>/<y^2>^2 = 3."""
    ratio = np.mean(positions ** 4) / np.mean(positions ** 2) ** 2
    assert ratio == pytest.approx(3.0, rel=0.05)


def test_simulation_is_reproducible_under_a_fixed_seed():
    np.random.seed(7)
    first = run_ground_state_simulation(5.0, 20, 30, 10)
    np.random.seed(7)
    second = run_ground_state_simulation(5.0, 20, 30, 10)
    assert first == pytest.approx(second)


def test_cancelling_aborts_and_returns_nothing():
    result = run_ground_state_simulation(5.0, 20, 100, 10,
                                         cancel_flag=lambda: True)
    assert result is None


def test_progress_callback_is_driven_to_completion():
    seen = []
    run_ground_state_simulation(5.0, 20, 100, 20,
                                progress_callback=lambda s, t: seen.append((s, t)))

    assert seen, 'the progress callback was never called'
    steps, totals = zip(*seen)
    assert set(totals) == {120}
    assert list(steps) == sorted(steps)
    assert steps[0] == 0


def test_a_coarse_lattice_undershoots_the_continuum_variance():
    """Finite-step bias, approaching 1/2 from below as eta decreases.

    Both runs are at the same temperature; only the lattice spacing
    differs, eta = 0.1 against eta = 2. The continuum limit is reached
    from below, which is the same sign as the negative eta^2 coefficient
    fitted for the energy in the report.
    """
    np.random.seed(11)
    fine = run_ground_state_simulation(20.0, 200, 200, 50)
    np.random.seed(11)
    coarse = run_ground_state_simulation(20.0, 10, 200, 50)

    assert fine.var() == pytest.approx(0.5, abs=0.05)
    assert coarse.var() < fine.var() - 0.05
