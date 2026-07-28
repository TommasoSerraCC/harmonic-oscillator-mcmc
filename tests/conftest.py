"""Shared fixtures for the test suite.

The tests never touch the real ``data/`` directory: every input is built
on the fly in a temporary directory, either from analytically known
series or from short runs of the Python re-implementation of the chain.
"""

import numpy as np
import pytest

from helpers import make_raw_array


@pytest.fixture
def rng():
    """A seeded random generator, so failures are reproducible."""
    return np.random.default_rng(20260211)


@pytest.fixture
def raw_dataset(tmp_path, rng):
    """Write a synthetic raw chain and return its description.

    Returns
    -------
    dict
        Keys ``path``, ``array``, ``nsteps``, ``ncorr``, ``nt``, ``eta``,
        ``xi`` and ``expected_gap``.
    """
    nsteps, ncorr, xi = 2000, 8, 3.0
    nt = 2 * ncorr
    eta = 0.5
    array = make_raw_array(rng, nsteps, ncorr, xi)

    path = tmp_path / f'raw_data_nt{nt}.dat'
    np.savetxt(path, array, fmt='%20.12e')

    return {
        'path': str(path),
        'array': array,
        'nsteps': nsteps,
        'ncorr': ncorr,
        'nt': nt,
        'eta': eta,
        'xi': xi,
        'expected_gap': 1.0 / (xi * eta),
    }
