"""
Pure-Python MCMC simulation for ground state position histogram.

Implements the same heat-bath + over-relaxation algorithm as the Fortran code,
using NumPy vectorization with checkerboard decomposition for performance.

Only collects y positions (no correlators/energy), as only the position
distribution is needed for the ground state histogram.
"""

import numpy as np


def run_ground_state_simulation(bhw, nt, nsteps, therm_steps,
                                progress_callback=None, cancel_flag=None):
    """
    Run MCMC simulation and return all y-positions for the histogram.

    Parameters
    ----------
    bhw : float
        beta * hbar * omega.
    nt : int
        Number of time slices.
    nsteps : int
        Number of MCMC configurations to collect (after thermalization).
    therm_steps : int
        Number of thermalization steps to discard.
    progress_callback : callable or None
        Called as progress_callback(current_step, total_steps) periodically.
    cancel_flag : callable or None
        If provided, checked periodically; if it returns True, simulation
        is aborted and None is returned.

    Returns
    -------
    all_positions : ndarray of shape (nsteps * nt,) or None if cancelled.
        All y-positions collected from all configurations.
    """
    eta = float(bhw) / nt
    alpha = eta / 2.0 + 1.0 / eta
    sigma = 1.0 / np.sqrt(2.0 * alpha)
    mu_coeff = 1.0 / (2.0 * alpha * eta)

    # Initialize path (cold start)
    y = np.zeros(nt, dtype=np.float64)

    # Precompute neighbour index arrays for checkerboard decomposition
    even = np.arange(0, nt, 2)
    odd = np.arange(1, nt, 2)
    left_even = (even - 1) % nt
    right_even = (even + 1) % nt
    left_odd = (odd - 1) % nt
    right_odd = (odd + 1) % nt

    total_steps = therm_steps + nsteps
    all_positions = np.empty(nsteps * nt, dtype=np.float64)
    progress_interval = max(1, total_steps // 100)

    for step in range(total_steps):
        if cancel_flag is not None and cancel_flag():
            return None

        # --- total_update: 10 × (1 HB + 5 OR) ---
        for _ in range(10):
            # Heat-bath sweep (checkerboard)
            mu_e = (y[left_even] + y[right_even]) * mu_coeff
            y[even] = np.random.normal(mu_e, sigma)
            mu_o = (y[left_odd] + y[right_odd]) * mu_coeff
            y[odd] = np.random.normal(mu_o, sigma)

            # Over-relaxation sweeps (checkerboard)
            for _ in range(5):
                mu_e = (y[left_even] + y[right_even]) * mu_coeff
                y[even] = 2.0 * mu_e - y[even]
                mu_o = (y[left_odd] + y[right_odd]) * mu_coeff
                y[odd] = 2.0 * mu_o - y[odd]

        # Collect data after thermalization
        if step >= therm_steps:
            idx = step - therm_steps
            all_positions[idx * nt:(idx + 1) * nt] = y

        if progress_callback and step % progress_interval == 0:
            progress_callback(step, total_steps)

    return all_positions
