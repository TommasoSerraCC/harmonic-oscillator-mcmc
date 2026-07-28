"""End-to-end analysis of one raw Markov chain.

:func:`analyze_dataset` turns a ``raw_data_nt*.dat`` file produced by the
Fortran simulation into the table of results consumed by the GUI and by
the plotting scripts. The file names and column layouts written here are
part of the contract with :mod:`ui.core.data_manager`; see the *Data
formats* page of the documentation.
"""

import os

import numpy as np

from analysis import statistics as st

__all__ = [
    "CORRELATOR_NAMES",
    "OBSERVABLE_NAMES",
    "analyze_dataset",
    "analyze_energy_only",
    "split_raw_columns",
]

#: Primary observables, in the order written by the Fortran code.
OBSERVABLE_NAMES = ["y", "y2", "y3", "A", "E"]

#: Connected correlators, in the order written by the Fortran code.
CORRELATOR_NAMES = ["yc", "y2c", "y3c", "Ac"]

#: Observable matching each correlator, used to subtract ``<O>^2``.
_CORRELATOR_SOURCE = ["y", "y2", "y3", "A"]

#: Default jackknife and blocking block length.
DEFAULT_BLOCK = 1000

#: Default number of lags kept for the autocorrelation function.
DEFAULT_MAX_LAG = 50


def split_raw_columns(data):
    """Split a raw chain into primary observables and correlators.

    Parameters
    ----------
    data : ndarray
        Array of shape ``(nsteps, 5 + 4 * ncorr)`` as loaded from a
        ``raw_data_nt*.dat`` file.

    Returns
    -------
    observables : dict of str to ndarray
        Series for ``y``, ``y2``, ``y3``, ``A`` and ``E``.
    correlators : dict of str to ndarray
        Series of shape ``(nsteps, ncorr)`` for each correlator.

    Raises
    ------
    ValueError
        If the column count is not consistent with the expected layout.
    """
    data = np.atleast_2d(np.asarray(data, dtype=float))
    ncols = data.shape[1]
    if ncols < 5 or (ncols - 5) % 4 != 0:
        raise ValueError(
            f"unexpected raw layout: {ncols} columns, expected 5 + 4*ncorr"
        )
    observables = {name: data[:, i] for i, name in enumerate(OBSERVABLE_NAMES)}
    tail = data[:, 5:]
    correlators = {
        name: tail[:, i::4] for i, name in enumerate(CORRELATOR_NAMES)
    }
    return observables, correlators


def analyze_dataset(raw_path, outdir, eta, skip=0,
                    block_size=DEFAULT_BLOCK, max_lag=DEFAULT_MAX_LAG,
                    verbose=False):
    """Analyse one raw chain and write every result table.

    Parameters
    ----------
    raw_path : str
        Path to the ``raw_data_nt*.dat`` file.
    outdir : str
        Directory to write the results into; created if missing.
    eta : float
        Lattice spacing ``beta * hbar * omega / N_t``.
    skip : int, optional
        Leading measurements discarded as thermalisation.
    block_size : int, optional
        Block length for the reported blocking and jackknife errors.
    max_lag : int, optional
        Lags kept for the autocorrelation function.
    verbose : bool, optional
        Print a short progress report.

    Returns
    -------
    dict
        Summary with keys ``nsteps``, ``ncorr``, ``means``, ``errors``,
        ``tau_int``, ``tau_exp`` and ``files``.

    Raises
    ------
    ValueError
        If, after ``skip``, too few measurements remain to form blocks.
    """
    data = np.loadtxt(raw_path)
    data = np.atleast_2d(data)[skip:]
    nsteps = data.shape[0]
    if nsteps < 4 * block_size:
        raise ValueError(
            f"only {nsteps} measurements left after skip={skip}; "
            f"need at least {4 * block_size} for block_size={block_size}"
        )

    observables, correlators = split_raw_columns(data)
    ncorr = correlators[CORRELATOR_NAMES[0]].shape[1]
    os.makedirs(outdir, exist_ok=True)

    if verbose:
        print(f"skip={skip}, nsteps={nsteps}, ncorr={ncorr}")

    means = {name: float(series.mean())
             for name, series in observables.items()}

    # ---- blocking scan on the primary observables -------------------
    kmax = max(st.max_block_exponent(nsteps), 0)
    k_values = 2 ** np.arange(kmax + 1)
    obs_sigmas = np.column_stack([
        st.blocking_scan(observables[name], kmax)[1]
        for name in OBSERVABLE_NAMES
    ])

    # ---- blocking scan of the jackknife error on the correlators ----
    check_points = st.select_check_points(ncorr)
    jack_sigmas = np.full((kmax + 1, len(CORRELATOR_NAMES) * len(check_points)),
                          np.nan)
    for ic, cname in enumerate(CORRELATOR_NAMES):
        source = observables[_CORRELATOR_SOURCE[ic]]
        for jn, sep in enumerate(check_points):
            col = ic * len(check_points) + jn
            for p, k in enumerate(k_values):
                if nsteps // int(k) < 4:
                    continue
                _, sigma = st.jackknife_connected(
                    correlators[cname][:, sep], source, int(k))
                jack_sigmas[p, col] = sigma

    files = {}

    header = "k  " + "  ".join(f"sigma_{n}" for n in OBSERVABLE_NAMES)
    path = os.path.join(outdir, "blocking_observables.dat")
    np.savetxt(path, np.column_stack([k_values, obs_sigmas]),
               header=header, fmt="%16.8e")
    files["blocking_observables"] = path

    header = "k  " + "  ".join(
        f"sigma_{c}_n{sep}" for c in CORRELATOR_NAMES for sep in check_points)
    path = os.path.join(outdir, "blocking_jackknife_correlators.dat")
    np.savetxt(path, np.column_stack([k_values, jack_sigmas]),
               header=header, fmt="%16.8e")
    files["blocking_jackknife_correlators"] = path

    # ---- connected correlators --------------------------------------
    columns = [np.arange(1, ncorr + 1)]
    header = "n"
    for ic, cname in enumerate(CORRELATOR_NAMES):
        source = observables[_CORRELATOR_SOURCE[ic]]
        mean, sigma = st.jackknife_connected_all_n(
            correlators[cname], source, block_size)
        columns.extend([mean, sigma])
        header += f"  {cname}_mean  {cname}_err"
    path = os.path.join(outdir, "connected_correlators.dat")
    np.savetxt(path, np.column_stack(columns), header=header, fmt="%16.8e")
    files["connected_correlators"] = path

    # ---- effective energy gaps --------------------------------------
    columns = [np.arange(1, ncorr)]
    header = "n"
    for ic, cname in enumerate(CORRELATOR_NAMES):
        source = observables[_CORRELATOR_SOURCE[ic]]
        gap, sigma = st.jackknife_energy_gap(
            correlators[cname], source, block_size, eta)
        columns.extend([gap, sigma])
        header += f"  DE_{cname}  DE_{cname}_err"
    path = os.path.join(outdir, "energy_gaps.dat")
    np.savetxt(path, np.column_stack(columns), header=header, fmt="%16.8e")
    files["energy_gaps"] = path

    # ---- autocorrelations -------------------------------------------
    lag_count = min(max_lag, nsteps)
    acfs = {name: st.autocorrelation(observables[name], lag_count)
            for name in OBSERVABLE_NAMES}
    tau_int = {name: st.integrated_time(acf) for name, acf in acfs.items()}

    path = os.path.join(outdir, "tau_int.dat")
    with open(path, "w") as handle:
        handle.write("# observable  tau_int\n")
        for name in OBSERVABLE_NAMES:
            handle.write(f"{name}  {tau_int[name]:.4f}\n")
    files["tau_int"] = path

    amplitude, tau_exp = st.fit_exponential_autocorrelation(
        acfs["y2"], fallback_tau=tau_int["y2"])
    path = os.path.join(outdir, "tau_exp_fit.dat")
    np.savetxt(
        path,
        np.column_stack([np.arange(lag_count), acfs["y2"]]),
        header=(f"lag  C_{{y^2}}\n"
                f"# fit_params: A={amplitude:.8e}  tau={tau_exp:.8e}"),
        fmt="%16.8e",
    )
    files["tau_exp_fit"] = path

    # ---- observable means and errors --------------------------------
    k_index = min(int(np.log2(block_size)), kmax)
    errors = {name: float(obs_sigmas[k_index, i])
              for i, name in enumerate(OBSERVABLE_NAMES)}

    path = os.path.join(outdir, "observables.dat")
    with open(path, "w") as handle:
        handle.write("# observable  mean  sigma_blocking\n")
        for name in OBSERVABLE_NAMES:
            handle.write(f"{name}  {means[name]:.8e}  {errors[name]:.8e}\n")
            if verbose:
                print(f"  {name} = {means[name]:.6f} +/- {errors[name]:.6f}")
    files["observables"] = path

    if verbose:
        print(f"tau_exp (from y^2) = {tau_exp:.2f}")
        print(f"results written to {outdir}")

    return {
        "nsteps": nsteps,
        "ncorr": ncorr,
        "means": means,
        "errors": errors,
        "tau_int": tau_int,
        "tau_exp": tau_exp,
        "files": files,
    }


def analyze_energy_only(raw_path, outdir, skip=0, block_size=None):
    """Analyse a chain recorded in ``energy_only`` mode.

    These runs store a single column, so only the mean energy and its
    blocking error can be extracted. A marker file is written so the GUI
    can tell such results apart from complete ones.

    Parameters
    ----------
    raw_path : str
        Path to the ``raw_energy_nt*.dat`` file.
    outdir : str
        Directory to write the results into; created if missing.
    skip : int, optional
        Leading measurements discarded as thermalisation.
    block_size : int, optional
        Block length. Defaults to the largest power of two that leaves at
        least four blocks, capped at :data:`DEFAULT_BLOCK`.

    Returns
    -------
    mean : float
        Mean energy.
    error : float
        Its blocking uncertainty.
    """
    data = np.loadtxt(raw_path).ravel()[skip:]
    if data.size < 8:
        raise ValueError(
            f"only {data.size} measurements left after skip={skip}")

    if block_size is None:
        kmax = max(st.max_block_exponent(data.size), 0)
        block_size = min(DEFAULT_BLOCK, 2 ** kmax)

    mean = float(data.mean())
    error = st.blocking_sigma(data, block_size)

    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "observables.dat"), "w") as handle:
        handle.write("# observable  mean  sigma_blocking\n")
        handle.write(f"E  {mean:.8e}  {error:.8e}\n")
    with open(os.path.join(outdir, "energy_only.marker"), "w") as handle:
        handle.write("Energy-only dataset\n")

    return mean, error
