"""Tests for the end-to-end analysis in :mod:`analysis.pipeline`.

The synthetic chain has a known correlation length, so the pipeline is
checked on physics content and not only on the shape of its output.
"""

import os

import numpy as np
import pytest

from analysis import pipeline
from analysis.pipeline import CORRELATOR_NAMES, OBSERVABLE_NAMES
from helpers import make_raw_array

BLOCK = 50
MAX_LAG = 20

EXPECTED_FILES = [
    'observables.dat',
    'blocking_observables.dat',
    'blocking_jackknife_correlators.dat',
    'connected_correlators.dat',
    'energy_gaps.dat',
    'tau_int.dat',
    'tau_exp_fit.dat',
]


@pytest.fixture
def analysed(raw_dataset, tmp_path):
    """Run the pipeline once and return ``(summary, outdir, dataset)``."""
    outdir = tmp_path / 'out'
    summary = pipeline.analyze_dataset(
        raw_dataset['path'], str(outdir), raw_dataset['eta'],
        block_size=BLOCK, max_lag=MAX_LAG)
    return summary, outdir, raw_dataset


# --------------------------------------------------------------------
#  Column splitting
# --------------------------------------------------------------------

def test_split_raw_columns_separates_observables_and_correlators(rng):
    array = make_raw_array(rng, 50, 6, xi=2.0)
    observables, correlators = pipeline.split_raw_columns(array)

    assert set(observables) == set(OBSERVABLE_NAMES)
    assert set(correlators) == set(CORRELATOR_NAMES)
    assert observables['y'] == pytest.approx(array[:, 0])
    assert observables['E'] == pytest.approx(array[:, 4])
    for j, name in enumerate(CORRELATOR_NAMES):
        assert correlators[name].shape == (50, 6)
        assert correlators[name] == pytest.approx(array[:, 5 + j::4])


@pytest.mark.parametrize('ncols', [3, 8, 10])
def test_split_raw_columns_rejects_an_inconsistent_layout(ncols):
    with pytest.raises(ValueError, match='unexpected raw layout'):
        pipeline.split_raw_columns(np.zeros((10, ncols)))


# --------------------------------------------------------------------
#  Output files
# --------------------------------------------------------------------

def test_pipeline_writes_every_expected_file(analysed):
    _, outdir, _ = analysed
    for name in EXPECTED_FILES:
        assert (outdir / name).is_file(), f'{name} was not written'


def test_pipeline_reports_the_files_it_wrote(analysed):
    summary, outdir, _ = analysed
    for path in summary['files'].values():
        assert os.path.isfile(path)
        assert os.path.dirname(path) == str(outdir)


def test_pipeline_creates_a_missing_output_directory(raw_dataset, tmp_path):
    outdir = tmp_path / 'deeply' / 'nested'
    pipeline.analyze_dataset(raw_dataset['path'], str(outdir),
                             raw_dataset['eta'], block_size=BLOCK,
                             max_lag=MAX_LAG)
    assert (outdir / 'observables.dat').is_file()


# --------------------------------------------------------------------
#  Physics content
# --------------------------------------------------------------------

def test_reported_means_match_the_raw_columns(analysed):
    summary, _, dataset = analysed
    for i, name in enumerate(OBSERVABLE_NAMES):
        assert summary['means'][name] == pytest.approx(
            dataset['array'][:, i].mean())


def test_connected_correlators_follow_the_input_decay(analysed):
    _, outdir, dataset = analysed
    data = np.loadtxt(outdir / 'connected_correlators.dat')

    separations = data[:, 0]
    assert separations == pytest.approx(np.arange(1, dataset['ncorr'] + 1))

    expected = np.exp(-separations / dataset['xi'])
    for i in range(len(CORRELATOR_NAMES)):
        assert data[:, 1 + 2 * i] == pytest.approx(expected, rel=0.02)
        assert np.all(data[:, 2 + 2 * i] > 0.0)


def test_energy_gaps_recover_the_input_correlation_length(analysed):
    _, outdir, dataset = analysed
    data = np.loadtxt(outdir / 'energy_gaps.dat')

    assert data.shape[0] == dataset['ncorr'] - 1
    for i in range(len(CORRELATOR_NAMES)):
        assert data[:, 1 + 2 * i] == pytest.approx(
            dataset['expected_gap'], rel=0.02)


def test_errors_are_positive_and_smaller_than_the_signal(analysed):
    summary, _, _ = analysed
    for name in OBSERVABLE_NAMES:
        assert summary['errors'][name] > 0.0
    assert summary['errors']['E'] < 0.01


# --------------------------------------------------------------------
#  Table shapes
# --------------------------------------------------------------------

def test_blocking_tables_have_one_row_per_block_size(analysed):
    _, outdir, dataset = analysed
    expected_rows = int(np.log2(dataset['nsteps'] // 4)) + 1

    observables = np.loadtxt(outdir / 'blocking_observables.dat')
    assert observables.shape == (expected_rows, 1 + len(OBSERVABLE_NAMES))
    assert list(observables[:, 0]) == [2 ** p for p in range(expected_rows)]

    jackknife = np.loadtxt(outdir / 'blocking_jackknife_correlators.dat')
    assert jackknife.shape[0] == expected_rows
    assert (jackknife.shape[1] - 1) % len(CORRELATOR_NAMES) == 0


def test_autocorrelation_table_has_the_requested_lags(analysed):
    _, outdir, _ = analysed
    data = np.loadtxt(outdir / 'tau_exp_fit.dat')
    assert data.shape == (MAX_LAG, 2)
    assert data[0, 1] == pytest.approx(1.0)


def test_tau_int_lists_every_observable(analysed):
    summary, outdir, _ = analysed
    text = (outdir / 'tau_int.dat').read_text().splitlines()
    assert text[0].startswith('#')
    assert [line.split()[0] for line in text[1:]] == OBSERVABLE_NAMES
    assert set(summary['tau_int']) == set(OBSERVABLE_NAMES)


# --------------------------------------------------------------------
#  Thermalisation and error handling
# --------------------------------------------------------------------

def test_skip_discards_the_leading_measurements(raw_dataset, tmp_path):
    summary = pipeline.analyze_dataset(
        raw_dataset['path'], str(tmp_path / 'skipped'), raw_dataset['eta'],
        skip=500, block_size=BLOCK, max_lag=MAX_LAG)

    assert summary['nsteps'] == raw_dataset['nsteps'] - 500
    assert summary['means']['E'] == pytest.approx(
        raw_dataset['array'][500:, 4].mean())


def test_pipeline_refuses_a_chain_that_is_too_short(raw_dataset, tmp_path):
    with pytest.raises(ValueError, match='need at least'):
        pipeline.analyze_dataset(
            raw_dataset['path'], str(tmp_path / 'short'), raw_dataset['eta'],
            skip=1950, block_size=BLOCK, max_lag=MAX_LAG)


# --------------------------------------------------------------------
#  Energy-only mode
# --------------------------------------------------------------------

def test_energy_only_reports_mean_and_error(tmp_path, rng):
    values = rng.normal(0.5, 0.02, size=4000)
    raw = tmp_path / 'raw_energy_nt10.dat'
    np.savetxt(raw, values, fmt='%20.12e')

    outdir = tmp_path / 'eout'
    mean, error = pipeline.analyze_energy_only(str(raw), str(outdir))

    assert mean == pytest.approx(values.mean())
    assert error == pytest.approx(values.std(ddof=1) / np.sqrt(values.size),
                                  rel=0.5)
    assert (outdir / 'energy_only.marker').is_file()

    written = (outdir / 'observables.dat').read_text().splitlines()
    assert written[1].split()[0] == 'E'
    assert float(written[1].split()[1]) == pytest.approx(mean)


def test_energy_only_honours_skip(tmp_path, rng):
    values = np.concatenate([np.full(200, 5.0), rng.normal(0.5, 0.02, 3000)])
    raw = tmp_path / 'raw_energy_nt10.dat'
    np.savetxt(raw, values, fmt='%20.12e')

    mean, _ = pipeline.analyze_energy_only(str(raw), str(tmp_path / 'o'),
                                           skip=200)
    assert mean == pytest.approx(values[200:].mean())


def test_energy_only_refuses_an_empty_chain(tmp_path):
    raw = tmp_path / 'raw_energy_nt10.dat'
    np.savetxt(raw, np.arange(10.0))
    with pytest.raises(ValueError, match='measurements left'):
        pipeline.analyze_energy_only(str(raw), str(tmp_path / 'o'), skip=9)
