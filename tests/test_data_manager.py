"""Tests for :mod:`ui.core.data_manager`.

These exercise the contract between the analysis pipeline and the GUI:
the directory naming conventions, the loaders that parse every result
table, and the saved-fit archive. A real analysis is run into a
temporary tree, so a change to either side that breaks the other is
caught here.
"""

import numpy as np
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from analysis import pipeline
from helpers import make_raw_array
from ui.core import data_manager as dm

BHW = 10
NSTEP = 2000
NT = 16
THERM = 0


@pytest.fixture
def workspace(tmp_path, monkeypatch, rng):
    """A temporary data/results/plots tree with one analysed dataset."""
    data_dir = tmp_path / 'data'
    results_dir = tmp_path / 'results'
    plots_dir = tmp_path / 'plots'

    monkeypatch.setattr(dm, 'DATA_DIR', str(data_dir))
    monkeypatch.setattr(dm, 'RESULTS_DIR', str(results_dir))
    monkeypatch.setattr(dm, 'PLOTS_DIR', str(plots_dir))

    raw_dir = data_dir / f'bhw{BHW}_nstep{NSTEP}'
    raw_dir.mkdir(parents=True)
    raw_path = raw_dir / f'raw_data_nt{NT}.dat'
    np.savetxt(raw_path, make_raw_array(rng, NSTEP, NT // 2, xi=3.0),
               fmt='%20.12e')

    return {
        'root': tmp_path,
        'raw_path': raw_path,
        'eta': float(BHW) / NT,
        'outdir': results_dir / f'bhw{BHW}_nstep{NSTEP}'
                  / f'nt{NT}_therm{THERM}',
    }


def _analyse(workspace):
    pipeline.analyze_dataset(str(workspace['raw_path']),
                             str(workspace['outdir']), workspace['eta'],
                             block_size=50, max_lag=20)


# --------------------------------------------------------------------
#  Name parsing
# --------------------------------------------------------------------

@pytest.mark.parametrize('name, expected', [
    ('bhw10_nstep1000000', (10, 1000000)),
    ('bhw3_nstep500000', (3, 500000)),
    ('not_a_dataset', None),
    ('bhwX_nstep10', None),
])
def test_parse_basedir(name, expected):
    assert dm._parse_basedir(name) == expected


@pytest.mark.parametrize('name, expected', [
    ('nt200_therm50000', (200, 50000)),
    ('nt4_therm0', (4, 0)),
    ('nt200', None),
])
def test_parse_resdir(name, expected):
    assert dm._parse_resdir(name) == expected


# --------------------------------------------------------------------
#  Discovery
# --------------------------------------------------------------------

def test_scan_data_sets_finds_the_raw_directory(workspace):
    assert dm.scan_data_sets() == [(BHW, NSTEP)]


def test_get_available_nt_lists_the_raw_files(workspace):
    assert dm.get_available_nt(BHW, NSTEP) == [NT]


def test_scan_data_sets_is_empty_without_a_data_directory(tmp_path,
                                                          monkeypatch):
    monkeypatch.setattr(dm, 'DATA_DIR', str(tmp_path / 'missing'))
    assert dm.scan_data_sets() == []


def test_unanalysed_raw_files_are_reported_until_analysed(workspace):
    assert dm.get_unanalyzed_raw(BHW, NSTEP, THERM) == [NT]
    _analyse(workspace)
    assert dm.get_unanalyzed_raw(BHW, NSTEP, THERM) == []


def test_results_are_discovered_after_the_analysis(workspace):
    assert dm.scan_results() == []
    _analyse(workspace)
    assert dm.scan_results() == [(BHW, NSTEP, NT, THERM)]
    assert dm.results_exist(BHW, NSTEP, NT, THERM)
    assert dm.get_available_therm(BHW, NSTEP, NT) == [THERM]
    assert dm.get_all_available_therms(BHW, NSTEP) == [THERM]


# --------------------------------------------------------------------
#  Loaders
# --------------------------------------------------------------------

def test_load_observables_returns_value_and_error(workspace):
    _analyse(workspace)
    obs = dm.load_observables(BHW, NSTEP, NT, THERM)

    assert set(obs) == {'y', 'y2', 'y3', 'A', 'E'}
    for mean, error in obs.values():
        assert np.isfinite(mean)
        assert error > 0.0


def test_load_connected_correlators_exposes_every_correlator(workspace):
    _analyse(workspace)
    data = dm.load_connected_correlators(BHW, NSTEP, NT, THERM)

    assert list(data['n']) == list(range(1, NT // 2 + 1))
    for name in ['yc', 'y2c', 'y3c', 'Ac']:
        assert data[name].shape == (NT // 2,)
        assert np.all(data[name + '_err'] > 0.0)


def test_load_energy_gaps_has_one_row_fewer(workspace):
    _analyse(workspace)
    data = dm.load_energy_gaps(BHW, NSTEP, NT, THERM)
    assert data['n'].size == NT // 2 - 1
    assert np.all(np.isfinite(data['yc']))


def test_load_blocking_observables(workspace):
    _analyse(workspace)
    data = dm.load_blocking_observables(BHW, NSTEP, NT, THERM)
    assert data['k'][0] == 1
    for name in ['y', 'y2', 'y3', 'A', 'E']:
        assert data[name].shape == data['k'].shape


def test_load_blocking_jackknife_recovers_the_check_points(workspace):
    _analyse(workspace)
    data = dm.load_blocking_jackknife(BHW, NSTEP, NT, THERM)

    assert data['n_check'] == [1, 2, 4]
    for name in ['yc', 'y2c', 'y3c', 'Ac']:
        assert set(data[name]) == set(data['n_check'])
        for series in data[name].values():
            assert series.shape == data['k'].shape


def test_load_tau_int_and_tau_exp_fit(workspace):
    _analyse(workspace)

    tau_int = dm.load_tau_int(BHW, NSTEP, NT, THERM)
    assert set(tau_int) == {'y', 'y2', 'y3', 'A', 'E'}

    fit = dm.load_tau_exp_fit(BHW, NSTEP, NT, THERM)
    assert fit['acf'][0] == pytest.approx(1.0)
    assert fit['tau'] > 0.0
    assert fit['lag'].size == fit['acf'].size


def test_collect_observables_vs_eta(workspace):
    _analyse(workspace)
    entries = dm.collect_observables_vs_eta(BHW, NSTEP, THERM)
    assert len(entries) == 1
    nt_value, therm_value, obs = entries[0]
    assert (nt_value, therm_value) == (NT, THERM)
    assert 'E' in obs


# --------------------------------------------------------------------
#  Energy-only datasets
# --------------------------------------------------------------------

def test_energy_only_dataset_is_flagged_and_analysed(workspace, rng):
    raw = workspace['raw_path'].parent / 'raw_energy_nt8.dat'
    np.savetxt(raw, rng.normal(0.5, 0.02, size=4000), fmt='%20.12e')

    assert dm.get_available_energy_only_nt(BHW, NSTEP) == [8]
    assert dm.get_unanalyzed_energy_only(BHW, NSTEP, THERM) == [8]

    mean, error = dm.analyze_energy_only(BHW, NSTEP, 8, THERM)
    assert mean == pytest.approx(0.5, abs=0.01)
    assert error > 0.0

    assert dm.is_energy_only(BHW, NSTEP, 8, THERM)
    assert dm.get_unanalyzed_energy_only(BHW, NSTEP, THERM) == []


def test_energy_vs_temperature_collects_matching_eta(workspace):
    _analyse(workspace)
    etas = dm.get_available_etas_across_bhw()
    assert etas == [pytest.approx(float(BHW) / NT)]

    points = dm.collect_energy_vs_temperature(float(BHW) / NT)
    assert len(points) == 1
    assert points[0]['bhw'] == BHW
    assert points[0]['inv_bhw'] == pytest.approx(1.0 / BHW)


# --------------------------------------------------------------------
#  Saved-fit archive
# --------------------------------------------------------------------

def test_fit_entries_round_trip(tmp_path):
    base = tmp_path / 'archive'
    base.mkdir()
    payload = {'a': 0.5, 'a_err': 0.01, 'points': np.arange(3),
               'n': np.int64(7), 'x': np.float64(1.5)}

    figure = Figure()
    FigureCanvasAgg(figure)
    tag = dm.save_fit_entry(str(base), 'continuum_E', payload, figure)

    assert dm.list_fit_entries(str(base)) == [tag]

    loaded, png = dm.load_fit_entry(str(base), tag)
    assert loaded['a'] == 0.5
    assert loaded['points'] == [0, 1, 2]
    assert loaded['n'] == 7
    assert png is not None

    dm.delete_fit_entry(str(base), tag)
    assert dm.list_fit_entries(str(base)) == []


def test_list_fit_entries_without_an_archive(tmp_path):
    assert dm.list_fit_entries(str(tmp_path / 'nothing')) == []
