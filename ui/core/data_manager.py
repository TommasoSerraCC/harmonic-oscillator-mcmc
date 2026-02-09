import numpy as np
import os
import re
import glob
import json
from datetime import datetime


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA_DIR = os.path.join(ROOT, 'data')
RESULTS_DIR = os.path.join(ROOT, 'results')
PLOTS_DIR = os.path.join(ROOT, 'plots')


def _parse_basedir(name):
    m = re.match(r'bhw(\d+)_nstep(\d+)', name)
    if m:
        return int(m.group(1)), int(m.group(2))
    return None


def _parse_resdir(name):
    m = re.match(r'nt(\d+)_therm(\d+)', name)
    if m:
        return int(m.group(1)), int(m.group(2))
    return None


def scan_data_sets():
    sets = []
    if not os.path.isdir(DATA_DIR):
        return sets
    for d in sorted(os.listdir(DATA_DIR)):
        p = _parse_basedir(d)
        if p:
            sets.append(p)
    return sets


def get_available_nt(bhw, nstep):
    path = os.path.join(DATA_DIR, f'bhw{bhw}_nstep{nstep}')
    if not os.path.isdir(path):
        return []
    nts = []
    for f in sorted(os.listdir(path)):
        m = re.match(r'raw_data_nt(\d+)\.dat', f)
        if m:
            nts.append(int(m.group(1)))
    return sorted(nts)


def scan_results():
    entries = []
    if not os.path.isdir(RESULTS_DIR):
        return entries
    for bd in sorted(os.listdir(RESULTS_DIR)):
        p = _parse_basedir(bd)
        if not p:
            continue
        bhw, nstep = p
        bdpath = os.path.join(RESULTS_DIR, bd)
        for rd in sorted(os.listdir(bdpath)):
            rp = _parse_resdir(rd)
            if rp:
                entries.append((bhw, nstep, rp[0], rp[1]))
    return entries


def get_available_therm(bhw, nstep, nt):
    bdpath = os.path.join(RESULTS_DIR, f'bhw{bhw}_nstep{nstep}')
    if not os.path.isdir(bdpath):
        return []
    therms = []
    for rd in os.listdir(bdpath):
        m = re.match(r'nt(\d+)_therm(\d+)', rd)
        if m and int(m.group(1)) == nt:
            therms.append(int(m.group(2)))
    return sorted(therms)


def _resdir(bhw, nstep, nt, therm):
    return os.path.join(RESULTS_DIR, f'bhw{bhw}_nstep{nstep}', f'nt{nt}_therm{therm}')


def _plotdir(bhw, nstep, nt, therm):
    d = os.path.join(PLOTS_DIR, f'bhw{bhw}_nstep{nstep}', f'nt{nt}_therm{therm}')
    os.makedirs(d, exist_ok=True)
    return d


def results_exist(bhw, nstep, nt, therm):
    return os.path.isdir(_resdir(bhw, nstep, nt, therm))


def load_observables(bhw, nstep, nt, therm):
    path = os.path.join(_resdir(bhw, nstep, nt, therm), 'observables.dat')
    obs = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            obs[parts[0]] = (float(parts[1]), float(parts[2]))
    return obs


def load_connected_correlators(bhw, nstep, nt, therm):
    path = os.path.join(_resdir(bhw, nstep, nt, therm), 'connected_correlators.dat')
    data = np.loadtxt(path)
    n = data[:, 0].astype(int)
    corr_names = ['yc', 'y2c', 'y3c', 'Ac']
    result = {'n': n}
    for i, cn in enumerate(corr_names):
        result[cn] = data[:, 1 + 2*i]
        result[cn + '_err'] = data[:, 2 + 2*i]
    return result


def load_energy_gaps(bhw, nstep, nt, therm):
    path = os.path.join(_resdir(bhw, nstep, nt, therm), 'energy_gaps.dat')
    data = np.loadtxt(path, ndmin=2)
    n = data[:, 0].astype(int)
    corr_names = ['yc', 'y2c', 'y3c', 'Ac']
    result = {'n': n}
    for i, cn in enumerate(corr_names):
        result[cn] = data[:, 1 + 2*i]
        result[cn + '_err'] = data[:, 2 + 2*i]
    return result


def load_blocking_observables(bhw, nstep, nt, therm):
    path = os.path.join(_resdir(bhw, nstep, nt, therm), 'blocking_observables.dat')
    data = np.loadtxt(path)
    obs_names = ['y', 'y2', 'y3', 'A', 'E']
    result = {'k': data[:, 0]}
    for i, name in enumerate(obs_names):
        result[name] = data[:, 1 + i]
    return result


def load_blocking_jackknife(bhw, nstep, nt, therm):
    path = os.path.join(_resdir(bhw, nstep, nt, therm), 'blocking_jackknife_correlators.dat')
    data = np.loadtxt(path)
    k = data[:, 0]
    n_check_vals = []
    with open(path) as f:
        header = f.readline()
    for token in header.split():
        if token.startswith('sigma_') and '_n' in token:
            nv = int(token.split('_n')[-1])
            if nv not in n_check_vals:
                n_check_vals.append(nv)
    corr_names = ['yc', 'y2c', 'y3c', 'Ac']
    npc = len(n_check_vals)
    result = {'k': k, 'n_check': n_check_vals}
    for ic, cn in enumerate(corr_names):
        result[cn] = {}
        for jn, nv in enumerate(n_check_vals):
            col = 1 + ic * npc + jn
            result[cn][nv] = data[:, col]
    return result


def load_tau_exp_fit(bhw, nstep, nt, therm):
    path = os.path.join(_resdir(bhw, nstep, nt, therm), 'tau_exp_fit.dat')
    data = np.loadtxt(path)
    A_exp, tau_exp = 1.0, 5.0
    with open(path) as f:
        for line in f:
            if 'fit_params' in line:
                for p in line.split():
                    if p.startswith('A='):
                        A_exp = float(p[2:])
                    elif p.startswith('tau='):
                        tau_exp = float(p[4:])
                break
    return {'lag': data[:, 0], 'acf': data[:, 1], 'A': A_exp, 'tau': tau_exp}


def load_correlation_length(bhw, nstep, nt, therm):
    path = os.path.join(_resdir(bhw, nstep, nt, therm), 'correlation_length.dat')
    if not os.path.isfile(path):
        return None
    result = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            result[parts[0]] = {
                'xi': float(parts[1]), 'xi_err': float(parts[2]),
                'DeltaE': float(parts[3]), 'DeltaE_err': float(parts[4])
            }
    return result


def load_tau_int(bhw, nstep, nt, therm):
    path = os.path.join(_resdir(bhw, nstep, nt, therm), 'tau_int.dat')
    result = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            result[parts[0]] = float(parts[1])
    return result


def collect_observables_vs_eta(bhw, nstep, therm=None):
    bdpath = os.path.join(RESULTS_DIR, f'bhw{bhw}_nstep{nstep}')
    if not os.path.isdir(bdpath):
        return None
    entries = []
    for rd in sorted(os.listdir(bdpath)):
        rp = _parse_resdir(rd)
        if not rp:
            continue
        nt_v, therm_v = rp
        if therm is not None and therm_v != therm:
            continue
        obs = load_observables(bhw, nstep, nt_v, therm_v)
        entries.append((nt_v, therm_v, obs))
    return entries


def load_fit_log(bhw, nstep, nt, therm):
    path = os.path.join(_resdir(bhw, nstep, nt, therm), 'fit_log.dat')
    if not os.path.isfile(path):
        return []
    entries = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            entries.append(line.strip())
    return entries


def append_fit_log(bhw, nstep, nt, therm, entry_str):
    rd = _resdir(bhw, nstep, nt, therm)
    os.makedirs(rd, exist_ok=True)
    path = os.path.join(rd, 'fit_log.dat')
    is_new = not os.path.isfile(path)
    with open(path, 'a') as f:
        if is_new:
            f.write("# timestamp  type  operator  n_min  n_max  value  error  chi2red  extra\n")
        f.write(entry_str + '\n')


def list_saved_plots(bhw, nstep, nt, therm):
    pd = os.path.join(PLOTS_DIR, f'bhw{bhw}_nstep{nstep}', f'nt{nt}_therm{therm}')
    if not os.path.isdir(pd):
        return []
    return sorted(glob.glob(os.path.join(pd, '*.png')))


def get_unanalyzed_raw(bhw, nstep, therm):
    raw_nts = set(get_available_nt(bhw, nstep))
    analyzed = set()
    bdpath = os.path.join(RESULTS_DIR, f'bhw{bhw}_nstep{nstep}')
    if os.path.isdir(bdpath):
        for rd in os.listdir(bdpath):
            m = re.match(r'nt(\d+)_therm(\d+)', rd)
            if m and int(m.group(2)) == therm:
                analyzed.add(int(m.group(1)))
    return sorted(raw_nts - analyzed)


def get_all_available_therms(bhw, nstep):
    bdpath = os.path.join(RESULTS_DIR, f'bhw{bhw}_nstep{nstep}')
    if not os.path.isdir(bdpath):
        return []
    therms = set()
    for rd in os.listdir(bdpath):
        m = re.match(r'nt(\d+)_therm(\d+)', rd)
        if m:
            therms.add(int(m.group(2)))
    return sorted(therms)


# ===== Energy vs Temperature (cross-bhw) =====

def get_available_etas_across_bhw():
    """Get all distinct eta values across all results."""
    etas = set()
    for bhw, nstep, nt, therm in scan_results():
        eta = round(float(bhw) / nt, 6)
        etas.add(eta)
    return sorted(etas)


def collect_energy_vs_temperature(eta_target, therm=None, tolerance=0.001):
    """Collect energy for fixed eta across different bhw."""
    results = []
    seen = set()
    for bhw, nstep, nt, therm_v in scan_results():
        eta = float(bhw) / nt
        if abs(eta - eta_target) > tolerance:
            continue
        if therm is not None and therm_v != therm:
            continue
        key = (bhw, nt, therm_v)
        if key in seen:
            continue
        seen.add(key)
        try:
            obs = load_observables(bhw, nstep, nt, therm_v)
        except Exception:
            continue
        if 'E' in obs:
            results.append({
                'bhw': bhw, 'nstep': nstep, 'nt': nt, 'therm': therm_v,
                'inv_bhw': 1.0 / bhw,
                'E': obs['E'][0], 'E_err': obs['E'][1],
                'eta': eta
            })
    results.sort(key=lambda x: x['inv_bhw'])
    return results


# ===== Fit save/load/delete =====

def _fits_dir(base_dir):
    d = os.path.join(base_dir, 'fits')
    os.makedirs(d, exist_ok=True)
    return d


def save_fit_entry(base_dir, fit_type, fit_data, fig):
    """Save a fit entry (JSON + PNG). Returns tag."""
    fd = _fits_dir(base_dir)
    ts = datetime.now().strftime('%Y%m%d_%H%M%S')
    tag = f'{fit_type}_{ts}'
    # Save JSON (convert numpy to float)
    clean = {}
    for k, v in fit_data.items():
        if isinstance(v, np.ndarray):
            clean[k] = v.tolist()
        elif isinstance(v, (np.floating,)):
            clean[k] = float(v)
        elif isinstance(v, (np.integer,)):
            clean[k] = int(v)
        else:
            clean[k] = v
    with open(os.path.join(fd, f'{tag}.json'), 'w') as f:
        json.dump(clean, f, indent=2)
    fig.savefig(os.path.join(fd, f'{tag}.png'), dpi=150, bbox_inches='tight')
    return tag


def list_fit_entries(base_dir):
    """List saved fit entry tags."""
    fd = os.path.join(base_dir, 'fits')
    if not os.path.isdir(fd):
        return []
    tags = []
    for f in sorted(os.listdir(fd)):
        if f.endswith('.json'):
            tags.append(f[:-5])
    return tags


def load_fit_entry(base_dir, tag):
    """Load a saved fit entry. Returns (data_dict, png_path_or_None)."""
    fd = os.path.join(base_dir, 'fits')
    json_path = os.path.join(fd, f'{tag}.json')
    png_path = os.path.join(fd, f'{tag}.png')
    data = None
    if os.path.isfile(json_path):
        with open(json_path) as f:
            data = json.load(f)
    return data, png_path if os.path.isfile(png_path) else None


def delete_fit_entry(base_dir, tag):
    """Delete a saved fit entry."""
    fd = os.path.join(base_dir, 'fits')
    for ext in ('.json', '.png'):
        p = os.path.join(fd, f'{tag}{ext}')
        if os.path.isfile(p):
            os.remove(p)


def energy_vs_temp_dir():
    d = os.path.join(RESULTS_DIR, 'energy_vs_temp')
    os.makedirs(d, exist_ok=True)
    return d


# ===== Energy-only data support =====

def get_available_energy_only_nt(bhw, nstep):
    """Scan for raw_energy_nt*.dat files (energy-only acquisitions)."""
    path = os.path.join(DATA_DIR, f'bhw{bhw}_nstep{nstep}')
    if not os.path.isdir(path):
        return []
    nts = []
    for f in sorted(os.listdir(path)):
        m = re.match(r'raw_energy_nt(\d+)\.dat', f)
        if m:
            nts.append(int(m.group(1)))
    return sorted(nts)


def is_energy_only(bhw, nstep, nt, therm):
    """Check if a result set is energy-only."""
    return os.path.isfile(os.path.join(_resdir(bhw, nstep, nt, therm), 'energy_only.marker'))


def get_unanalyzed_energy_only(bhw, nstep, therm):
    """Get energy-only raw files not yet analyzed."""
    raw_nts = set(get_available_energy_only_nt(bhw, nstep))
    analyzed = set()
    bdpath = os.path.join(RESULTS_DIR, f'bhw{bhw}_nstep{nstep}')
    if os.path.isdir(bdpath):
        for rd in os.listdir(bdpath):
            m = re.match(r'nt(\d+)_therm(\d+)', rd)
            if m and int(m.group(2)) == therm:
                # Consider analyzed if marker exists (energy-only result)
                if is_energy_only(bhw, nstep, int(m.group(1)), therm):
                    analyzed.add(int(m.group(1)))
    return sorted(raw_nts - analyzed)


def analyze_energy_only(bhw, nstep, nt, skip):
    """Analyze energy-only raw data: blocking for mean and uncertainty of E."""
    raw_path = os.path.join(DATA_DIR, f'bhw{bhw}_nstep{nstep}', f'raw_energy_nt{nt}.dat')
    data = np.loadtxt(raw_path)
    data = data[skip:]
    nsteps_eff = len(data)

    mean_E = data.mean()

    # Blocking to estimate uncertainty
    kmax = int(np.log2(nsteps_eff // 4))

    def blocking_sigma(x, k):
        nblocks = len(x) // k
        trimmed = x[:nblocks * k].reshape(nblocks, k)
        block_means = trimmed.mean(axis=1)
        return block_means.std(ddof=1) / np.sqrt(nblocks)

    block_size = min(1000, 2**kmax)
    err_E = blocking_sigma(data, block_size)

    # Save results
    outdir = _resdir(bhw, nstep, nt, skip)
    os.makedirs(outdir, exist_ok=True)

    with open(os.path.join(outdir, 'observables.dat'), 'w') as f:
        f.write("# observable  mean  sigma_blocking\n")
        f.write(f"E  {mean_E:.8e}  {err_E:.8e}\n")

    # Write marker
    with open(os.path.join(outdir, 'energy_only.marker'), 'w') as f:
        f.write("Energy-only dataset\n")

    return mean_E, err_E
