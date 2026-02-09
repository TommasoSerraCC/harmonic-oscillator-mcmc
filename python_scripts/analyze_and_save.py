import numpy as np
import argparse
import os
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser()
parser.add_argument('--nt', type=int, required=True)
parser.add_argument('--skip', type=int, default=0, help='thermalization steps to skip')
parser.add_argument('--bhw', type=int, required=True)
parser.add_argument('--nstep', type=int, default=1000000)
args = parser.parse_args()

nt = args.nt
skip = args.skip
ncorr = nt // 2
bhw = float(args.bhw)
eta = bhw / nt

basedir = f'bhw{args.bhw}_nstep{args.nstep}'

# ========== Load data ==========
data = np.loadtxt(f'data/{basedir}/raw_data_nt{nt}.dat')
data = data[skip:]
nsteps = data.shape[0]
print(f"nt={nt}, skip={skip}, nsteps={nsteps}")

y, y2, y3, A, E = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]
corr_data = data[:, 5:]
yc  = corr_data[:, 0::4]  # shape (nsteps, ncorr)
y2c = corr_data[:, 1::4]
y3c = corr_data[:, 2::4]
Ac  = corr_data[:, 3::4]

obs_names = ['y', 'y2', 'y3', 'A', 'E']
obs_arrays = [y, y2, y3, A, E]
obs_means = np.array([d.mean() for d in obs_arrays])
for name, m in zip(obs_names, obs_means):
    print(f"  <{name}> = {m:.6f}")

# ========== Output directory ==========
outdir = f'results/{basedir}/nt{nt}_therm{skip}'
os.makedirs(outdir, exist_ok=True)

# ========== Blocking for simple observables ==========
# blocking: divide data in blocks of size k, compute block means,
# then std of block means / sqrt(nblocks) = sigma of the mean
# We save sigma(mean) vs k for k = 2^0, 2^1, ..., 2^kmax

def blocking_all_k(x, kmax):
    """Returns array of sigma(mean) for k = 2^0 .. 2^kmax, fully vectorized."""
    n = len(x)
    sigmas = np.zeros(kmax + 1)
    for p in range(kmax + 1):
        k = 2**p
        nblocks = n // k
        trimmed = x[:nblocks * k].reshape(nblocks, k)
        block_means = trimmed.mean(axis=1)
        sigmas[p] = block_means.std(ddof=1) / np.sqrt(nblocks)
    return sigmas

kmax = int(np.log2(nsteps // 4))  # leave at least 4 blocks
k_values = 2**np.arange(kmax + 1)

print(f"\nBlocking: kmax={kmax}, max block size={2**kmax}")

obs_block_sigmas = np.zeros((kmax + 1, 5))
for i, d in enumerate(obs_arrays):
    obs_block_sigmas[:, i] = blocking_all_k(d, kmax)

# ========== Jackknife for connected correlators ==========
# Connected correlator: C_conn(n) = <O(0)O(n)> - <O>^2
# Jackknife: remove block j, recompute C_conn on the rest, collect, compute std

def jackknife_connected(corr_col, obs_col, block_size):
    """Jackknife estimate of mean and std for a single connected correlator.
    corr_col: raw correlator values, shape (nsteps,)
    obs_col:  raw observable values, shape (nsteps,)
    block_size: int
    Returns: (mean_jack, sigma_jack)
    """
    n = len(corr_col)
    nblocks = n // block_size
    nt = nblocks * block_size
    corr_t = corr_col[:nt].reshape(nblocks, block_size)
    obs_t = obs_col[:nt].reshape(nblocks, block_size)
    
    total_corr = corr_t.sum(axis=1)   # sum per block
    total_obs = obs_t.sum(axis=1)
    
    sum_corr_all = total_corr.sum()
    sum_obs_all = total_obs.sum()
    
    # leave-one-block-out means
    corr_jack = (sum_corr_all - total_corr) / (nt - block_size)
    obs_jack = (sum_obs_all - total_obs) / (nt - block_size)
    
    conn_jack = corr_jack - obs_jack**2
    
    mean_jack = conn_jack.mean()
    sigma_jack = np.sqrt((nblocks - 1) / nblocks * np.sum((conn_jack - mean_jack)**2))
    return mean_jack, sigma_jack


def jackknife_connected_all_n(corr_matrix, obs_col, block_size):
    """Jackknife for all n values at once. 
    corr_matrix: shape (nsteps, ncorr)
    obs_col: shape (nsteps,)
    Returns: means (ncorr,), sigmas (ncorr,)
    """
    n = len(obs_col)
    nblocks = n // block_size
    nt = nblocks * block_size
    nc = corr_matrix.shape[1]
    
    corr_t = corr_matrix[:nt].reshape(nblocks, block_size, nc)
    obs_t = obs_col[:nt].reshape(nblocks, block_size)
    
    total_corr = corr_t.sum(axis=1)   # (nblocks, nc)
    total_obs = obs_t.sum(axis=1)     # (nblocks,)
    
    sum_corr_all = total_corr.sum(axis=0)  # (nc,)
    sum_obs_all = total_obs.sum()
    
    # leave-one-block-out
    corr_jack = (sum_corr_all[None, :] - total_corr) / (nt - block_size)  # (nblocks, nc)
    obs_jack = (sum_obs_all - total_obs) / (nt - block_size)               # (nblocks,)
    
    conn_jack = corr_jack - obs_jack[:, None]**2  # (nblocks, nc)
    
    mean_jack = conn_jack.mean(axis=0)
    sigma_jack = np.sqrt((nblocks - 1) / nblocks * np.sum((conn_jack - mean_jack[None, :])**2, axis=0))
    return mean_jack, sigma_jack


# Choose 3 representative n values for the adaptive blocking scan
n_check = sorted(set([1, ncorr // 4, ncorr // 2]))
n_check = [min(v, ncorr - 1) for v in n_check]
n_check = sorted(set(n_check))
if len(n_check) < 3 and ncorr > 3:
    candidates = [1, ncorr // 3, 2 * ncorr // 3]
    n_check = sorted(set([min(v, ncorr - 1) for v in candidates]))
n_check = n_check[:3]
print(f"Jackknife blocking scan at n = {n_check}")

corr_names = ['yc', 'y2c', 'y3c', 'Ac']
corr_matrices = [yc, y2c, y3c, Ac]
obs_for_corr = [y, y2, y3, A]  # the simple obs matching each correlator

# Jackknife blocking scan: 4 correlators x 3 n values = 12 columns
jack_block_sigmas = np.zeros((kmax + 1, len(corr_names) * len(n_check)))

for ic, (cmat, obs_col) in enumerate(zip(corr_matrices, obs_for_corr)):
    for jn, nv in enumerate(n_check):
        col_idx = ic * len(n_check) + jn
        for p in range(kmax + 1):
            k = 2**p
            if nsteps // k < 4:
                jack_block_sigmas[p, col_idx] = np.nan
                continue
            _, sig = jackknife_connected(cmat[:, nv], obs_col, k)
            jack_block_sigmas[p, col_idx] = sig

# ========== Save blocking results ==========
# File: blocking_observables.dat
# columns: k  sigma_y  sigma_y2  sigma_y3  sigma_A  sigma_E
header_block = "k  sigma_y  sigma_y2  sigma_y3  sigma_A  sigma_E"
block_out = np.column_stack([k_values, obs_block_sigmas])
np.savetxt(f'{outdir}/blocking_observables.dat', block_out, header=header_block, fmt='%16.8e')

# File: blocking_jackknife_correlators.dat
cols_jack = []
for cn in corr_names:
    for nv in n_check:
        cols_jack.append(f"sigma_{cn}_n{nv}")
header_jack = "k  " + "  ".join(cols_jack)
jack_out = np.column_stack([k_values, jack_block_sigmas])
np.savetxt(f'{outdir}/blocking_jackknife_correlators.dat', jack_out, header=header_jack, fmt='%16.8e')

print("Blocking files saved.")

# ========== Jackknife with default block_size=1000 for all n ==========
default_block = 1000

corr_conn_means = {}
corr_conn_sigmas = {}
for cn, cmat, obs_col in zip(corr_names, corr_matrices, obs_for_corr):
    m, s = jackknife_connected_all_n(cmat, obs_col, default_block)
    corr_conn_means[cn] = m
    corr_conn_sigmas[cn] = s

# Save connected correlators: n  yc_mean  yc_err  y2c_mean  y2c_err  ...
n_arr = np.arange(1, ncorr + 1)
cols = [n_arr]
header_corr = "n"
for cn in corr_names:
    cols.append(corr_conn_means[cn])
    cols.append(corr_conn_sigmas[cn])
    header_corr += f"  {cn}_mean  {cn}_err"
corr_out = np.column_stack(cols)
np.savetxt(f'{outdir}/connected_correlators.dat', corr_out, header=header_corr, fmt='%16.8e')

# ========== Energy gaps with error propagation ==========
# DeltaE(n) = log(C(n) / C(n+1)) / eta
# Error propagation: sigma_DE = (1/eta) * sqrt( (sigC_n/C_n)^2 + (sigC_{n+1}/C_{n+1})^2 )
# (C_n and C_{n+1} are correlated, but we use jackknife directly on the ratio for better estimate)

def jackknife_energy_gap(corr_matrix, obs_col, block_size, eta):
    """Jackknife for energy gaps log(C(n)/C(n+1))/eta for all n."""
    n = len(obs_col)
    nblocks = n // block_size
    nt = nblocks * block_size
    nc = corr_matrix.shape[1]
    
    corr_t = corr_matrix[:nt].reshape(nblocks, block_size, nc)
    obs_t = obs_col[:nt].reshape(nblocks, block_size)
    
    total_corr = corr_t.sum(axis=1)
    total_obs = obs_t.sum(axis=1)
    sum_corr_all = total_corr.sum(axis=0)
    sum_obs_all = total_obs.sum()
    
    corr_jack = (sum_corr_all[None, :] - total_corr) / (nt - block_size)
    obs_jack = (sum_obs_all - total_obs) / (nt - block_size)
    conn_jack = corr_jack - obs_jack[:, None]**2  # (nblocks, nc)
    
    # energy gap for each jackknife sample
    ratio_jack = conn_jack[:, :-1] / conn_jack[:, 1:]  # (nblocks, nc-1)
    # protect against log of negative/zero
    ratio_jack = np.where(ratio_jack > 0, ratio_jack, np.nan)
    gap_jack = np.log(ratio_jack) / eta  # (nblocks, nc-1)
    
    mean_gap = np.nanmean(gap_jack, axis=0)
    sigma_gap = np.sqrt((nblocks - 1) / nblocks * np.nansum((gap_jack - mean_gap[None, :])**2, axis=0))
    return mean_gap, sigma_gap

n_gap = np.arange(1, ncorr)
gap_cols = [n_gap]
header_gap = "n"
for cn, cmat, obs_col in zip(corr_names, corr_matrices, obs_for_corr):
    mg, sg = jackknife_energy_gap(cmat, obs_col, default_block, eta)
    gap_cols.append(mg)
    gap_cols.append(sg)
    header_gap += f"  DE_{cn}  DE_{cn}_err"

gap_out = np.column_stack(gap_cols)
np.savetxt(f'{outdir}/energy_gaps.dat', gap_out, header=header_gap, fmt='%16.8e')
print("Energy gaps saved.")

# ========== Autocorrelation ==========
def autocorr(x, max_lag):
    n = len(x)
    xm = x - x.mean()
    var = np.dot(xm, xm)
    acf = np.array([np.dot(xm[:n-t], xm[t:]) / var * n / (n - t) for t in range(max_lag)])
    return acf

max_lag = 50
tau_ints = np.zeros(5)
C_obs = {}  # Autocorrelation functions: C_obs[name] = C_F(k)
for i, (name, d) in enumerate(zip(obs_names, obs_arrays)):
    C_obs[name] = autocorr(d, max_lag)
    tau_ints[i] = C_obs[name][1:].sum()

# Save tau_int
tau_out = np.column_stack([obs_names, tau_ints])
with open(f'{outdir}/tau_int.dat', 'w') as f:
    f.write("# observable  tau_int\n")
    for name, ti in zip(obs_names, tau_ints):
        f.write(f"{name}  {ti:.4f}\n")

# Fit tau_exp on y2
def exp_decay(t, A, tau):
    return A * np.exp(-t / tau)

C_y2 = C_obs['y2']
idx_zero = max(np.argmax(C_y2[1:] < 0.01) + 1, 10)
t_fit = np.arange(1, idx_zero)
try:
    popt, _ = curve_fit(exp_decay, t_fit, C_y2[1:idx_zero], p0=[1, 5], bounds=([0, 0.1], [np.inf, 100]))
    A_exp, tau_exp = popt
except:
    A_exp, tau_exp = 1.0, tau_ints[1]

print(f"  tau_exp (fit y$^2$) = {tau_exp:.2f}")

# Save tau_exp fit data: lag, C_{y^2} (first max_lag points), and fit params
lag_save = np.arange(max_lag)
C_y2_save = C_y2[:max_lag]
np.savetxt(f'{outdir}/tau_exp_fit.dat', np.column_stack([lag_save, C_y2_save]),
           header=f"lag  C_{{y^2}}\n# fit_params: A={A_exp:.8e}  tau={tau_exp:.8e}", fmt='%16.8e')

# ========== Save observable means + blocking errors ==========
# Use the blocking plateau value: take the sigma at a reasonable block size
# We pick block_size = 1000 (2^~10) from the blocking table
k_default_idx = min(int(np.log2(default_block)), kmax)
obs_errs = obs_block_sigmas[k_default_idx, :]

with open(f'{outdir}/observables.dat', 'w') as f:
    f.write("# observable  mean  sigma_blocking\n")
    for name, m, s in zip(obs_names, obs_means, obs_errs):
        f.write(f"{name}  {m:.8e}  {s:.8e}\n")
        print(f"  {name} = {m:.6f} +/- {s:.6f}")

print(f"\nAll results saved in {outdir}/")
