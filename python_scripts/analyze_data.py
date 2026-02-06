import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser()
parser.add_argument('--nt', type=int, required=True)
parser.add_argument('--skip', type=int, default=0, help='steps to skip (thermalization)')
args = parser.parse_args()

nt = args.nt
skip = args.skip
ncorr = nt // 2
ncols = 5 + 4 * ncorr
bhw = 10.0
eta = bhw / nt

data = np.loadtxt(f'data/raw_data_nt{nt}.dat')
data = data[skip:]
nsteps = data.shape[0]

print(f"nt={nt}, skip={skip}, nsteps={nsteps}")

y, y2, y3, A, E = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]

corr_data = data[:, 5:]
yc = corr_data[:, 0::4]
y2c = corr_data[:, 1::4]
y3c = corr_data[:, 2::4]
Ac = corr_data[:, 3::4]

obs_names = ['y', 'y2', 'y3', 'A', 'E']
obs_data = [y, y2, y3, A, E]
for name, d in zip(obs_names, obs_data):
    print(f"  <{name}> = {d.mean():.6f}")

corr_means = {
    'yc': yc.mean(axis=0) - y.mean()**2,
    'y2c': y2c.mean(axis=0) - y2.mean()**2,
    'y3c': y3c.mean(axis=0) - y3.mean()**2,
    'Ac': Ac.mean(axis=0) - A.mean()**2
}

# Plot 1: osservabili (primi 10000 step dopo skip)
fig1, axes1 = plt.subplots(5, 1, figsize=(10, 10), sharex=True)
for ax, name, d in zip(axes1, obs_names, obs_data):
    ax.plot(d[:10000], lw=0.5)
    ax.set_ylabel(name)
    ax.axhline(d.mean(), color='r', ls='--', lw=1)
axes1[-1].set_xlabel('step MCMC')
fig1.suptitle(f'Osservabili (nt={nt}, skip={skip})')
fig1.tight_layout()
fig1.savefig(f'plots/obs_nt{nt}_skip{skip}.png', dpi=150)

# Plot 2: energy gaps
fig2, ax2 = plt.subplots(figsize=(10, 6))
n_vals = np.arange(1, ncorr)
corr_labels = ['yc', 'y2c', 'y3c', 'Ac']
corr_arrays = [corr_means['yc'], corr_means['y2c'], corr_means['y3c'], corr_means['Ac']]
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

for C, label, col in zip(corr_arrays, corr_labels, colors):
    m_eff = np.log(C[:-1] / C[1:]) / eta
    ax2.plot(n_vals, m_eff, marker='o', ms=3, color=col, label=label)

ax2.set_xlabel('n')
ax2.set_ylabel(r'$\Delta E(n)$')
ax2.set_ylim(0, 5)
for yline in (1.0, 2.0, 3.0):
    ax2.axhline(yline, color='k', ls='--', lw=0.5)
fig2.suptitle(f'Energy gaps ($N_{{t}}={nt}$)')
fig2.tight_layout()
fig2.savefig(f'plots/meff_nt{nt}_skip{skip}.png', dpi=150)

# Autocorrelazione
def autocorr(x, max_lag):
    n = len(x)
    x = x - x.mean()
    var = np.sum(x**2)
    return np.array([np.sum(x[:n-t] * x[t:]) / var for t in range(max_lag)])

def tau_int(rho):
    tau = 0.0
    for k in range(1, len(rho)//2):
        tau += rho[k]
    return tau

max_lag = 500
print("\n=== Autocorrelazione ===")
tau_ints, rhos = [], []
for name, d in zip(obs_names, obs_data):
    rho = autocorr(d, max_lag)
    rhos.append(rho)
    ti = tau_int(rho)
    tau_ints.append(ti)
    print(f"  tau_int({name}) = {ti:.2f}")

def exp_decay(t, A, tau):
    return A * np.exp(-t / tau)

rho_y2 = rhos[1]
idx_zero = max(np.argmax(rho_y2[1:] < 0.01) + 1, 10)
t_fit = np.arange(1, idx_zero)
try:
    popt, _ = curve_fit(exp_decay, t_fit, rho_y2[1:idx_zero], p0=[1, 5], bounds=([0, 0.1], [np.inf, 100]))
    tau_exp = popt[1]
except:
    tau_exp = tau_ints[1]
print(f"  tau_exp (fit y2) = {tau_exp:.2f}")

# Plot 3: autocorrelazione
fig3, axes3 = plt.subplots(2, 3, figsize=(12, 7))
axes3 = axes3.flatten()
for i, (name, rho) in enumerate(zip(obs_names, rhos)):
    ax = axes3[i]
    ax.plot(rho[:100], lw=0.8)
    ax.axvline(tau_ints[i], color='r', ls='--', label=f'tau_int={tau_ints[i]:.1f}')
    ax.set_xlabel('lag')
    ax.set_ylabel(r'$\rho(t)$')
    ax.set_title(name)
    ax.legend(fontsize=8)

ax = axes3[5]
ax.plot(np.arange(1, 51), rho_y2[1:51], 'b-', lw=0.8, label='dati y2')
t_plot = np.linspace(1, 50, 100)
ax.plot(t_plot, exp_decay(t_plot, popt[0] if 'popt' in dir() else 1, tau_exp), 'r--', lw=1.5, label=f'tau_exp={tau_exp:.1f}')
ax.set_xlabel('lag')
ax.set_ylabel(r'$\rho(t)$')
ax.set_title('Fit esponenziale (y2)')
ax.legend(fontsize=8)

fig3.suptitle(f'Autocorrelazione (nt={nt}, skip={skip})')
fig3.tight_layout()
fig3.savefig(f'plots/autocorr_nt{nt}_skip{skip}.png', dpi=150)

plt.show()
