import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--nt', type=int, required=True)
parser.add_argument('--skip', type=int, default=0, help='thermalization steps')
parser.add_argument('--gap_nmax', type=int, default=0, help='max n for energy gap plot (0=auto)')
args = parser.parse_args()

nt = args.nt
skip = args.skip
ncorr = nt // 2
bhw = 10.0
eta = bhw / nt

resdir = f'results/nt{nt}_therm{skip}'
plotdir = f'plots/nt{nt}_therm{skip}'
os.makedirs(plotdir, exist_ok=True)

# ========== 1. Tau_exp fit plot ==========
tau_data = np.loadtxt(f'{resdir}/tau_exp_fit.dat')
lag = tau_data[:, 0]
rho_y2 = tau_data[:, 1]

# Read fit params from header
A_exp, tau_exp = 1.0, 5.0
with open(f'{resdir}/tau_exp_fit.dat') as f:
    for line in f:
        if 'fit_params' in line:
            parts = line.split()
            for p in parts:
                if p.startswith('A='):
                    A_exp = float(p[2:])
                elif p.startswith('tau='):
                    tau_exp = float(p[4:])
            break

fig1, ax1 = plt.subplots(figsize=(8, 5))
ax1.plot(lag[1:], rho_y2[1:], 'b-', lw=0.8, label='dati y2')
t_plot = np.linspace(1, lag[-1], 200)
ax1.plot(t_plot, A_exp * np.exp(-t_plot / tau_exp), 'r--', lw=1.5,
         label=f'fit $\\tau_{{exp}}$={tau_exp:.1f}')
ax1.set_xlabel('lag')
ax1.set_ylabel(r'$\rho(t)$')
ax1.set_title(f'Fit esponenziale y2 (nt={nt})')
ax1.legend(fontsize=9)
fig1.tight_layout()
fig1.savefig(f'{plotdir}/tau_exp_fit.png', dpi=150)

# ========== 2. Energy gaps ==========
gap_data = np.loadtxt(f'{resdir}/energy_gaps.dat')
n_gap = gap_data[:, 0].astype(int)
corr_labels = ['yc', 'y2c', 'y3c', 'Ac']
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

nmax = args.gap_nmax if args.gap_nmax > 0 else len(n_gap)
sel = slice(0, nmax)

fig2, ax2 = plt.subplots(figsize=(10, 6))
for i, (label, col) in enumerate(zip(corr_labels, colors)):
    de = gap_data[sel, 1 + 2*i]
    de_err = gap_data[sel, 2 + 2*i]
    ax2.errorbar(n_gap[sel], de, yerr=de_err, fmt='o', ms=3, color=col,
                 label=label, capsize=2, lw=0.8)
ax2.set_xlabel('n')
ax2.set_ylabel(r'$\Delta E(n)$')
ax2.set_ylim(0, 5)
for yline in (1.0, 2.0, 3.0):
    ax2.axhline(yline, color='k', ls='--', lw=0.5)
ax2.legend(fontsize=9)
fig2.suptitle(f'Energy gaps (nt={nt}, skip={skip})')
fig2.tight_layout()
fig2.savefig(f'{plotdir}/energy_gaps.png', dpi=150)

# ========== 3. Blocking plateau plots ==========
# 3a: simple observables
block_obs = np.loadtxt(f'{resdir}/blocking_observables.dat')
k_vals = block_obs[:, 0]
obs_names = ['y', 'y2', 'y3', 'A', 'E']

fig3, axes3 = plt.subplots(2, 3, figsize=(12, 7))
axes3 = axes3.flatten()
for i, name in enumerate(obs_names):
    ax = axes3[i]
    ax.plot(np.log2(k_vals), block_obs[:, 1 + i], 'o-', ms=3, lw=0.8)
    ax.set_xlabel(r'$\log_2(k)$')
    ax.set_ylabel(r'$\sigma(\langle ' + name + r' \rangle)$')
    ax.set_title(f'Blocking: {name}')
axes3[5].axis('off')
fig3.suptitle(f'Blocking observables (nt={nt}, skip={skip})')
fig3.tight_layout()
fig3.savefig(f'{plotdir}/blocking_observables.png', dpi=150)

# 3b: jackknife correlators
block_jack = np.loadtxt(f'{resdir}/blocking_jackknife_correlators.dat')
k_vals_j = block_jack[:, 0]

# Read n_check values from header
with open(f'{resdir}/blocking_jackknife_correlators.dat') as f:
    header_line = f.readline()
n_check_vals = []
for token in header_line.split():
    if token.startswith('sigma_') and '_n' in token:
        nv = int(token.split('_n')[-1])
        if nv not in n_check_vals:
            n_check_vals.append(nv)

corr_names = ['yc', 'y2c', 'y3c', 'Ac']
n_per_corr = len(n_check_vals)

fig4, axes4 = plt.subplots(2, 2, figsize=(10, 8))
axes4 = axes4.flatten()
for ic, cn in enumerate(corr_names):
    ax = axes4[ic]
    for jn, nv in enumerate(n_check_vals):
        col_idx = 1 + ic * n_per_corr + jn
        ax.plot(np.log2(k_vals_j), block_jack[:, col_idx], 'o-', ms=3, lw=0.8,
                label=f'n={nv}')
    ax.set_xlabel(r'$\log_2(k)$')
    ax.set_ylabel(r'$\sigma$')
    ax.set_title(f'Jackknife blocking: {cn}')
    ax.legend(fontsize=8)
fig4.suptitle(f'Jackknife blocking correlators (nt={nt}, skip={skip})')
fig4.tight_layout()
fig4.savefig(f'{plotdir}/blocking_jackknife_correlators.png', dpi=150)

plt.show()
print(f"Plots saved in {plotdir}/")
