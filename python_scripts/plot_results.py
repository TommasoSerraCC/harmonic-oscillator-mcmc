import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--nt', type=int, required=True)
parser.add_argument('--skip', type=int, default=0, help='thermalization steps')
parser.add_argument('--bhw', type=int, required=True)
parser.add_argument('--nstep', type=int, default=1000000)
parser.add_argument('-te', action='store_true', help='plot tau_exp fit')
parser.add_argument('-eg', action='store_true', help='plot energy gaps')
parser.add_argument('-bo', action='store_true', help='plot blocking observables')
parser.add_argument('-jc', action='store_true', help='plot jackknife correlators')
parser.add_argument('--nmax', type=int, default=0, help='max n for energy gap plot (0=auto)')
parser.add_argument('--save', action='store_true', help='save plots to file')
args = parser.parse_args()

# If no plot flag is given, show all
if not (args.te or args.eg or args.bo or args.jc):
    args.te = args.eg = args.bo = args.jc = True

nt = args.nt
skip = args.skip
ncorr = nt // 2
bhw = float(args.bhw)
eta = bhw / nt

basedir = f'bhw{args.bhw}_nstep{args.nstep}'
resdir = f'results/{basedir}/nt{nt}_therm{skip}'
plotdir = f'plots/{basedir}/nt{nt}_therm{skip}'
os.makedirs(plotdir, exist_ok=True)

# ========== 1. Tau_exp fit plot ==========
if args.te:
    tau_data = np.loadtxt(f'{resdir}/tau_exp_fit.dat')
    lag = tau_data[:, 0]
    acf_y2 = tau_data[:, 1]

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
    ax1.plot(lag[1:], acf_y2[1:], 'o', color='black', ms=1.2, markerfacecolor='none', markeredgewidth=0.6, label=r'dati $C_{y^2}$')
    t_plot = np.linspace(1, lag[-1], 200)
    ax1.plot(t_plot, A_exp * np.exp(-t_plot / tau_exp), 'r-', lw=1.0,
             label=f'fit $\\tau_{{exp}}$={tau_exp:.1f}')
    ax1.set_xlabel('$n$', fontsize=16)
    ax1.set_ylabel(r'$C_{y^2}(n)$', fontsize=16)
    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax1.legend(fontsize=11)
    fig1.tight_layout()
    if args.save:
        fig1.savefig(f'{plotdir}/tau_exp_fit.png', dpi=150)
        print(f"Saved {plotdir}/tau_exp_fit.png")

# ========== 2. Energy gaps ==========
if args.eg:
    gap_data = np.loadtxt(f'{resdir}/energy_gaps.dat', ndmin=2)
    n_gap = gap_data[:, 0].astype(int)
    corr_labels = [r'$y$', r'$y^2$', r'$y^3$', r'$A$']
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

    nmax = args.nmax if args.nmax > 0 else len(n_gap)
    sel = slice(0, nmax)

    fig2, ax2 = plt.subplots(figsize=(10, 6))
    for i, (label, col) in enumerate(zip(corr_labels, colors)):
        de = gap_data[sel, 1 + 2*i]
        de_err = gap_data[sel, 2 + 2*i]
        ax2.errorbar(n_gap[sel] * eta, de, yerr=de_err, fmt='o', ms=1.2, color=col,
                     label=label, capsize=0.8, elinewidth=0.6, markerfacecolor='none', markeredgewidth=0.6)
    ax2.set_xlabel(r'$n\eta$', fontsize=16)
    ax2.set_ylabel(r'$\Delta E(n)$', fontsize=16)
    ax2.tick_params(axis='both', which='major', labelsize=14)
    ax2.set_ylim(0, 5)
    for yline in (1.0, 2.0, 3.0):
        ax2.axhline(yline, color='k', ls='--', lw=0.5)
    ax2.legend(fontsize=11)
    fig2.tight_layout()
    if args.save:
        fig2.savefig(f'{plotdir}/energy_gaps.png', dpi=150)
        print(f"Saved {plotdir}/energy_gaps.png")

# ========== 3. Blocking observables ==========
if args.bo:
    block_obs = np.loadtxt(f'{resdir}/blocking_observables.dat')
    k_vals = block_obs[:, 0]
    obs_labels = [r'$y$', r'$y^2$', r'$y^3$', r'$A$', r'$E$']
    obs_sigma_labels = [r'$\sigma(\langle y \rangle)$', r'$\sigma(\langle y^2 \rangle)$',
                        r'$\sigma(\langle y^3 \rangle)$', r'$\sigma(\langle A \rangle)$',
                        r'$\sigma(\langle E \rangle)$']

    fig3, axes3 = plt.subplots(2, 3, figsize=(12, 7))
    axes3 = axes3.flatten()
    for i, (label, slabel) in enumerate(zip(obs_labels, obs_sigma_labels)):
        ax = axes3[i]
        ax.plot(np.log2(k_vals), block_obs[:, 1 + i], '-o', ms=1.2, markerfacecolor='none', markeredgewidth=0.6)
        ax.set_xlabel(r'$\log_2(k)$', fontsize=14)
        ax.set_ylabel(slabel, fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=12)
    axes3[5].axis('off')
    fig3.tight_layout()
    if args.save:
        fig3.savefig(f'{plotdir}/blocking_observables.png', dpi=150)
        print(f"Saved {plotdir}/blocking_observables.png")

# ========== 4. Jackknife blocking correlators ==========
if args.jc:
    block_jack = np.loadtxt(f'{resdir}/blocking_jackknife_correlators.dat')
    k_vals_j = block_jack[:, 0]

    with open(f'{resdir}/blocking_jackknife_correlators.dat') as f:
        header_line = f.readline()
    n_check_vals = []
    for token in header_line.split():
        if token.startswith('sigma_') and '_n' in token:
            nv = int(token.split('_n')[-1])
            if nv not in n_check_vals:
                n_check_vals.append(nv)

    corr_labels = [r'$C_y$', r'$C_{y^2}$', r'$C_{y^3}$', r'$C_A$']
    corr_sigma_labels = [r'$\sigma(C_y)$', r'$\sigma(C_{y^2})$', r'$\sigma(C_{y^3})$', r'$\sigma(C_A)$']
    n_per_corr = len(n_check_vals)

    fig4, axes4 = plt.subplots(2, 2, figsize=(10, 8))
    axes4 = axes4.flatten()
    for ic, (clabel, slabel) in enumerate(zip(corr_labels, corr_sigma_labels)):
        ax = axes4[ic]
        for jn, nv in enumerate(n_check_vals):
            col_idx = 1 + ic * n_per_corr + jn
            ax.plot(np.log2(k_vals_j), block_jack[:, col_idx], '-o', ms=1.2, markerfacecolor='none', markeredgewidth=0.6,
                    label=f'n={nv}')
        ax.set_xlabel(r'$\log_2(k)$', fontsize=14)
        ax.set_ylabel(slabel, fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.legend(fontsize=10)
    fig4.tight_layout()
    if args.save:
        fig4.savefig(f'{plotdir}/blocking_jackknife_correlators.png', dpi=150)
        print(f"Saved {plotdir}/blocking_jackknife_correlators.png")

plt.show()
