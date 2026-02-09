import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser()
parser.add_argument('--nt', type=int, required=True)
parser.add_argument('--skip', type=int, default=0)
parser.add_argument('--bhw', type=int, required=True)
parser.add_argument('--nstep', type=int, default=1000000)
parser.add_argument('--nfit', type=int, default=0, help='max n for fit range (0=auto)')
parser.add_argument('--save', action='store_true')
args = parser.parse_args()

nt = args.nt
skip = args.skip
bhw = float(args.bhw)
eta = bhw / nt

basedir = f'bhw{args.bhw}_nstep{args.nstep}'
resdir = f'results/{basedir}/nt{nt}_therm{skip}'
plotdir = f'plots/{basedir}/nt{nt}_therm{skip}'
os.makedirs(plotdir, exist_ok=True)

data = np.loadtxt(f'{resdir}/connected_correlators.dat')
n = data[:, 0].astype(int)

corr_labels = [r'$C_y$', r'$C_{y^2}$', r'$C_{y^3}$', r'$C_A$']
corr_names = ['yc', 'y2c', 'y3c', 'Ac']

def exp_model(n, A, xi):
    return A * np.exp(-n / xi)

fig, axes = plt.subplots(2, 2, figsize=(10, 8))
axes = axes.flatten()

results = []

for i, (label, cname) in enumerate(zip(corr_labels, corr_names)):
    ax = axes[i]
    mean = data[:, 1 + 2*i]
    err = data[:, 2 + 2*i]

    # select positive values for fit
    pos = mean > 0
    n_pos = n[pos]
    m_pos = mean[pos]
    e_pos = err[pos]

    if args.nfit > 0:
        sel = n_pos <= args.nfit
        n_pos, m_pos, e_pos = n_pos[sel], m_pos[sel], e_pos[sel]

    # auto range: cut where signal/noise < 2
    sn = m_pos / e_pos
    good = sn > 2
    if good.sum() < 3:
        good = np.ones(len(n_pos), dtype=bool)
        good[min(10, len(good)):] = False
    n_fit = n_pos[good]
    m_fit = m_pos[good]
    e_fit = e_pos[good]

    try:
        popt, pcov = curve_fit(exp_model, n_fit, m_fit, p0=[m_fit[0], 5.0],
                               sigma=e_fit, absolute_sigma=True,
                               bounds=([0, 0.1], [np.inf, np.inf]))
        A_fit, xi_fit = popt
        A_err, xi_err = np.sqrt(np.diag(pcov))
        de = 1.0 / (xi_fit * eta)
        de_err = xi_err / (xi_fit**2 * eta)
    except:
        A_fit, xi_fit, A_err, xi_err = np.nan, np.nan, np.nan, np.nan
        de, de_err = np.nan, np.nan

    results.append((cname, xi_fit, xi_err, de, de_err))

    ax.errorbar(n, mean, yerr=err, fmt='o', ms=1.2, capsize=0.8, elinewidth=0.6, markerfacecolor='none', markeredgewidth=0.6, label='data')
    if not np.isnan(xi_fit):
        nplot = np.linspace(1, n_fit[-1] * 1.5, 200)
        ax.plot(nplot, exp_model(nplot, A_fit, xi_fit), 'r-', lw=1.0,
                label=f'$\\xi={xi_fit:.2f}\\pm{xi_err:.2f}$')
    ax.set_xlabel('$n$')
    ax.set_ylabel(label)
    ax.set_yscale('log')
    ax.legend(fontsize=8)

fig.suptitle(f'$N_t={nt}$, $\\eta={eta:.4f}$', fontsize=11)
fig.tight_layout()
if args.save:
    fig.savefig(f'{plotdir}/correlation_length_fit.png', dpi=150)
    print(f"Saved {plotdir}/correlation_length_fit.png")

# save fit results
with open(f'{resdir}/correlation_length.dat', 'w') as f:
    f.write("# corr  xi  xi_err  DeltaE  DeltaE_err\n")
    for cname, xi, xi_e, de, de_e in results:
        f.write(f"{cname}  {xi:.8e}  {xi_e:.8e}  {de:.8e}  {de_e:.8e}\n")
print(f"Correlation lengths saved in {resdir}/correlation_length.dat")
for cname, xi, xi_e, de, de_e in results:
    print(f"  {cname}: xi={xi:.3f}+/-{xi_e:.3f}, DeltaE={de:.4f}+/-{de_e:.4f}")

if not args.save:
    plt.show()
