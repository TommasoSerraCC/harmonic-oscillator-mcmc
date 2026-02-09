import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--nt', type=int, required=True)
parser.add_argument('--skip', type=int, default=0)
parser.add_argument('--bhw', type=int, required=True)
parser.add_argument('--nstep', type=int, default=1000000)
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

corr_labels = [r'$C_y(n)$', r'$C_{y^2}(n)$', r'$C_{y^3}(n)$', r'$C_A(n)$']
corr_names = ['yc', 'y2c', 'y3c', 'Ac']

fig, axes = plt.subplots(2, 2, figsize=(10, 8))
axes = axes.flatten()
for i, (label, cname) in enumerate(zip(corr_labels, corr_names)):
    ax = axes[i]
    mean = data[:, 1 + 2*i]
    err = data[:, 2 + 2*i]
    ax.errorbar(n, mean, yerr=err, fmt='o', ms=3, capsize=2, elinewidth=1.0)
    ax.set_xlabel('$n$')
    ax.set_ylabel(label)
    ax.set_yscale('log')
fig.suptitle(f'$N_t={nt}$, $\\eta={eta:.4f}$', fontsize=11)
fig.tight_layout()
if args.save:
    fig.savefig(f'{plotdir}/correlators.png', dpi=150)
    print(f"Saved {plotdir}/correlators.png")

if not args.save:
    plt.show()
