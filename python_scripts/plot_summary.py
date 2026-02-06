import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument('--skip', type=int, nargs='+', required=True,
                    help='thermalization skip values to use (one per nt, or one for all)')
parser.add_argument('--nt', type=int, nargs='+', default=None,
                    help='specific nt values (default: auto-detect from results/)')
args = parser.parse_args()

# Find available nt directories
all_dirs = os.listdir('results') if os.path.isdir('results') else []

if args.nt:
    nt_vals = sorted(args.nt)
else:
    nt_vals = []
    for d in all_dirs:
        m = re.match(r'nt(\d+)_therm(\d+)', d)
        if m:
            nt_vals.append(int(m.group(1)))
    nt_vals = sorted(set(nt_vals))

# Expand skip list
if len(args.skip) == 1:
    skip_vals = [args.skip[0]] * len(nt_vals)
else:
    skip_vals = args.skip

print(f"nt values: {nt_vals}")
print(f"skip values: {skip_vals}")

obs_names = ['y', 'y2', 'y3', 'A', 'E']

# Collect data
nt_arr = []
means = {n: [] for n in obs_names}
errs = {n: [] for n in obs_names}

for ntv, sk in zip(nt_vals, skip_vals):
    fpath = f'results/nt{ntv}_therm{sk}/observables.dat'
    if not os.path.isfile(fpath):
        print(f"  missing: {fpath}")
        continue
    with open(fpath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            name = parts[0]
            if name in obs_names:
                means[name].append(float(parts[1]))
                errs[name].append(float(parts[2]))
    nt_arr.append(ntv)

nt_arr = np.array(nt_arr)
eta_arr = 10.0 / nt_arr

# Plot 5 observables vs eta (or nt)
fig, axes = plt.subplots(3, 2, figsize=(11, 10))
axes = axes.flatten()

for i, name in enumerate(obs_names):
    ax = axes[i]
    m = np.array(means[name])
    e = np.array(errs[name])
    ax.errorbar(eta_arr, m, yerr=e, fmt='o-', ms=4, capsize=3, lw=0.8)
    ax.set_xlabel(r'$\eta = \beta\hbar\omega / N_t$')
    ax.set_ylabel(f'$\\langle {name} \\rangle$')
    ax.set_title(f'<{name}> vs $\\eta$')
    ax.grid(True, alpha=0.3)

axes[5].axis('off')
fig.suptitle(r'Observables vs $\eta$ ($\beta\hbar\omega=10$)')
fig.tight_layout()

os.makedirs('plots', exist_ok=True)
fig.savefig('plots/observables_vs_eta.png', dpi=150)
plt.show()
print("Saved plots/observables_vs_eta.png")
