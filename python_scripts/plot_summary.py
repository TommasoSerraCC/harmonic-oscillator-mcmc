import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument('--bhw', type=int, required=True)
parser.add_argument('--nstep', type=int, default=1000000)
parser.add_argument('--skip', type=int, default=None,
                    help='thermalization skip value (default: auto-detect all)')
parser.add_argument('-oy', action='store_true', help='plot <y> vs eta')
parser.add_argument('-o2', action='store_true', help='plot <y^2> vs eta')
parser.add_argument('-o3', action='store_true', help='plot <y^3> vs eta')
parser.add_argument('-oA', action='store_true', help='plot <A> vs eta')
parser.add_argument('-oE', action='store_true', help='plot <E> vs eta')
parser.add_argument('--save', action='store_true', help='save plots to file')
args = parser.parse_args()

# If no plot flag is given, show all
show_flags = [args.oy, args.o2, args.o3, args.oA, args.oE]
if not any(show_flags):
    show_flags = [True, True, True, True, True]

# Scan results directory for this bhw/nstep/ntherm
basedir = f'bhw{args.bhw}_nstep{args.nstep}'
resbase = f'results/{basedir}'
all_dirs = os.listdir(resbase) if os.path.isdir(resbase) else []

# Collect (nt, skip) pairs
data_map = {}  # skip -> list of (nt, filepath)
for d in all_dirs:
    m = re.match(r'nt(\d+)_therm(\d+)', d)
    if m:
        ntv = int(m.group(1))
        sk = int(m.group(2))
        fpath = f'{resbase}/{d}/observables.dat'
        if os.path.isfile(fpath):
            if sk not in data_map:
                data_map[sk] = []
            data_map[sk].append((ntv, fpath))

# Filter by skip if specified
if args.skip is not None:
    if args.skip in data_map:
        data_map = {args.skip: data_map[args.skip]}
    else:
        print(f"No data found for skip={args.skip}")
        exit(1)

if not data_map:
    print("No data found in results/")
    exit(1)

print(f"Found {len(data_map)} skip value(s): {sorted(data_map.keys())}")

obs_names = ['y', 'y2', 'y3', 'A', 'E']
obs_mean_labels = [r'$\langle y \rangle$', r'$\langle y^2 \rangle$',
                   r'$\langle y^3 \rangle$', r'$\langle A \rangle$',
                   r'$\langle E \rangle$']
os.makedirs(f'plots/{basedir}', exist_ok=True)

# Process each skip value separately
for skip_val in sorted(data_map.keys()):
    entries = sorted(data_map[skip_val])
    nt_arr = []
    means = {n: [] for n in obs_names}
    errs = {n: [] for n in obs_names}
    
    print(f"\nProcessing skip={skip_val}:")
    for ntv, fpath in entries:
        print(f"  nt={ntv}")
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
    eta_arr = float(args.bhw) / nt_arr
    
    # Plot each observable separately
    for i, (name, meanlabel, do_plot) in enumerate(
            zip(obs_names, obs_mean_labels, show_flags)):
        if not do_plot:
            continue
        m = np.array(means[name])
        e = np.array(errs[name])
        fig, ax = plt.subplots(figsize=(7, 5))
        if name in ('y2', 'E'):
            ax.errorbar(eta_arr**2, m, yerr=e, fmt='o', ms=4, capsize=0.8, elinewidth=0.6, markerfacecolor='none', markeredgewidth=0.6)
            ax.set_xlabel(r'$\eta^2$', fontsize=16)
        else:
            ax.errorbar(eta_arr, m, yerr=e, fmt='o', ms=4, capsize=0.8, elinewidth=0.6, markerfacecolor='none', markeredgewidth=0.6)
            ax.set_xlabel(r'$\eta$', fontsize=16)
        ax.set_ylabel(meanlabel, fontsize=16)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        if args.save:
            outfile = f'plots/{basedir}/{name}_vs_eta_therm{skip_val}.png'
            fig.savefig(outfile, dpi=150)
            print(f"  Saved {outfile}")

plt.show()
print("\nDone.")
