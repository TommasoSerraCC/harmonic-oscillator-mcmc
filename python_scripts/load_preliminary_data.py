import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--nt', type=int, default=50, choices=[50, 100, 200])
args = parser.parse_args()

nt = args.nt
ncorr = nt // 2
ncols_expected = 5 + 4 * ncorr

data = np.loadtxt(f'preliminary_data/raw_data_nt{nt}.dat')

# Verifica formato dati
assert data.shape[0] == 100000, f"Righe attese: 100000, trovate: {data.shape[0]}"
assert data.shape[1] == ncols_expected, f"Colonne attese: {ncols_expected}, trovate: {data.shape[1]}"

# Estrai osservabili
y = data[:, 0]
y2 = data[:, 1]
y3 = data[:, 2]
A = data[:, 3]
E = data[:, 4]

# Estrai correlatori (shape: nsteps x ncorr)
corr_data = data[:, 5:]
yc = corr_data[:, 0::4]
y2c = corr_data[:, 1::4]
y3c = corr_data[:, 2::4]
Ac = corr_data[:, 3::4]

assert yc.shape[1] == ncorr, f"ncorr atteso: {ncorr}, trovato: {yc.shape[1]}"

# Calcolo medie osservabili
obs_names = ['y', 'y2', 'y3', 'A', 'E']
obs_means = [y.mean(), y2.mean(), y3.mean(), A.mean(), E.mean()]
print(f"nt={nt}, ncorr={ncorr}, shape={data.shape}")
for name, mean in zip(obs_names, obs_means):
    print(f"  <{name}> = {mean:.6f}")

# Calcolo medie correlatori
corr_means = {
    'yc': yc.mean(axis=0) - y.mean()**2,
    'y2c': y2c.mean(axis=0) - y2.mean()**2,
    'y3c': y3c.mean(axis=0) - y3.mean()**2,
    'Ac': Ac.mean(axis=0) - A.mean()**2
}

# Plot 1: andamento temporale delle 5 osservabili (primi 1000 step)
fig1, axes1 = plt.subplots(5, 1, figsize=(10, 10), sharex=True)
obs_data = [y, y2, y3, A, E]
for ax, name, d in zip(axes1, obs_names, obs_data):
    ax.plot(d[:1000], lw=0.5)
    ax.set_ylabel(name)
    ax.axhline(d.mean(), color='r', ls='--', lw=1)
axes1[-1].set_xlabel('step MCMC')
fig1.suptitle(f'Osservabili (nt={nt})')
fig1.tight_layout()
fig1.savefig(f'plots/observables_nt{nt}.png', dpi=150)

# Plot 2: massa efficace m_eff(n) = log(C(n)/C(n+1))
bhw = 5.0
eta = bhw / nt
n_vals = np.arange(1, ncorr)  # n da 1 a ncorr-1

# Plot 2: tutti e quattro i gap sullo stesso grafico
fig2, ax2 = plt.subplots(figsize=(10, 6))
corr_labels = ['yc', 'y2c', 'y3c', 'Ac']
corr_arrays = [corr_means['yc'], corr_means['y2c'], corr_means['y3c'], corr_means['Ac']]
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

for C, label, col in zip(corr_arrays, corr_labels, colors):
    # m_eff(n) = log(C(n)/C(n+1)) / eta
    m_eff = np.log(C[:-1] / C[1:]) / eta
    ax2.plot(n_vals, m_eff, marker='o', ms=3, lw=1, color=col, label=label)

# Etichette asse
ax2.set_xlabel('n')
ax2.set_ylabel(r'$\Delta E(n)$')

# Mostra solo i dati con valori di Delta E tra 0 e 5
ax2.set_ylim(0, 5)

# Rimuovo box/legend visibile (nessuna didascalia a riquadro)
if ax2.get_legend() is not None:
    ax2.get_legend().remove()

# Linee orizzontali sottili, tratteggiate, nere in corrispondenza di 1,2,3
for yline in (1.0, 2.0, 3.0):
    ax2.axhline(yline, color='k', ls='--', lw=0.5)

fig2.suptitle(f'Energy gaps (nt={nt})')
fig2.tight_layout()
fig2.savefig(f'plots/meff_nt{nt}.png', dpi=150)

# === AUTOCORRELAZIONE ===
def autocorr(x, max_lag):
    """Calcola autocorrelazione normalizzata fino a max_lag."""
    n = len(x)
    x = x - x.mean()
    var = np.sum(x**2)
    rho = np.array([np.sum(x[:n-t] * x[t:]) / var for t in range(max_lag)])
    return rho

def tau_int(rho):
    """Calcola tau_int usando formula (F): tau_int = sum_{k=1}^{inf} C_F(k)"""
    # Somma fino a quando rho diventa trascurabile o negativo
    # Formula: tau_int = sum_{k=1}^{inf} rho(k)
    # In pratica, sommiamo fino a un cutoff ragionevole
    tau = 0.0
    for k in range(1, int(len(rho)/2)):
        tau += rho[k]
    return tau

max_lag = 500

print("\n=== Autocorrelazione ===")
obs_data = [y, y2, y3, A, E]
tau_ints = []
rhos = []

for name, d in zip(obs_names, obs_data):
    rho = autocorr(d, max_lag)
    rhos.append(rho)
    ti = tau_int(rho)
    tau_ints.append(ti)
    print(f"  tau_int({name}) = {ti:.2f}")

# Fit esponenziale per tau_exp (uso y2 come riferimento)
from scipy.optimize import curve_fit
def exp_decay(t, A, tau):
    return A * np.exp(-t / tau)

# Trova primo zero-crossing per limitare il fit
rho_y2 = rhos[1]
idx_zero = np.argmax(rho_y2[1:] < 0.01) + 1
idx_zero = max(idx_zero, 10)  # almeno 10 punti
t_fit = np.arange(1, idx_zero)
try:
    popt, _ = curve_fit(exp_decay, t_fit, rho_y2[1:idx_zero], p0=[1, 5], bounds=([0, 0.1], [np.inf, 100]))
    A_exp, tau_exp = popt
except:
    A_exp = 1.0
    tau_exp = tau_ints[1]  # fallback
print(f"  tau_exp (fit y2) = {tau_exp:.2f}")

# Plot 3: funzioni di autocorrelazione
fig3, axes3 = plt.subplots(2, 3, figsize=(12, 7))
axes3 = axes3.flatten()

for i, (name, rho) in enumerate(zip(obs_names, rhos)):
    ax = axes3[i]
    ax.plot(rho[:100], lw=0.8)
    ax.axvline(tau_ints[i], color='r', ls='--', label=f'tau_int={tau_ints[i]:.1f}')
    ax.set_xlabel('lag')
    ax.set_ylabel(r'$\rho(t)$')
    ax.set_title(name)
#    ax.set_yscale('log')
    ax.legend(fontsize=8)

# Ultimo pannello: fit esponenziale su y2
ax = axes3[5]
ax.plot(np.arange(1, 51), rho_y2[1:51], 'b-', lw=0.8, label='dati y2')
t_plot = np.linspace(1, 50, 100)
ax.plot(t_plot, exp_decay(t_plot, A_exp, tau_exp), 'r--', lw=1.5, label=f'fit tau_exp={tau_exp:.1f}')
ax.set_xlabel('lag')
ax.set_ylabel(r'$\rho(t)$')
ax.set_title('Fit esponenziale (y2)')
#ax.set_yscale('log')
ax.legend(fontsize=8)

fig3.suptitle(f'Autocorrelazione (nt={nt})')
fig3.tight_layout()
fig3.savefig(f'plots/autocorr_nt{nt}.png', dpi=150)

plt.show()
