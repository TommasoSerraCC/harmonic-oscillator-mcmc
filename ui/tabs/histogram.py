"""
GUI tab for ground state histogram: runs MCMC simulation, plots position
distribution and compares with the theoretical ground state |psi_0(x)|^2.

This tab is independent of the rest of the analysis pipeline.
"""

import tkinter as tk
from tkinter import ttk, messagebox
import threading
import os
import numpy as np

from ui.core import data_manager as dm
from ui.core.plotting import PlotFrame, apply_grid, FIT_LINE_KW
from ui.core.ground_state_sim import run_ground_state_simulation


class HistogramTab(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent)

        self._cancel = False
        self._running = False

        # ── Controls ──────────────────────────────────────────────
        ctrl = tk.LabelFrame(self, text='Simulazione stato fondamentale')
        ctrl.pack(fill=tk.X, padx=8, pady=4)

        r1 = tk.Frame(ctrl)
        r1.pack(fill=tk.X, padx=4, pady=2)

        tk.Label(r1, text='βℏω:').pack(side=tk.LEFT, padx=2)
        self.bhw_var = tk.StringVar(value='10')
        tk.Entry(r1, textvariable=self.bhw_var, width=6).pack(side=tk.LEFT, padx=2)

        tk.Label(r1, text='N_t:').pack(side=tk.LEFT, padx=6)
        self.nt_var = tk.StringVar(value='50')
        tk.Entry(r1, textvariable=self.nt_var, width=6).pack(side=tk.LEFT, padx=2)

        tk.Label(r1, text='Configurazioni:').pack(side=tk.LEFT, padx=6)
        self.nstep_var = tk.StringVar(value='50000')
        tk.Entry(r1, textvariable=self.nstep_var, width=10).pack(side=tk.LEFT, padx=2)

        tk.Label(r1, text='Termalizzazione:').pack(side=tk.LEFT, padx=6)
        self.therm_var = tk.StringVar(value='10000')
        tk.Entry(r1, textvariable=self.therm_var, width=8).pack(side=tk.LEFT, padx=2)

        r2 = tk.Frame(ctrl)
        r2.pack(fill=tk.X, padx=4, pady=4)

        tk.Label(r2, text='Bins:').pack(side=tk.LEFT, padx=2)
        self.bins_var = tk.StringVar(value='200')
        tk.Entry(r2, textvariable=self.bins_var, width=5).pack(side=tk.LEFT, padx=2)

        self.run_btn = tk.Button(r2, text='Avvia simulazione',
                                 command=self._run)
        self.run_btn.pack(side=tk.LEFT, padx=12)

        self.cancel_btn = tk.Button(r2, text='Annulla', state='disabled',
                                    command=self._request_cancel)
        self.cancel_btn.pack(side=tk.LEFT, padx=4)

        tk.Button(r2, text='Salva plot', command=self._save).pack(
            side=tk.LEFT, padx=8)

        # ── Info label ───────────────────────────────────────────
        self.info_var = tk.StringVar()
        tk.Label(self, textvariable=self.info_var, fg='darkgreen',
                 anchor='w').pack(fill=tk.X, padx=8)

        # ── Status / progress ────────────────────────────────────
        self.status_var = tk.StringVar()
        tk.Label(self, textvariable=self.status_var, fg='blue',
                 anchor='w').pack(fill=tk.X, padx=8)

        self.progress = ttk.Progressbar(self, mode='determinate')
        self.progress.pack(fill=tk.X, padx=8, pady=2)

        # ── Plot ─────────────────────────────────────────────────
        self.plot = PlotFrame(self, figsize=(10, 6))
        self.plot.pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

        self._last_params = None

    # ── actions ───────────────────────────────────────────────────

    def _validate_params(self):
        try:
            bhw = float(self.bhw_var.get())
            nt = int(self.nt_var.get())
            nsteps = int(self.nstep_var.get())
            therm = int(self.therm_var.get())
            bins = int(self.bins_var.get())
            assert bhw > 0 and nt > 1 and nsteps > 0 and therm >= 0 and bins > 0
            return bhw, nt, nsteps, therm, bins
        except Exception:
            messagebox.showerror('Errore',
                                 'Parametri non validi. Controllare i campi.')
            return None

    def _run(self):
        params = self._validate_params()
        if params is None:
            return
        if self._running:
            return
        bhw, nt, nsteps, therm, bins = params
        self._cancel = False
        self._running = True
        self.run_btn.config(state='disabled')
        self.cancel_btn.config(state='normal')
        self.progress['value'] = 0
        self.status_var.set('Simulazione in corso…')

        threading.Thread(target=self._sim_thread,
                         args=(bhw, nt, nsteps, therm, bins),
                         daemon=True).start()

    def _request_cancel(self):
        self._cancel = True
        self.status_var.set('Annullamento in corso…')

    def _sim_thread(self, bhw, nt, nsteps, therm, bins):
        def progress_cb(step, total):
            pct = step / total * 100
            self.after(0, lambda p=pct: self.progress.configure(value=p))
            phase = 'termalizzazione' if step < therm else 'raccolta dati'
            self.after(0, lambda s=step, t=total, ph=phase:
                       self.status_var.set(
                           f'{ph}: {s}/{t} ({s/t*100:.0f}%)'))

        positions = run_ground_state_simulation(
            bhw, nt, nsteps, therm,
            progress_callback=progress_cb,
            cancel_flag=lambda: self._cancel)

        self.after(0, lambda: self._sim_done(positions, bhw, nt, nsteps,
                                             therm, bins))

    def _sim_done(self, positions, bhw, nt, nsteps, therm, bins):
        self._running = False
        self.run_btn.config(state='normal')
        self.cancel_btn.config(state='disabled')
        self.progress['value'] = 100

        if positions is None:
            self.status_var.set('Simulazione annullata.')
            return

        # Store params AND statistics
        mean_x = np.mean(positions)
        mean_x2 = np.mean(positions**2)
        mean_x4 = np.mean(positions**4)
        variance = np.var(positions)
        binder = mean_x4 / (mean_x2**2) if mean_x2 != 0 else 0
        
        self._last_params = {
            'bhw': bhw, 'nt': nt, 'nsteps': nsteps, 'therm': therm,
            'bins': bins, 'n_samples': len(positions),
            'mean': mean_x, 'variance': variance, 'binder': binder
        }

        self.status_var.set(
            f'Completato — {len(positions)} campioni '
            f'({nsteps} configurazioni × {nt} siti)')

        # Update info label with statistics
        self.info_var.set(
            f'μ = {self._last_params["mean"]:.4f}   '
            f'σ² = {self._last_params["variance"]:.4f}   '
            f'U = {self._last_params["binder"]:.4f}   '
            f'Campioni: {len(positions)}   '
            f'βℏω = {bhw}   N_t = {nt}')

        self._draw_histogram(positions, bhw, nt, nsteps, therm, bins)

    # ── plotting ──────────────────────────────────────────────────

    def _draw_histogram(self, positions, bhw, nt, nsteps, therm, bins):
        self.plot.set_layout(1, 1)
        ax = self.plot.axes[0]

        ax.hist(positions, bins=bins, density=True, alpha=0.7,
                color='tab:blue', edgecolor='black', linewidth=0.3,
                label=r'simulazione $\beta\hbar\omega={}, N_t={}$'.format(
                    bhw, nt))

        # Theoretical |psi_0(x)|^2 = exp(-x^2) / sqrt(pi)
        xmin, xmax = positions.min(), positions.max()
        x_th = np.linspace(xmin, xmax, 1000)
        psi2 = np.exp(-x_th ** 2) / np.sqrt(np.pi)
        ax.plot(x_th, psi2, **FIT_LINE_KW,
                label=r'$|\psi_0(x)|^2$ (teoria)')

        ax.set_xlabel(r'$y$')
        ax.set_ylabel(r'Densità di probabilità')
        ax.set_title('Istogramma stato fondamentale')
        ax.legend(loc='upper left', fontsize=9)

        apply_grid(ax)
        self.plot.draw()

    # ── save ──────────────────────────────────────────────────────

    def _save(self):
        if self._last_params is None:
            messagebox.showinfo('Info',
                                'Eseguire prima una simulazione.')
            return
        p = self._last_params
        save_dir = os.path.join(
            dm.PLOTS_DIR, 'istogrammi',
            f'bhw{p["bhw"]}_nt{p["nt"]}_nstep{p["nsteps"]}'
            f'_therm{p["therm"]}')
        os.makedirs(save_dir, exist_ok=True)
        
        # Save plot
        path = os.path.join(save_dir, 'ground_state_histogram.png')
        self.plot.save(path)
        
        # Save info file with statistics and parameters
        info_path = os.path.join(save_dir, 'histogram_info.txt')
        with open(info_path, 'w', encoding='utf-8') as f:
            f.write('=== Istogramma stato fondamentale ===\n\n')
            f.write('Parametri simulazione:\n')
            f.write(f'  βℏω = {p["bhw"]}\n')
            f.write(f'  N_t = {p["nt"]}\n')
            f.write(f'  Configurazioni = {p["nsteps"]}\n')
            f.write(f'  Termalizzazione = {p["therm"]}\n')
            f.write(f'  Bins = {p["bins"]}\n')
            f.write(f'\nStatistiche:\n')
            f.write(f'  Campioni totali = {p["n_samples"]}\n')
            f.write(f'  Media (μ) = {p["mean"]:.6f}\n')
            f.write(f'  Varianza (σ²) = {p["variance"]:.6f}\n')
            f.write(f'  Dev. standard (σ) = {np.sqrt(p["variance"]):.6f}\n')
            f.write(f'  Cumulante di Binder (U) = {p["binder"]:.6f}\n')
        
        messagebox.showinfo('Salvato', 
                           f'Plot salvato: {path}\nInfo salvate: {info_path}')
