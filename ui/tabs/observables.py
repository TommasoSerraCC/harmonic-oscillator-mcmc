import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import os
from scipy.optimize import curve_fit

from ui.core import data_manager as dm
from ui.core.plotting import (PlotFrame, ParamSelector, SavedFitsPanel,
                               apply_grid, plot_residuals, FIT_LINE_KW,
                               DATA_DOT_KW, format_value_with_uncertainty)


def _theory_value(bhw):
    x = np.exp(-float(bhw))
    return 0.5 * (1 + x) / (1 - x)


class ObservablesTab(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        self.summary_frame = SummarySubtab(self.notebook)
        self.blocking_frame = BlockingSubtab(self.notebook)
        self.autocorr_frame = AutocorrSubtab(self.notebook)

        self.notebook.add(self.summary_frame, text='Medie vs \u03b7')
        self.notebook.add(self.blocking_frame, text='Errore statistico')
        self.notebook.add(self.autocorr_frame, text='Autocorrelazione')


# =====================================================================
#  Summary (means vs eta) with continuum-limit fit
# =====================================================================

class SummarySubtab(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent)

        ctrl = tk.Frame(self)
        ctrl.pack(fill=tk.X, padx=8, pady=4)

        tk.Label(ctrl, text='\u03b2\u210f\u03c9:').pack(side=tk.LEFT, padx=2)
        self.bhw_var = tk.StringVar()
        self.bhw_cb = tk.OptionMenu(ctrl, self.bhw_var, '')
        self.bhw_cb.config(width=6)
        self.bhw_cb.pack(side=tk.LEFT, padx=2)

        tk.Label(ctrl, text='therm:').pack(side=tk.LEFT, padx=2)
        self.therm_var = tk.StringVar()
        self.therm_cb = tk.OptionMenu(ctrl, self.therm_var, '')
        self.therm_cb.config(width=8)
        self.therm_cb.pack(side=tk.LEFT, padx=2)

        self.obs_var = tk.StringVar(value='y2')
        tk.Label(ctrl, text='Osservabile:').pack(side=tk.LEFT, padx=6)
        for name in ['y', 'y2', 'y3', 'A', 'E']:
            tk.Radiobutton(ctrl, text=name, variable=self.obs_var, value=name).pack(side=tk.LEFT)

        tk.Button(ctrl, text='Plot', command=self._plot).pack(side=tk.LEFT, padx=8)

        r2 = tk.Frame(self)
        r2.pack(fill=tk.X, padx=8, pady=2)
        tk.Label(r2, text='Fit lineare \u03b7\u00b2 (E, y\u00b2):  range Nt min:').pack(side=tk.LEFT)
        self.fit_ntmin_var = tk.StringVar(value='')
        tk.Entry(r2, textvariable=self.fit_ntmin_var, width=5).pack(side=tk.LEFT, padx=2)
        tk.Label(r2, text='max:').pack(side=tk.LEFT)
        self.fit_ntmax_var = tk.StringVar(value='')
        tk.Entry(r2, textvariable=self.fit_ntmax_var, width=5).pack(side=tk.LEFT, padx=2)
        tk.Button(r2, text='Fit', command=self._fit_continuum).pack(side=tk.LEFT, padx=8)
        tk.Button(r2, text='Salva fit', command=self._save_fit).pack(side=tk.LEFT, padx=4)

        self._show_res_var = tk.BooleanVar(value=False)
        self._res_btn = tk.Checkbutton(r2, text='Residui',
                                        variable=self._show_res_var,
                                        command=self._toggle_residuals)
        self._res_btn.pack(side=tk.LEFT, padx=4)
        self._res_btn.config(state='disabled')

        self.info_var = tk.StringVar()
        tk.Label(self, textvariable=self.info_var, fg='darkgreen', anchor='w').pack(fill=tk.X, padx=8)

        bottom = tk.PanedWindow(self, orient=tk.HORIZONTAL)
        bottom.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        self.plot = PlotFrame(bottom, figsize=(9, 6))
        bottom.add(self.plot, stretch='always')

        self.fits_panel = SavedFitsPanel(bottom, self._get_fits_dir)
        bottom.add(self.fits_panel, width=260)

        self.bhw_var.trace_add('write', lambda *_: self._refresh_therms())
        self._refresh_bhw()
        self._fit_result = None
        self._last_data = None
        self._fit_draw_data = None

    def _get_fits_dir(self):
        bhw_s = self.bhw_var.get()
        if not bhw_s:
            return None
        nstep = self.nstep_map.get(bhw_s)
        if not nstep:
            return None
        return os.path.join(dm.RESULTS_DIR, f'bhw{bhw_s}_nstep{nstep}')

    def _refresh_bhw(self):
        sets = dm.scan_data_sets()
        self.nstep_map = {str(b): n for b, n in sets}
        menu = self.bhw_cb['menu']
        menu.delete(0, 'end')
        for bhw, _ in sets:
            menu.add_command(label=str(bhw), command=tk._setit(self.bhw_var, str(bhw)))
        if sets:
            self.bhw_var.set(str(sets[0][0]))

    def _refresh_therms(self):
        bhw_s = self.bhw_var.get()
        if not bhw_s:
            return
        nstep = self.nstep_map.get(bhw_s)
        if not nstep:
            return
        therms = dm.get_all_available_therms(int(bhw_s), nstep)
        menu = self.therm_cb['menu']
        menu.delete(0, 'end')
        for t in therms:
            menu.add_command(label=str(t), command=tk._setit(self.therm_var, str(t)))
        if therms:
            self.therm_var.set(str(therms[0]))
        self.fits_panel.refresh()

    def _gather_data(self):
        bhw_s = self.bhw_var.get()
        if not bhw_s:
            return None
        bhw = int(bhw_s)
        nstep = self.nstep_map.get(bhw_s)
        if not nstep:
            return None
        therm_s = self.therm_var.get()
        therm = int(therm_s) if therm_s else None
        entries = dm.collect_observables_vs_eta(bhw, nstep, therm)
        if not entries:
            return None
        obs = self.obs_var.get()
        nt_arr, m_arr, e_arr = [], [], []
        for nt, th, obs_dict in entries:
            if obs in obs_dict:
                nt_arr.append(nt)
                m_arr.append(obs_dict[obs][0])
                e_arr.append(obs_dict[obs][1])
        if not nt_arr:
            return None
        nt_arr = np.array(nt_arr)
        eta_arr = float(bhw) / nt_arr
        self._last_data = {
            'bhw': bhw, 'nstep': nstep, 'obs': obs,
            'nt': nt_arr, 'eta': eta_arr,
            'mean': np.array(m_arr), 'err': np.array(e_arr)
        }
        return self._last_data

    def _plot(self):
        d = self._gather_data()
        if d is None:
            messagebox.showinfo('Info', 'Nessun dato trovato')
            return
        obs = d['obs']
        self.plot.set_layout(1, 1)
        ax = self.plot.axes[0]
        bhw = d['bhw']
        labels = {'y': r'$\langle y \rangle$', 'y2': r'$\langle y^2 \rangle$',
                  'y3': r'$\langle y^3 \rangle$', 'A': r'$\langle A \rangle$',
                  'E': r'$\langle E \rangle$'}
        if obs in ('y2', 'E'):
            ax.errorbar(d['eta']**2, d['mean'], yerr=d['err'],
                        **DATA_DOT_KW)
            ax.set_xlabel(r'$\eta^2$')
            tv = _theory_value(bhw)
            ax.axhline(tv, color='k', ls='--', lw=0.8,
                       label=f'teoria = {tv:.6f}')
            ax.legend(fontsize=8)
        else:
            ax.errorbar(d['eta'], d['mean'], yerr=d['err'],
                        **DATA_DOT_KW)
            ax.set_xlabel(r'$\eta$')
            ax.axhline(0, color='k', ls='--', lw=0.8)
        ax.set_ylabel(labels.get(obs, obs))
        apply_grid(ax)
        self.plot.draw()
        self._fit_result = None
        self.info_var.set('')

    def _fit_continuum(self):
        d = self._gather_data()
        if d is None:
            return
        obs = d['obs']
        if obs not in ('y2', 'E'):
            messagebox.showinfo('Info', 'Fit disponibile solo per E e y\u00b2')
            return

        ntmin = int(self.fit_ntmin_var.get()) if self.fit_ntmin_var.get().strip() else 0
        ntmax = int(self.fit_ntmax_var.get()) if self.fit_ntmax_var.get().strip() else 9999
        sel = (d['nt'] >= ntmin) & (d['nt'] <= ntmax)
        eta2 = d['eta'][sel]**2
        m = d['mean'][sel]
        e = d['err'][sel]
        if len(eta2) < 2:
            messagebox.showwarning('Fit', 'Troppo pochi punti')
            return

        def linear(x, a, b):
            return a + b * x

        popt, pcov = curve_fit(linear, eta2, m, sigma=e, absolute_sigma=True)
        a_fit, b_fit = popt
        perr = np.sqrt(np.diag(pcov))
        corr_matrix = pcov / np.outer(perr, perr)
        residuals = (m - linear(eta2, *popt)) / e
        chi2 = np.sum(residuals**2)
        ndof = len(eta2) - 2
        chi2red = chi2 / ndof if ndof > 0 else np.nan

        bhw = d['bhw']
        tv = _theory_value(bhw)
        self._fit_result = {
            'obs': obs, 'a': float(a_fit), 'a_err': float(perr[0]),
            'b': float(b_fit), 'b_err': float(perr[1]),
            'chi2red': float(chi2red), 'ndof': int(ndof),
            'corr_ab': float(corr_matrix[0, 1]),
            'theory': float(tv), 'bhw': bhw, 'nstep': d['nstep'],
            'therm': self.therm_var.get(),
            'ntmin': ntmin, 'ntmax': ntmax
        }

        self._fit_draw_data = {
            'obs': obs, 'popt': popt,
            'eta_all': d['eta'], 'mean_all': d['mean'], 'err_all': d['err'],
            'eta2_sel': eta2, 'm_sel': m, 'e_sel': e,
            'residuals': residuals, 'tv': tv,
            'linear': linear,
        }

        self._show_res_var.set(False)
        self._res_btn.config(state='normal')
        self._draw_fit_full()

        v_a, e_a = format_value_with_uncertainty(a_fit, perr[0])
        v_b, e_b = format_value_with_uncertainty(b_fit, perr[1])
        self.info_var.set(
            f'a = {v_a} \u00b1 {e_a}   b = {v_b} \u00b1 {e_b}   '
            f'\u03c7\u00b2/ndof = {chi2red:.3f} ({ndof})   '
            f'corr(a,b) = {corr_matrix[0,1]:.4f}   '
            f'teoria = {tv:.6f}'
        )

    def _draw_fit_full(self):
        """Draw full data with fit line (no residuals)."""
        d = self._fit_draw_data
        obs = d['obs']
        labels = {'y2': r'$\langle y^2 \rangle$', 'E': r'$\langle E \rangle$'}
        self.plot.set_layout(1, 1)
        ax = self.plot.axes[0]
        ax.errorbar(d['eta_all']**2, d['mean_all'], yerr=d['err_all'],
                     **DATA_DOT_KW, label='dati')
        x_plot = np.linspace(0, max(d['eta_all']**2) * 1.1, 200)
        ax.plot(x_plot, d['linear'](x_plot, *d['popt']), **FIT_LINE_KW,
                label='fit')
        ax.axhline(d['tv'], color='k', ls='--', lw=0.8,
                   label=f'teoria = {d["tv"]:.6f}')
        ax.set_xlabel(r'$\eta^2$')
        ax.set_ylabel(labels.get(obs, obs))
        ax.legend(fontsize=8)
        apply_grid(ax)
        self.plot.draw()

    def _draw_fit_with_residuals(self):
        """Draw data + residuals, both trimmed to fit range."""
        d = self._fit_draw_data
        obs = d['obs']
        labels = {'y2': r'$\langle y^2 \rangle$', 'E': r'$\langle E \rangle$'}
        ax_main, ax_res = self.plot.set_layout_with_residuals(height_ratios=(3, 1))
        eta2 = d['eta2_sel']
        ax_main.errorbar(eta2, d['m_sel'], yerr=d['e_sel'],
                         **DATA_DOT_KW, label='dati')
        x_plot = np.linspace(eta2.min(), eta2.max(), 200)
        ax_main.plot(x_plot, d['linear'](x_plot, *d['popt']), **FIT_LINE_KW,
                     label='fit')
        ax_main.axhline(d['tv'], color='k', ls='--', lw=0.8,
                        label=f'teoria = {d["tv"]:.6f}')
        ax_main.set_ylabel(labels.get(obs, obs))
        ax_main.legend(fontsize=8)
        apply_grid(ax_main)
        margin = (eta2.max() - eta2.min()) * 0.05 if len(eta2) > 1 else 0.01
        xlim = (eta2.min() - margin, eta2.max() + margin)
        ax_main.set_xlim(xlim)
        plot_residuals(ax_res, eta2, d['residuals'])
        ax_res.set_xlabel(r'$\eta^2$')
        ax_res.set_xlim(xlim)
        self.plot.draw()

    def _toggle_residuals(self):
        if not self._fit_draw_data:
            return
        if self._show_res_var.get():
            self._draw_fit_with_residuals()
        else:
            self._draw_fit_full()

    def _save_fit(self):
        if not self._fit_result:
            messagebox.showinfo('Info', 'Eseguire prima un fit')
            return
        r = self._fit_result
        base = os.path.join(dm.RESULTS_DIR, f'bhw{r["bhw"]}_nstep{r["nstep"]}')
        tag = dm.save_fit_entry(base, f'continuum_{r["obs"]}', r, self.plot.fig)
        messagebox.showinfo('Salvato', f'Fit salvato: {tag}')
        self.fits_panel.refresh()


# =====================================================================
#  Blocking / Jackknife error analysis
# =====================================================================

class BlockingSubtab(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent)

        ctrl = tk.Frame(self)
        ctrl.pack(fill=tk.X, padx=8, pady=4)

        self.selector = ParamSelector(ctrl, dm)
        self.selector.pack(side=tk.LEFT)

        self.mode_var = tk.StringVar(value='simple')
        tk.Label(ctrl, text='  Tipo:').pack(side=tk.LEFT, padx=2)
        tk.Radiobutton(ctrl, text='Osservabili primarie (blocking)',
                       variable=self.mode_var, value='simple').pack(side=tk.LEFT)
        tk.Radiobutton(ctrl, text='Correlatori connessi (jackknife a blocchi)',
                       variable=self.mode_var, value='jack').pack(side=tk.LEFT)

        r2 = tk.Frame(self)
        r2.pack(fill=tk.X, padx=8, pady=2)
        tk.Label(r2, text='k_max (log\u2082):').pack(side=tk.LEFT, padx=2)
        self.kmax_var = tk.StringVar(value='')
        tk.Entry(r2, textvariable=self.kmax_var, width=5).pack(side=tk.LEFT, padx=2)
        tk.Label(r2, text='  (vuoto = tutti)').pack(side=tk.LEFT)

        tk.Button(r2, text='Plot', command=self._plot).pack(side=tk.LEFT, padx=12)
        tk.Button(r2, text='Salva', command=self._save).pack(side=tk.LEFT, padx=2)

        self.plot = PlotFrame(self, figsize=(10, 6))
        self.plot.pack(fill=tk.BOTH, expand=True)

    def _plot(self):
        p = self.selector.get_params()
        if not p:
            return
        if not dm.results_exist(**{k: p[k] for k in ('bhw', 'nstep', 'nt', 'therm')}):
            messagebox.showinfo('Info', 'Risultati non disponibili')
            return
        kmax_s = self.kmax_var.get().strip()
        if self.mode_var.get() == 'simple':
            self._plot_simple(p, kmax_s)
        else:
            self._plot_jack(p, kmax_s)

    def _plot_simple(self, p, kmax_s):
        bd = dm.load_blocking_observables(p['bhw'], p['nstep'], p['nt'], p['therm'])
        k = bd['k']
        log2k = np.log2(k)
        sel = log2k <= float(kmax_s) if kmax_s else np.ones(len(k), dtype=bool)

        obs_names = ['y', 'y2', 'y3', 'A', 'E']
        labels = [r'$\sigma(\langle y \rangle)$', r'$\sigma(\langle y^2 \rangle)$',
                  r'$\sigma(\langle y^3 \rangle)$', r'$\sigma(\langle A \rangle)$',
                  r'$\sigma(\langle E \rangle)$']

        self.plot.set_layout(2, 3)
        for i, (name, lab) in enumerate(zip(obs_names, labels)):
            ax = self.plot.axes[i]
            ax.plot(log2k[sel], bd[name][sel], '-o', ms=2)
            ax.set_xlabel(r'$\log_2(k)$')
            ax.set_ylabel(lab)
            apply_grid(ax)
        if len(self.plot.axes) > 5:
            self.plot.axes[5].axis('off')
        self.plot.draw()

    def _plot_jack(self, p, kmax_s):
        bd = dm.load_blocking_jackknife(p['bhw'], p['nstep'], p['nt'], p['therm'])
        k = bd['k']
        log2k = np.log2(k)
        sel = log2k <= float(kmax_s) if kmax_s else np.ones(len(k), dtype=bool)

        corr_names = ['yc', 'y2c', 'y3c', 'Ac']
        labels = [r'$\sigma(C_y)$', r'$\sigma(C_{y^2})$',
                  r'$\sigma(C_{y^3})$', r'$\sigma(C_A)$']

        self.plot.set_layout(2, 2)
        for ic, (cn, lab) in enumerate(zip(corr_names, labels)):
            ax = self.plot.axes[ic]
            for nv in bd['n_check']:
                ax.plot(log2k[sel], bd[cn][nv][sel], '-o', ms=2, label=f'n={nv}')
            ax.set_xlabel(r'$\log_2(k)$')
            ax.set_ylabel(lab)
            ax.legend(fontsize=8)
            apply_grid(ax)
        self.plot.draw()

    def _save(self):
        p = self.selector.get_params()
        if not p:
            return
        pd = dm._plotdir(p['bhw'], p['nstep'], p['nt'], p['therm'])
        mode = self.mode_var.get()
        fname = f'blocking_{"observables" if mode == "simple" else "jackknife"}.png'
        path = os.path.join(pd, fname)
        self.plot.save(path)
        messagebox.showinfo('Salvato', path)


# =====================================================================
#  Autocorrelation (tau_exp) with variable-window fit
# =====================================================================

class AutocorrSubtab(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent)

        ctrl = tk.Frame(self)
        ctrl.pack(fill=tk.X, padx=8, pady=4)
        self.selector = ParamSelector(ctrl, dm)
        self.selector.pack(side=tk.LEFT)
        tk.Button(ctrl, text='Plot', command=self._plot).pack(side=tk.LEFT, padx=8)

        r2 = tk.Frame(self)
        r2.pack(fill=tk.X, padx=8, pady=2)
        tk.Label(r2, text='Fit range n:').pack(side=tk.LEFT, padx=4)
        self.nmin_var = tk.StringVar(value='1')
        self.nmax_var = tk.StringVar(value='')
        tk.Entry(r2, textvariable=self.nmin_var, width=5).pack(side=tk.LEFT, padx=1)
        tk.Label(r2, text='-').pack(side=tk.LEFT)
        tk.Entry(r2, textvariable=self.nmax_var, width=5).pack(side=tk.LEFT, padx=1)
        tk.Label(r2, text='(vuoto=auto)').pack(side=tk.LEFT, padx=2)
        tk.Button(r2, text='Fit exp', command=self._fit).pack(side=tk.LEFT, padx=8)
        tk.Button(r2, text='Salva fit', command=self._save_fit).pack(side=tk.LEFT, padx=4)

        self._show_res_var = tk.BooleanVar(value=False)
        self._res_btn = tk.Checkbutton(r2, text='Residui',
                                        variable=self._show_res_var,
                                        command=self._toggle_residuals)
        self._res_btn.pack(side=tk.LEFT, padx=4)
        self._res_btn.config(state='disabled')

        self.info_var = tk.StringVar()
        tk.Label(self, textvariable=self.info_var, fg='darkgreen', anchor='w').pack(fill=tk.X, padx=8)

        bottom = tk.PanedWindow(self, orient=tk.HORIZONTAL)
        bottom.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        self.plot = PlotFrame(bottom, figsize=(8, 5))
        bottom.add(self.plot, stretch='always')

        self.fits_panel = SavedFitsPanel(bottom, self._get_fits_dir)
        bottom.add(self.fits_panel, width=260)

        self._fit_result = None
        self._fit_draw_data = None

    def _get_fits_dir(self):
        p = self.selector.get_params()
        if not p:
            return None
        return dm._resdir(p['bhw'], p['nstep'], p['nt'], p['therm'])

    def _plot(self):
        p = self.selector.get_params()
        if not p or not dm.results_exist(p['bhw'], p['nstep'], p['nt'], p['therm']):
            messagebox.showinfo('Info', 'Risultati non disponibili')
            return
        tau_data = dm.load_tau_exp_fit(p['bhw'], p['nstep'], p['nt'], p['therm'])
        tau_int = dm.load_tau_int(p['bhw'], p['nstep'], p['nt'], p['therm'])

        self.plot.set_layout(1, 1)
        ax = self.plot.axes[0]
        lag = tau_data['lag']
        acf = tau_data['acf']
        ax.plot(lag[1:], acf[1:], 'o', color='black', ms=1.5,
                label=r'dati $C_{y^2}$')
        t_plot = np.linspace(1, lag[-1], 200)
        A, tau = tau_data['A'], tau_data['tau']
        ax.plot(t_plot, A * np.exp(-t_plot / tau), **FIT_LINE_KW,
                label=f'fit $\\tau_{{exp}}$={tau:.1f}')
        ax.set_xlabel('$n$')
        ax.set_ylabel(r'$C_{y^2}(n)$')
        ax.legend(fontsize=9)
        apply_grid(ax)
        self.plot.draw()

        info_parts = [f'\u03c4_exp={tau:.2f}']
        for name in ['y', 'y2', 'y3', 'A', 'E']:
            if name in tau_int:
                info_parts.append(f'\u03c4_int({name})={tau_int[name]:.2f}')
        self.info_var.set('   '.join(info_parts))
        self._fit_result = None
        self.fits_panel.refresh()

    def _fit(self):
        p = self.selector.get_params()
        if not p or not dm.results_exist(p['bhw'], p['nstep'], p['nt'], p['therm']):
            messagebox.showinfo('Info', 'Risultati non disponibili')
            return
        tau_data = dm.load_tau_exp_fit(p['bhw'], p['nstep'], p['nt'], p['therm'])
        lag = tau_data['lag']
        acf = tau_data['acf']

        nmin = int(self.nmin_var.get()) if self.nmin_var.get().strip() else 1
        nmax_s = self.nmax_var.get().strip()

        # Select positive data points
        pos = (acf > 0) & (lag >= nmin)
        if nmax_s:
            pos &= (lag <= int(nmax_s))
        else:
            # Auto: use points where acf has dropped by at least a factor of 5
            pos &= acf > acf[1] * 0.01 if len(acf) > 1 else pos
        lf = lag[pos]
        af = acf[pos]

        if len(lf) < 2:
            messagebox.showwarning('Fit', 'Troppo pochi punti per il fit')
            return

        # Estimate errors as proportional to values (simple approach)
        ef = np.abs(af) * 0.1 + 1e-10

        def exp_model(n, A, tau):
            return A * np.exp(-n / tau)

        try:
            popt, pcov = curve_fit(exp_model, lf, af, p0=[af[0], 5.0],
                                   sigma=ef, absolute_sigma=True,
                                   bounds=([0, 0.1], [np.inf, np.inf]))
            A_fit, tau_fit = popt
            A_err, tau_err = np.sqrt(np.diag(pcov))
            corr_mat = pcov / np.outer(np.sqrt(np.diag(pcov)), np.sqrt(np.diag(pcov)))
            residuals = (af - exp_model(lf, *popt)) / ef
            chi2 = np.sum(residuals**2)
            ndof = len(lf) - 2
            chi2red = chi2 / ndof if ndof > 0 else np.nan
        except Exception as e:
            messagebox.showerror('Fit', str(e))
            return

        self._fit_result = {
            'A': float(A_fit), 'A_err': float(A_err),
            'tau': float(tau_fit), 'tau_err': float(tau_err),
            'chi2red': float(chi2red), 'ndof': int(ndof),
            'corr_A_tau': float(corr_mat[0, 1]),
            'n_min': int(lf[0]), 'n_max': int(lf[-1]),
            'bhw': p['bhw'], 'nstep': p['nstep'],
            'nt': p['nt'], 'therm': p['therm']
        }

        self._fit_draw_data = {
            'lag': lag, 'acf': acf,
            'lf': lf, 'af': af, 'ef': ef,
            'popt': popt, 'residuals': residuals,
        }

        self._show_res_var.set(False)
        self._res_btn.config(state='normal')
        self._draw_fit_full()

        v_A, e_A = format_value_with_uncertainty(A_fit, A_err)
        v_tau, e_tau = format_value_with_uncertainty(tau_fit, tau_err)
        self.info_var.set(
            f'A = {v_A} \u00b1 {e_A}   '
            f'\u03c4_exp = {v_tau} \u00b1 {e_tau}   '
            f'\u03c7\u00b2/ndof = {chi2red:.3f} ({ndof})   '
            f'corr(A,\u03c4) = {corr_mat[0,1]:.4f}   '
            f'range: n \u2208 [{lf[0]}, {lf[-1]}]'
        )

    def _draw_fit_full(self):
        """Draw full data with fit line (no residuals)."""
        d = self._fit_draw_data
        def _exp(n, A, tau):
            return A * np.exp(-n / tau)
        self.plot.set_layout(1, 1)
        ax = self.plot.axes[0]
        ax.plot(d['lag'][1:], d['acf'][1:], 'o', color='black', ms=1.5, label='dati')
        t_plot = np.linspace(d['lf'][0], d['lf'][-1] * 1.3, 200)
        ax.plot(t_plot, _exp(t_plot, *d['popt']), **FIT_LINE_KW, label='fit')
        ax.set_xlabel('$n$')
        ax.set_ylabel(r'$C_{y^2}(n)$')
        ax.legend(fontsize=9)
        apply_grid(ax)
        self.plot.draw()

    def _draw_fit_with_residuals(self):
        """Draw fit data + residuals, both trimmed to fit range."""
        d = self._fit_draw_data
        def _exp(n, A, tau):
            return A * np.exp(-n / tau)
        ax_main, ax_res = self.plot.set_layout_with_residuals(height_ratios=(3, 1))
        lf = d['lf']
        ax_main.plot(lf, d['af'], 'o', color='black', ms=1.5, label='dati')
        t_plot = np.linspace(lf[0], lf[-1], 200)
        ax_main.plot(t_plot, _exp(t_plot, *d['popt']), **FIT_LINE_KW, label='fit')
        ax_main.set_ylabel(r'$C_{y^2}(n)$')
        ax_main.legend(fontsize=9)
        apply_grid(ax_main)
        xlim = (lf[0] - 0.5, lf[-1] + 0.5)
        ax_main.set_xlim(xlim)
        plot_residuals(ax_res, lf, d['residuals'])
        ax_res.set_xlabel('$n$')
        ax_res.set_xlim(xlim)
        self.plot.draw()

    def _toggle_residuals(self):
        if not self._fit_draw_data:
            return
        if self._show_res_var.get():
            self._draw_fit_with_residuals()
        else:
            self._draw_fit_full()

    def _save_fit(self):
        if not self._fit_result:
            messagebox.showinfo('Info', 'Eseguire prima un fit')
            return
        r = self._fit_result
        base = dm._resdir(r['bhw'], r['nstep'], r['nt'], r['therm'])
        tag = dm.save_fit_entry(base, 'tau_exp', r, self.plot.fig)
        messagebox.showinfo('Salvato', f'Fit salvato: {tag}')
        self.fits_panel.refresh()
