import tkinter as tk
from tkinter import messagebox
import numpy as np
import os

from ui.core import data_manager as dm
from ui.core.plotting import (PlotFrame, ParamSelector, SavedFitsPanel,
                               apply_grid, plot_residuals, FIT_LINE_KW,
                               SMALL_MS, format_value_with_uncertainty)


CORR_NAMES = ['yc', 'y2c', 'y3c', 'Ac']
CORR_DISPLAY = {'yc': r'$y$', 'y2c': r'$y^2$', 'y3c': r'$y^3$', 'Ac': r'$A$'}
COLORS = {'yc': 'tab:blue', 'y2c': 'tab:orange', 'y3c': 'tab:green', 'Ac': 'tab:red'}
MARKERS_FOR_CORR = {'yc': 'o', 'y2c': 's', 'y3c': '^', 'Ac': 'D'}


class EnergyGapsTab(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent)

        ctrl = tk.Frame(self)
        ctrl.pack(fill=tk.X, padx=8, pady=4)
        self.selector = ParamSelector(ctrl, dm)
        self.selector.pack(side=tk.LEFT)

        tk.Label(ctrl, text='  Operatore:').pack(side=tk.LEFT, padx=4)
        self.corr_var = tk.StringVar(value='yc')
        for cn in CORR_NAMES:
            tk.Radiobutton(ctrl, text=cn, variable=self.corr_var, value=cn).pack(side=tk.LEFT)

        self.all_var = tk.BooleanVar(value=False)
        tk.Checkbutton(ctrl, text='Tutti', variable=self.all_var).pack(side=tk.LEFT, padx=4)

        r2 = tk.Frame(self)
        r2.pack(fill=tk.X, padx=8, pady=2)

        tk.Label(r2, text='n\u00b7\u03b7 range:').pack(side=tk.LEFT, padx=2)
        self.neta_min_var = tk.StringVar(value='')
        self.neta_max_var = tk.StringVar(value='')
        tk.Entry(r2, textvariable=self.neta_min_var, width=6).pack(side=tk.LEFT, padx=1)
        tk.Label(r2, text='-').pack(side=tk.LEFT)
        tk.Entry(r2, textvariable=self.neta_max_var, width=6).pack(side=tk.LEFT, padx=1)

        tk.Button(r2, text='Plot', command=self._plot).pack(side=tk.LEFT, padx=8)

        tk.Label(r2, text='  Fit costante n:').pack(side=tk.LEFT, padx=4)
        self.fit_nmin_var = tk.StringVar(value='')
        self.fit_nmax_var = tk.StringVar(value='')
        tk.Entry(r2, textvariable=self.fit_nmin_var, width=5).pack(side=tk.LEFT, padx=1)
        tk.Label(r2, text='-').pack(side=tk.LEFT)
        tk.Entry(r2, textvariable=self.fit_nmax_var, width=5).pack(side=tk.LEFT, padx=1)

        tk.Button(r2, text='Fit', command=self._fit).pack(side=tk.LEFT, padx=8)
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
        bottom.pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

        self.plot = PlotFrame(bottom)
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

    def _get_data(self):
        p = self.selector.get_params()
        if not p or not dm.results_exist(p['bhw'], p['nstep'], p['nt'], p['therm']):
            return None, None
        data = dm.load_energy_gaps(p['bhw'], p['nstep'], p['nt'], p['therm'])
        return p, data

    def _plot(self):
        p, data = self._get_data()
        if data is None:
            messagebox.showinfo('Info', 'Risultati non disponibili')
            return
        eta = float(p['bhw']) / p['nt']
        n = data['n']
        neta = n * eta

        neta_min = float(self.neta_min_var.get()) if self.neta_min_var.get() else 0
        neta_max = float(self.neta_max_var.get()) if self.neta_max_var.get() else neta[-1] * 1.1

        self.plot.set_layout(1, 1)
        ax = self.plot.axes[0]

        if self.all_var.get():
            plot_list = CORR_NAMES
        else:
            plot_list = [self.corr_var.get()]

        for cn in plot_list:
            de = data[cn]
            de_err = data[cn + '_err']
            sel = (neta >= neta_min) & (neta <= neta_max)
            marker = MARKERS_FOR_CORR[cn]
            ax.errorbar(neta[sel], de[sel], yerr=de_err[sel],
                        fmt=marker, ms=SMALL_MS, color=COLORS[cn],
                        label=CORR_DISPLAY[cn], capsize=2, elinewidth=0.8)

        ax.set_xlabel(r'$n\eta$')
        ax.set_ylabel(r'$\Delta E(n)$')
        ax.set_ylim(0, 5)
        for yline in (1.0, 2.0, 3.0):
            ax.axhline(yline, color='k', ls='--', lw=0.5)
        if len(plot_list) > 1 or self.all_var.get():
            ax.legend(fontsize=9)
        apply_grid(ax)
        self.plot.draw()
        self._fit_result = None
        self.info_var.set('')
        self.fits_panel.refresh()

    def _fit(self):
        p, data = self._get_data()
        if data is None:
            return
        cn = self.corr_var.get()
        eta = float(p['bhw']) / p['nt']
        n = data['n']
        de = data[cn]
        de_err = data[cn + '_err']

        nmin_s = self.fit_nmin_var.get().strip()
        nmax_s = self.fit_nmax_var.get().strip()
        if not nmin_s or not nmax_s:
            messagebox.showwarning('Fit', 'Specificare range n per fit costante')
            return
        nmin, nmax = int(nmin_s), int(nmax_s)

        sel = (n >= nmin) & (n <= nmax)
        valid = sel & np.isfinite(de) & np.isfinite(de_err) & (de_err > 0)
        nf = n[valid]
        df = de[valid]
        ef = de_err[valid]

        if len(nf) < 1:
            messagebox.showwarning('Fit', 'Nessun punto valido nel range')
            return

        w = 1.0 / ef**2
        wm = np.sum(w * df) / np.sum(w)
        wm_err = 1.0 / np.sqrt(np.sum(w))
        residuals = (df - wm) / ef
        chi2 = np.sum(residuals**2)
        ndof = len(nf) - 1
        chi2red = chi2 / ndof if ndof > 0 else np.nan

        self._fit_result = {
            'corr': cn, 'value': float(wm), 'error': float(wm_err),
            'chi2red': float(chi2red), 'ndof': int(ndof),
            'n_min': int(nmin), 'n_max': int(nmax),
            'bhw': p['bhw'], 'nstep': p['nstep'],
            'nt': p['nt'], 'therm': p['therm']
        }

        neta = n * eta
        neta_min = float(self.neta_min_var.get()) if self.neta_min_var.get() else 0
        neta_max = float(self.neta_max_var.get()) if self.neta_max_var.get() else neta[-1] * 1.1

        self._fit_draw_data = {
            'cn': cn, 'eta': eta, 'p': p,
            'n_all': n, 'neta_all': neta,
            'de_all': de, 'de_err_all': de_err,
            'neta_min': neta_min, 'neta_max': neta_max,
            'nf': nf, 'df': df, 'ef': ef,
            'wm': wm, 'wm_err': wm_err,
            'residuals': residuals,
        }

        self._show_res_var.set(False)
        self._res_btn.config(state='normal')
        self._draw_fit_full()

        v_de, e_de = format_value_with_uncertainty(wm, wm_err)
        self.info_var.set(
            f'\u0394E = {v_de} \u00b1 {e_de}   '
            f'\u03c7\u00b2/ndof = {chi2red:.3f} ({ndof})   '
            f'range: n \u2208 [{nmin}, {nmax}]'
        )

    def _draw_fit_full(self):
        """Draw full data with fit line overlaid (no residuals)."""
        d = self._fit_draw_data
        cn = d['cn']
        self.plot.set_layout(1, 1)
        ax = self.plot.axes[0]

        neta = d['neta_all']
        de = d['de_all']
        de_err = d['de_err_all']
        plot_sel = (neta >= d['neta_min']) & (neta <= d['neta_max'])
        marker = MARKERS_FOR_CORR[cn]
        ax.errorbar(neta[plot_sel], de[plot_sel], yerr=de_err[plot_sel],
                    fmt=marker, ms=SMALL_MS, color=COLORS[cn],
                    label=CORR_DISPLAY[cn], capsize=2, elinewidth=0.8)
        ax.axhline(d['wm'], **FIT_LINE_KW, label='fit')
        ax.axhspan(d['wm'] - d['wm_err'], d['wm'] + d['wm_err'],
                   alpha=0.15, color='r')
        ax.set_xlabel(r'$n\eta$')
        ax.set_ylabel(r'$\Delta E(n)$')
        ax.set_ylim(0, 5)
        for yline in (1.0, 2.0, 3.0):
            ax.axhline(yline, color='k', ls='--', lw=0.5)
        ax.legend(fontsize=9)
        apply_grid(ax)
        self.plot.draw()

    def _draw_fit_with_residuals(self):
        """Draw fit data + residuals, both trimmed to fit range."""
        d = self._fit_draw_data
        cn = d['cn']
        ax_main, ax_res = self.plot.set_layout_with_residuals(height_ratios=(3, 1))

        nf_eta = d['nf'] * d['eta']
        marker = MARKERS_FOR_CORR[cn]
        ax_main.errorbar(nf_eta, d['df'], yerr=d['ef'],
                         fmt=marker, ms=SMALL_MS, color=COLORS[cn],
                         label=CORR_DISPLAY[cn], capsize=2, elinewidth=0.8)
        ax_main.axhline(d['wm'], **FIT_LINE_KW, label='fit')
        ax_main.axhspan(d['wm'] - d['wm_err'], d['wm'] + d['wm_err'],
                        alpha=0.15, color='r')
        ax_main.set_ylabel(r'$\Delta E(n)$')
        ax_main.set_ylim(0, 5)
        for yline in (1.0, 2.0, 3.0):
            ax_main.axhline(yline, color='k', ls='--', lw=0.5)
        ax_main.legend(fontsize=9)
        apply_grid(ax_main)
        margin = d['eta'] * 0.5
        xlim = (nf_eta[0] - margin, nf_eta[-1] + margin)
        ax_main.set_xlim(xlim)

        plot_residuals(ax_res, nf_eta, d['residuals'])
        ax_res.set_xlabel(r'$n\eta$')
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
        tag = dm.save_fit_entry(base, f'energy_gap_{r["corr"]}', r, self.plot.fig)
        messagebox.showinfo('Salvato', f'Fit salvato: {tag}')
        self.fits_panel.refresh()
