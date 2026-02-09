import tkinter as tk
from tkinter import messagebox
import numpy as np
import os
from scipy.optimize import curve_fit

from ui.core import data_manager as dm
from ui.core.plotting import (PlotFrame, SavedFitsPanel,
                               apply_grid, plot_residuals, FIT_LINE_KW,
                               DATA_DOT_KW, format_value_with_uncertainty)


def _bose_model(inv_bhw, a, b):
    """E = a + b / (exp(1/inv_bhw) - 1)  =  a + b / (exp(bhw) - 1)

    Here inv_bhw = 1/(beta*hbar*omega) is the x-axis (temperature).
    bhw = 1/inv_bhw.
    """
    bhw = 1.0 / inv_bhw
    return a + b / (np.exp(bhw) - 1.0)


def _theory_energy(bhw):
    """Exact QHO energy <E> = 0.5 * coth(bhw/2) (in units of hbar*omega)."""
    return 0.5 / np.tanh(float(bhw) / 2.0)


class EnergyVsTempTab(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent)

        # --- Controls ---
        ctrl = tk.Frame(self)
        ctrl.pack(fill=tk.X, padx=8, pady=4)

        tk.Label(ctrl, text='\u03b7 (passo):').pack(side=tk.LEFT, padx=2)
        self.eta_var = tk.StringVar()
        self.eta_cb = tk.OptionMenu(ctrl, self.eta_var, '')
        self.eta_cb.config(width=10)
        self.eta_cb.pack(side=tk.LEFT, padx=2)

        tk.Label(ctrl, text='therm:').pack(side=tk.LEFT, padx=2)
        self.therm_var = tk.StringVar(value='')
        tk.Entry(ctrl, textvariable=self.therm_var, width=8).pack(side=tk.LEFT, padx=2)
        tk.Label(ctrl, text='(vuoto=qualsiasi)').pack(side=tk.LEFT, padx=2)

        tk.Button(ctrl, text='Plot', command=self._plot).pack(side=tk.LEFT, padx=8)
        tk.Button(ctrl, text='Aggiorna \u03b7', command=self._refresh_etas).pack(side=tk.LEFT, padx=4)

        r2 = tk.Frame(self)
        r2.pack(fill=tk.X, padx=8, pady=2)
        tk.Label(r2, text='Fit  a + b/(exp(\u03b2\u210f\u03c9)-1):').pack(side=tk.LEFT)
        tk.Label(r2, text='   range 1/\u03b2\u210f\u03c9 min:').pack(side=tk.LEFT)
        self.fit_xmin_var = tk.StringVar(value='')
        tk.Entry(r2, textvariable=self.fit_xmin_var, width=6).pack(side=tk.LEFT, padx=1)
        tk.Label(r2, text='max:').pack(side=tk.LEFT)
        self.fit_xmax_var = tk.StringVar(value='')
        tk.Entry(r2, textvariable=self.fit_xmax_var, width=6).pack(side=tk.LEFT, padx=1)
        tk.Button(r2, text='Fit', command=self._fit).pack(side=tk.LEFT, padx=8)
        tk.Button(r2, text='Salva fit', command=self._save_fit).pack(side=tk.LEFT, padx=4)

        self._show_res_var = tk.BooleanVar(value=False)
        self._res_btn = tk.Checkbutton(r2, text='Residui',
                                        variable=self._show_res_var,
                                        command=self._toggle_residuals)
        self._res_btn.pack(side=tk.LEFT, padx=4)
        self._res_btn.config(state='disabled')

        self.info_var = tk.StringVar()
        tk.Label(self, textvariable=self.info_var, fg='darkgreen', anchor='w',
                 wraplength=900, justify='left').pack(fill=tk.X, padx=8)

        bottom = tk.PanedWindow(self, orient=tk.HORIZONTAL)
        bottom.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        self.plot = PlotFrame(bottom, figsize=(9, 6))
        bottom.add(self.plot, stretch='always')

        self.fits_panel = SavedFitsPanel(bottom, self._get_fits_dir)
        bottom.add(self.fits_panel, width=260)

        self._fit_result = None
        self._last_data = None
        self._fit_draw_data = None
        self._refresh_etas()

    def _get_fits_dir(self):
        return dm.energy_vs_temp_dir()

    def _refresh_etas(self):
        etas = dm.get_available_etas_across_bhw()
        menu = self.eta_cb['menu']
        menu.delete(0, 'end')
        for e in etas:
            label = f'{e:.6g}'
            menu.add_command(label=label, command=tk._setit(self.eta_var, label))
        if etas:
            self.eta_var.set(f'{etas[0]:.6g}')
        self.fits_panel.refresh()

    def _gather_data(self):
        eta_s = self.eta_var.get().strip()
        if not eta_s:
            return None
        eta = float(eta_s)
        therm_s = self.therm_var.get().strip()
        therm = int(therm_s) if therm_s else None
        results = dm.collect_energy_vs_temperature(eta, therm)
        if not results:
            return None
        inv_bhw = np.array([r['inv_bhw'] for r in results])
        E = np.array([r['E'] for r in results])
        E_err = np.array([r['E_err'] for r in results])
        bhw_vals = np.array([r['bhw'] for r in results])
        self._last_data = {
            'eta': eta, 'inv_bhw': inv_bhw, 'E': E, 'E_err': E_err,
            'bhw_vals': bhw_vals, 'results': results
        }
        return self._last_data

    def _plot(self):
        d = self._gather_data()
        if d is None:
            messagebox.showinfo('Info', 'Nessun dato trovato per questa \u03b7')
            return

        self.plot.set_layout(1, 1)
        ax = self.plot.axes[0]

        ax.errorbar(d['inv_bhw'], d['E'], yerr=d['E_err'], **DATA_DOT_KW)

        # Theory curve
        x_th = np.linspace(max(d['inv_bhw'].min() * 0.5, 0.001),
                           d['inv_bhw'].max() * 1.2, 200)
        bhw_th = 1.0 / x_th
        E_th = np.array([_theory_energy(b) for b in bhw_th])
        ax.plot(x_th, E_th, 'k--', lw=0.8, label='teoria esatta')

        ax.set_xlabel(r'$1/(\beta\hbar\omega)$')
        ax.set_ylabel(r'$\langle E \rangle$')
        ax.set_title(f'$\\eta = {d["eta"]:.4g}$', fontsize=11)
        ax.legend(fontsize=8)
        apply_grid(ax)
        self.plot.draw()
        self._fit_result = None
        self.info_var.set(f'{len(d["inv_bhw"])} punti dati')
        self.fits_panel.refresh()

    def _fit(self):
        d = self._gather_data()
        if d is None:
            messagebox.showinfo('Info', 'Nessun dato')
            return

        xmin = float(self.fit_xmin_var.get()) if self.fit_xmin_var.get().strip() else 0
        xmax = float(self.fit_xmax_var.get()) if self.fit_xmax_var.get().strip() else 9999

        sel = (d['inv_bhw'] >= xmin) & (d['inv_bhw'] <= xmax)
        x = d['inv_bhw'][sel]
        y = d['E'][sel]
        yerr = d['E_err'][sel]

        if len(x) < 2:
            messagebox.showwarning('Fit', 'Troppo pochi punti')
            return

        try:
            popt, pcov = curve_fit(_bose_model, x, y, p0=[0.5, 1.0],
                                   sigma=yerr, absolute_sigma=True)
            a_fit, b_fit = popt
            perr = np.sqrt(np.diag(pcov))
            corr_matrix = pcov / np.outer(perr, perr)
            residuals = (y - _bose_model(x, *popt)) / yerr
            chi2 = np.sum(residuals**2)
            ndof = len(x) - 2
            chi2red = chi2 / ndof if ndof > 0 else np.nan
        except Exception as e:
            messagebox.showerror('Fit', str(e))
            return

        self._fit_result = {
            'eta': d['eta'],
            'a': float(a_fit), 'a_err': float(perr[0]),
            'b': float(b_fit), 'b_err': float(perr[1]),
            'chi2red': float(chi2red), 'ndof': int(ndof),
            'corr_ab': float(corr_matrix[0, 1]),
            'xmin': float(xmin), 'xmax': float(xmax),
            'n_points': int(len(x))
        }

        self._fit_draw_data = {
            'eta': d['eta'],
            'inv_bhw_all': d['inv_bhw'], 'E_all': d['E'], 'E_err_all': d['E_err'],
            'x_sel': x, 'y_sel': y, 'yerr_sel': yerr,
            'popt': popt, 'residuals': residuals,
        }

        self._show_res_var.set(False)
        self._res_btn.config(state='normal')
        self._draw_fit_full()

        v_a, e_a = format_value_with_uncertainty(a_fit, perr[0])
        v_b, e_b = format_value_with_uncertainty(b_fit, perr[1])
        self.info_var.set(
            f'a = {v_a} \u00b1 {e_a}   '
            f'b = {v_b} \u00b1 {e_b}   '
            f'\u03c7\u00b2/ndof = {chi2red:.3f} ({ndof})   '
            f'corr(a,b) = {corr_matrix[0,1]:.4f}   '
            f'cov = [[{pcov[0,0]:.2e}, {pcov[0,1]:.2e}], [{pcov[1,0]:.2e}, {pcov[1,1]:.2e}]]'
        )

    def _draw_fit_full(self):
        """Draw full data with fit line (no residuals)."""
        d = self._fit_draw_data
        self.plot.set_layout(1, 1)
        ax = self.plot.axes[0]
        ax.errorbar(d['inv_bhw_all'], d['E_all'], yerr=d['E_err_all'],
                     **DATA_DOT_KW, label='dati')
        x_plot = np.linspace(max(d['inv_bhw_all'].min() * 0.5, 0.001),
                             d['inv_bhw_all'].max() * 1.2, 200)
        ax.plot(x_plot, _bose_model(x_plot, *d['popt']), **FIT_LINE_KW,
                label='fit')
        bhw_th = 1.0 / x_plot
        E_th = np.array([_theory_energy(b) for b in bhw_th])
        ax.plot(x_plot, E_th, 'k--', lw=0.8, label='teoria esatta')
        ax.set_xlabel(r'$1/(\beta\hbar\omega)$')
        ax.set_ylabel(r'$\langle E \rangle$')
        ax.set_title(f'$\\eta = {d["eta"]:.4g}$', fontsize=11)
        ax.legend(fontsize=8)
        apply_grid(ax)
        self.plot.draw()

    def _draw_fit_with_residuals(self):
        """Draw data + residuals, both trimmed to fit range."""
        d = self._fit_draw_data
        ax_main, ax_res = self.plot.set_layout_with_residuals(height_ratios=(3, 1))
        x, y, yerr = d['x_sel'], d['y_sel'], d['yerr_sel']
        ax_main.errorbar(x, y, yerr=yerr, **DATA_DOT_KW, label='dati')
        x_plot = np.linspace(x.min(), x.max(), 200)
        ax_main.plot(x_plot, _bose_model(x_plot, *d['popt']), **FIT_LINE_KW,
                     label='fit')
        bhw_th = 1.0 / x_plot
        E_th = np.array([_theory_energy(b) for b in bhw_th])
        ax_main.plot(x_plot, E_th, 'k--', lw=0.8, label='teoria esatta')
        ax_main.set_ylabel(r'$\langle E \rangle$')
        ax_main.set_title(f'$\\eta = {d["eta"]:.4g}$', fontsize=11)
        ax_main.legend(fontsize=8)
        apply_grid(ax_main)
        margin = (x.max() - x.min()) * 0.05 if len(x) > 1 else 0.01
        xlim = (x.min() - margin, x.max() + margin)
        ax_main.set_xlim(xlim)
        plot_residuals(ax_res, x, d['residuals'])
        ax_res.set_xlabel(r'$1/(\beta\hbar\omega)$')
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
        base = dm.energy_vs_temp_dir()
        r = self._fit_result
        tag = dm.save_fit_entry(base, f'E_vs_T_eta{r["eta"]:.4g}', r, self.plot.fig)
        messagebox.showinfo('Salvato', f'Fit salvato: {tag}')
        self.fits_panel.refresh()
