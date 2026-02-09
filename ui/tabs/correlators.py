import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import os
from scipy.optimize import curve_fit

from ui.core import data_manager as dm
from ui.core.plotting import (PlotFrame, ParamSelector, SavedFitsPanel,
                               apply_grid, plot_residuals, FIT_LINE_KW,
                               DATA_DOT_KW, SMALL_MARKERS, SMALL_MS,
                               format_value_with_uncertainty)


class CorrelatorsTab(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        self.single = SingleCorrSubtab(self.notebook)
        self.multi = MultiNtSubtab(self.notebook)

        self.notebook.add(self.single, text='Singolo correlatore')
        self.notebook.add(self.multi, text='Confronto multi-Nt')


CORR_NAMES = ['yc', 'y2c', 'y3c', 'Ac']
CORR_LABELS = {'yc': r'$C_y(n)$', 'y2c': r'$C_{y^2}(n)$',
               'y3c': r'$C_{y^3}(n)$', 'Ac': r'$C_A(n)$'}


def exp_model(n, A, xi):
    return A * np.exp(-n / xi)


class SingleCorrSubtab(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent)

        ctrl = tk.Frame(self)
        ctrl.pack(fill=tk.X, padx=8, pady=4)
        self.selector = ParamSelector(ctrl, dm)
        self.selector.pack(side=tk.LEFT)

        tk.Label(ctrl, text='  Correlatore:').pack(side=tk.LEFT, padx=4)
        self.corr_var = tk.StringVar(value='yc')
        for cn in CORR_NAMES:
            tk.Radiobutton(ctrl, text=cn, variable=self.corr_var, value=cn).pack(side=tk.LEFT)

        btns = tk.Frame(self)
        btns.pack(fill=tk.X, padx=8, pady=2)
        tk.Button(btns, text='Plot', command=self._plot).pack(side=tk.LEFT, padx=4)

        tk.Label(btns, text='  Fit range n:').pack(side=tk.LEFT, padx=4)
        self.nmin_var = tk.StringVar(value='1')
        self.nmax_var = tk.StringVar(value='')
        tk.Entry(btns, textvariable=self.nmin_var, width=5).pack(side=tk.LEFT, padx=1)
        tk.Label(btns, text='-').pack(side=tk.LEFT)
        tk.Entry(btns, textvariable=self.nmax_var, width=5).pack(side=tk.LEFT, padx=1)
        tk.Label(btns, text='(vuoto=auto)').pack(side=tk.LEFT, padx=2)

        tk.Button(btns, text='Fit exp', command=self._fit).pack(side=tk.LEFT, padx=8)
        tk.Button(btns, text='Salva fit', command=self._save_fit).pack(side=tk.LEFT, padx=4)

        self._show_res_var = tk.BooleanVar(value=False)
        self._res_btn = tk.Checkbutton(btns, text='Residui',
                                        variable=self._show_res_var,
                                        command=self._toggle_residuals)
        self._res_btn.pack(side=tk.LEFT, padx=4)
        self._res_btn.config(state='disabled')

        self.info_var = tk.StringVar()
        tk.Label(self, textvariable=self.info_var, fg='darkgreen', anchor='w').pack(fill=tk.X, padx=8)

        bottom = tk.PanedWindow(self, orient=tk.HORIZONTAL)
        bottom.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

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
        data = dm.load_connected_correlators(p['bhw'], p['nstep'], p['nt'], p['therm'])
        return p, data

    def _plot(self):
        p, data = self._get_data()
        if data is None:
            messagebox.showinfo('Info', 'Risultati non disponibili')
            return
        cn = self.corr_var.get()
        eta = float(p['bhw']) / p['nt']

        self.plot.set_layout(1, 1)
        ax = self.plot.axes[0]
        ax.errorbar(data['n'], data[cn], yerr=data[cn + '_err'],
                     **DATA_DOT_KW)
        ax.set_xlabel('$n$')
        ax.set_ylabel(CORR_LABELS[cn])
        ax.set_yscale('log')
        ax.set_title(f'$N_t={p["nt"]}$, $\\eta={eta:.4f}$', fontsize=11)
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
        mean = data[cn]
        err = data[cn + '_err']

        nmin = int(self.nmin_var.get()) if self.nmin_var.get() else 1
        nmax_s = self.nmax_var.get().strip()

        pos = mean > 0
        n_pos, m_pos, e_pos = n[pos], mean[pos], err[pos]

        sel = n_pos >= nmin
        if nmax_s:
            sel &= n_pos <= int(nmax_s)
        else:
            sn = m_pos / e_pos
            auto = sn > 2
            sel &= auto

        nf, mf, ef = n_pos[sel], m_pos[sel], e_pos[sel]
        if len(nf) < 2:
            messagebox.showwarning('Fit', 'Troppo pochi punti per il fit')
            return

        try:
            popt, pcov = curve_fit(exp_model, nf, mf, p0=[mf[0], 5.0],
                                    sigma=ef, absolute_sigma=True,
                                    bounds=([0, 0.1], [np.inf, np.inf]))
            A_fit, xi_fit = popt
            A_err, xi_err = np.sqrt(np.diag(pcov))
            de = 1.0 / (xi_fit * eta)
            de_err = xi_err / (xi_fit**2 * eta)
            residuals = (mf - exp_model(nf, *popt)) / ef
            chi2 = np.sum(residuals**2)
            ndof = len(nf) - 2
            chi2red = chi2 / ndof if ndof > 0 else np.nan
            corr_mat = pcov / np.outer(np.sqrt(np.diag(pcov)), np.sqrt(np.diag(pcov)))
        except Exception as e:
            messagebox.showerror('Fit', str(e))
            return

        self._fit_result = {
            'corr': cn, 'A': float(A_fit), 'A_err': float(A_err),
            'xi': float(xi_fit), 'xi_err': float(xi_err),
            'DeltaE': float(de), 'DeltaE_err': float(de_err),
            'chi2red': float(chi2red), 'ndof': int(ndof),
            'n_min': int(nf[0]), 'n_max': int(nf[-1]),
            'corr_A_xi': float(corr_mat[0, 1]),
            'bhw': p['bhw'], 'nstep': p['nstep'],
            'nt': p['nt'], 'therm': p['therm']
        }

        self._fit_draw_data = {
            'cn': cn, 'eta': eta, 'p': p,
            'n_all': n, 'mean_all': mean, 'err_all': err,
            'nf': nf, 'mf': mf, 'ef': ef,
            'A_fit': A_fit, 'xi_fit': xi_fit,
            'popt': popt, 'residuals': residuals,
        }

        self._show_res_var.set(False)
        self._res_btn.config(state='normal')
        self._draw_fit_full()

        v_xi, e_xi = format_value_with_uncertainty(xi_fit, xi_err)
        v_de, e_de = format_value_with_uncertainty(de, de_err)
        self.info_var.set(
            f'\u03be = {v_xi} \u00b1 {e_xi}   '
            f'\u0394E = {v_de} \u00b1 {e_de}   '
            f'\u03c7\u00b2/ndof = {chi2red:.3f} ({ndof})   '
            f'corr(A,\u03be) = {corr_mat[0,1]:.4f}   '
            f'range: n \u2208 [{nf[0]}, {nf[-1]}]'
        )

    def _draw_fit_full(self):
        """Draw full data with fit line overlaid (no residuals)."""
        d = self._fit_draw_data
        cn = d['cn']
        self.plot.set_layout(1, 1)
        ax = self.plot.axes[0]
        ax.errorbar(d['n_all'], d['mean_all'], yerr=d['err_all'],
                     **DATA_DOT_KW, label='dati')
        nplot = np.linspace(d['nf'][0], d['nf'][-1] * 1.3, 200)
        ax.plot(nplot, exp_model(nplot, *d['popt']), **FIT_LINE_KW,
                label='fit')
        ax.set_xlabel('$n$')
        ax.set_ylabel(CORR_LABELS[cn])
        ax.set_yscale('log')
        ax.legend(fontsize=9)
        ax.set_title(f'$N_t={d["p"]["nt"]}$, $\\eta={d["eta"]:.4f}$', fontsize=11)
        apply_grid(ax)
        self.plot.draw()

    def _draw_fit_with_residuals(self):
        """Draw fit data + residuals, both trimmed to fit range."""
        d = self._fit_draw_data
        cn = d['cn']
        ax_main, ax_res = self.plot.set_layout_with_residuals(height_ratios=(3, 1))
        nf, mf, ef = d['nf'], d['mf'], d['ef']
        ax_main.errorbar(nf, mf, yerr=ef, **DATA_DOT_KW, label='dati')
        nplot = np.linspace(nf[0], nf[-1], 200)
        ax_main.plot(nplot, exp_model(nplot, *d['popt']), **FIT_LINE_KW,
                     label='fit')
        ax_main.set_ylabel(CORR_LABELS[cn])
        ax_main.set_yscale('log')
        ax_main.legend(fontsize=9)
        apply_grid(ax_main)
        xlim = (nf[0] - 0.5, nf[-1] + 0.5)
        ax_main.set_xlim(xlim)
        plot_residuals(ax_res, nf, d['residuals'])
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
        tag = dm.save_fit_entry(base, f'corr_length_{r["corr"]}', r, self.plot.fig)
        messagebox.showinfo('Salvato', f'Fit salvato: {tag}')
        self.fits_panel.refresh()


class MultiNtSubtab(tk.Frame):
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
        tk.Entry(ctrl, textvariable=self.therm_var, width=8).pack(side=tk.LEFT, padx=2)

        tk.Label(ctrl, text='  Correlatore:').pack(side=tk.LEFT, padx=4)
        self.corr_var = tk.StringVar(value='yc')
        for cn in CORR_NAMES:
            tk.Radiobutton(ctrl, text=cn, variable=self.corr_var, value=cn).pack(side=tk.LEFT)

        r2 = tk.Frame(self)
        r2.pack(fill=tk.X, padx=8, pady=2)
        tk.Label(r2, text='Nt da plottare (virgola):').pack(side=tk.LEFT, padx=2)
        self.nts_var = tk.StringVar()
        tk.Entry(r2, textvariable=self.nts_var, width=30).pack(side=tk.LEFT, padx=2)
        tk.Button(r2, text='Plot', command=self._plot).pack(side=tk.LEFT, padx=8)
        tk.Button(r2, text='Salva', command=self._save).pack(side=tk.LEFT, padx=4)

        self.plot = PlotFrame(self)
        self.plot.pack(fill=tk.BOTH, expand=True)

        self._refresh_bhw()
        self.bhw_var.trace_add('write', lambda *_: self._fill_nts())

    def _refresh_bhw(self):
        sets = dm.scan_data_sets()
        self.nstep_map = {str(b): n for b, n in sets}
        menu = self.bhw_cb['menu']
        menu.delete(0, 'end')
        for bhw, _ in sets:
            menu.add_command(label=str(bhw), command=tk._setit(self.bhw_var, str(bhw)))
        if sets:
            self.bhw_var.set(str(sets[0][0]))

    def _fill_nts(self):
        bhw = self.bhw_var.get()
        if not bhw:
            return
        nstep = self.nstep_map.get(bhw)
        nts = dm.get_available_nt(int(bhw), nstep)
        self.nts_var.set(','.join(str(n) for n in nts))

    def _plot(self):
        bhw_s = self.bhw_var.get()
        if not bhw_s:
            return
        bhw = int(bhw_s)
        nstep = self.nstep_map.get(bhw_s)
        therm_s = self.therm_var.get().strip()
        therm = int(therm_s) if therm_s else None
        cn = self.corr_var.get()

        nts_s = self.nts_var.get().strip()
        if not nts_s:
            return
        nts = [int(x.strip()) for x in nts_s.split(',') if x.strip()]

        self.plot.clear()
        ax = self.plot.ax
        for idx, nt in enumerate(nts):
            if therm is not None:
                therms = [therm]
            else:
                therms = dm.get_available_therm(bhw, nstep, nt)
            marker = SMALL_MARKERS[idx % len(SMALL_MARKERS)]
            for th in therms:
                if not dm.results_exist(bhw, nstep, nt, th):
                    continue
                data = dm.load_connected_correlators(bhw, nstep, nt, th)
                eta = float(bhw) / nt
                ax.errorbar(data['n'] * eta, data[cn], yerr=data[cn + '_err'],
                            fmt=marker, ms=SMALL_MS, capsize=2, elinewidth=0.8,
                            label=f'$N_t={nt}$')
        ax.set_xlabel(r'$n\eta$')
        ax.set_ylabel(CORR_LABELS[cn])
        ax.set_yscale('log')
        ax.legend(fontsize=8)
        apply_grid(ax)
        self.plot.draw()

    def _save(self):
        bhw_s = self.bhw_var.get()
        if not bhw_s:
            return
        nstep = self.nstep_map.get(bhw_s)
        cn = self.corr_var.get()
        outdir = os.path.join(dm.PLOTS_DIR, f'bhw{bhw_s}_nstep{nstep}')
        os.makedirs(outdir, exist_ok=True)
        path = os.path.join(outdir, f'multi_nt_{cn}.png')
        self.plot.save(path)
        messagebox.showinfo('Salvato', path)
