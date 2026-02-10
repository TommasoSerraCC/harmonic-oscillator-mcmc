import tkinter as tk
from tkinter import messagebox
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import os
import json
import math

try:
    from PIL import Image, ImageTk
    HAS_PIL = True
except ImportError:
    HAS_PIL = False

from ui.core import data_manager as dm

# ===== Style constants =====
GRID_KW = dict(alpha=0.25, linewidth=0.4, linestyle='-')
FIT_LINE_KW = dict(color='red', linestyle='-', linewidth=0.9)
DATA_DOT_KW = dict(fmt='o', color='black', ms=4, capsize=0.8, elinewidth=0.6, markerfacecolor='none', markeredgewidth=0.6)
RESIDUAL_DOT_KW = dict(marker='o', color='black', ms=4, linestyle='none', zorder=5)
SMALL_MARKERS = ['o', 's', '^', 'v', 'D', '<', '>', 'p', 'h', 'X', 'd', '*']
SMALL_MS = 4


def format_value_with_uncertainty(value, error, n_sig=2):
    """Format value +/- error, rounding to n_sig significant figures of error.
    Returns (val_str, err_str)."""
    if error == 0 or not math.isfinite(error) or not math.isfinite(value):
        return f'{value}', f'{error}'
    mag = math.floor(math.log10(abs(error)))
    decimals = -int(mag) + n_sig - 1
    rounded_err = round(error, decimals)
    rounded_val = round(value, decimals)
    if decimals > 0:
        return f'{rounded_val:.{decimals}f}', f'{rounded_err:.{decimals}f}'
    else:
        return f'{int(rounded_val)}', f'{int(rounded_err)}'


def apply_grid(ax):
    """Apply minimal background grid."""
    ax.grid(True, **GRID_KW)


def plot_residuals(ax, x, residuals):
    """Plot normalized residuals with standard styling."""
    ax.plot(x, residuals, **RESIDUAL_DOT_KW)
    ax.axhline(0, color='k', ls='--', lw=0.8)
    ax.set_ylabel('Residui norm.')
    apply_grid(ax)


class PlotFrame(tk.Frame):
    def __init__(self, parent, figsize=(8, 5), nrows=1, ncols=1):
        super().__init__(parent)
        self.fig = Figure(figsize=figsize, dpi=100)
        if nrows == 1 and ncols == 1:
            self.ax = self.fig.add_subplot(111)
            self.axes = [self.ax]
        else:
            self.axes = []
            for i in range(nrows * ncols):
                self.axes.append(self.fig.add_subplot(nrows, ncols, i + 1))
            self.ax = self.axes[0]
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self)
        self.toolbar.pack(side=tk.BOTTOM, fill=tk.X)

    def clear(self):
        for a in self.axes:
            a.clear()

    def draw(self):
        self.fig.tight_layout()
        self.canvas.draw()

    def save(self, path):
        self.fig.savefig(path, dpi=150, bbox_inches='tight')

    def set_layout(self, nrows, ncols):
        self.fig.clear()
        self.axes = []
        for i in range(nrows * ncols):
            self.axes.append(self.fig.add_subplot(nrows, ncols, i + 1))
        self.ax = self.axes[0]

    def set_layout_with_residuals(self, height_ratios=(3, 1)):
        """Create 2-row layout with fit plot + residuals at given height ratio."""
        self.fig.clear()
        gs = self.fig.add_gridspec(2, 1, height_ratios=height_ratios, hspace=0.35)
        ax_main = self.fig.add_subplot(gs[0])
        ax_res = self.fig.add_subplot(gs[1])
        self.axes = [ax_main, ax_res]
        self.ax = ax_main
        return ax_main, ax_res


class ParamSelector(tk.Frame):
    """Reusable frame for selecting bhw, nt, therm parameters."""
    def __init__(self, parent, dm, on_change=None, show_therm=True):
        super().__init__(parent)
        self.dm = dm
        self.on_change = on_change
        self.show_therm = show_therm

        self.bhw_var = tk.StringVar()
        self.nt_var = tk.StringVar()
        self.therm_var = tk.StringVar()
        self.nstep_map = {}

        row = tk.Frame(self)
        row.pack(fill=tk.X, pady=2)

        tk.Label(row, text='\u03b2\u210f\u03c9:').pack(side=tk.LEFT, padx=2)
        self.bhw_cb = tk.OptionMenu(row, self.bhw_var, '')
        self.bhw_cb.config(width=6)
        self.bhw_cb.pack(side=tk.LEFT, padx=2)

        tk.Label(row, text='Nt:').pack(side=tk.LEFT, padx=2)
        self.nt_cb = tk.OptionMenu(row, self.nt_var, '')
        self.nt_cb.config(width=6)
        self.nt_cb.pack(side=tk.LEFT, padx=2)

        if show_therm:
            tk.Label(row, text='therm:').pack(side=tk.LEFT, padx=2)
            self.therm_cb = tk.OptionMenu(row, self.therm_var, '')
            self.therm_cb.config(width=8)
            self.therm_cb.pack(side=tk.LEFT, padx=2)

        self.bhw_var.trace_add('write', self._bhw_changed)
        self.nt_var.trace_add('write', self._nt_changed)
        if show_therm:
            self.therm_var.trace_add('write', self._therm_changed)

        self.refresh()

    def refresh(self):
        sets = self.dm.scan_data_sets()
        self.nstep_map = {str(bhw): nstep for bhw, nstep in sets}
        menu = self.bhw_cb['menu']
        menu.delete(0, 'end')
        for bhw, _ in sets:
            menu.add_command(label=str(bhw), command=tk._setit(self.bhw_var, str(bhw)))
        if sets:
            self.bhw_var.set(str(sets[0][0]))

    def _bhw_changed(self, *_):
        bhw = self.bhw_var.get()
        if not bhw:
            return
        nstep = self.nstep_map.get(bhw)
        if not nstep:
            return
        nts = self.dm.get_available_nt(int(bhw), nstep)
        menu = self.nt_cb['menu']
        menu.delete(0, 'end')
        for nt in nts:
            menu.add_command(label=str(nt), command=tk._setit(self.nt_var, str(nt)))
        if nts:
            self.nt_var.set(str(nts[0]))

    def _nt_changed(self, *_):
        if not self.show_therm:
            self._notify()
            return
        bhw = self.bhw_var.get()
        nt = self.nt_var.get()
        if not bhw or not nt:
            return
        nstep = self.nstep_map.get(bhw)
        therms = self.dm.get_available_therm(int(bhw), nstep, int(nt))
        menu = self.therm_cb['menu']
        menu.delete(0, 'end')
        for t in therms:
            menu.add_command(label=str(t), command=tk._setit(self.therm_var, str(t)))
        if therms:
            self.therm_var.set(str(therms[0]))
        else:
            self.therm_var.set('')
            self._notify()

    def _therm_changed(self, *_):
        self._notify()

    def _notify(self):
        if self.on_change:
            self.on_change()

    def get_params(self):
        bhw = self.bhw_var.get()
        nt = self.nt_var.get()
        therm = self.therm_var.get() if self.show_therm else '0'
        if not bhw or not nt:
            return None
        nstep = self.nstep_map.get(bhw)
        return {
            'bhw': int(bhw), 'nstep': nstep, 'nt': int(nt),
            'therm': int(therm) if therm else 0
        }

    def get_nstep(self):
        bhw = self.bhw_var.get()
        return self.nstep_map.get(bhw)


class SavedFitsPanel(tk.LabelFrame):
    """Reusable panel for viewing/deleting saved fits."""
    def __init__(self, parent, get_fits_dir):
        super().__init__(parent, text='Fit salvati')
        self.get_fits_dir = get_fits_dir

        self.listbox = tk.Listbox(self, font=('Consolas', 8), height=8)
        sb = tk.Scrollbar(self, orient=tk.VERTICAL, command=self.listbox.yview)
        self.listbox.config(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.listbox.pack(fill=tk.BOTH, expand=True, padx=2, pady=2)
        self.listbox.bind('<Double-1>', lambda e: self._view())

        btn = tk.Frame(self)
        btn.pack(fill=tk.X, pady=2)
        tk.Button(btn, text='Visualizza', command=self._view).pack(side=tk.LEFT, padx=2)
        tk.Button(btn, text='Elimina', command=self._delete).pack(side=tk.LEFT, padx=2)
        tk.Button(btn, text='Aggiorna', command=self.refresh).pack(side=tk.LEFT, padx=2)

        self._tags = []

    def refresh(self):
        self.listbox.delete(0, tk.END)
        fd = self.get_fits_dir()
        if not fd:
            self._tags = []
            return
        self._tags = dm.list_fit_entries(fd)
        for t in self._tags:
            self.listbox.insert(tk.END, t)

    def _get_selected(self):
        sel = self.listbox.curselection()
        if not sel:
            return None
        return self._tags[sel[0]]

    def _view(self):
        tag = self._get_selected()
        if not tag:
            return
        fd = self.get_fits_dir()
        data, png_path = dm.load_fit_entry(fd, tag)

        win = tk.Toplevel(self)
        win.title(tag)
        win.geometry('850x650')

        if png_path and HAS_PIL:
            try:
                img = Image.open(png_path)
                max_w, max_h = 820, 500
                w, h = img.size
                ratio = min(max_w / w, max_h / h, 1.0)
                img = img.resize((int(w * ratio), int(h * ratio)), Image.LANCZOS)
                photo = ImageTk.PhotoImage(img)
                lbl = tk.Label(win, image=photo)
                lbl.image = photo
                lbl.pack(padx=4, pady=4)
            except Exception:
                pass

        if data:
            txt = tk.Text(win, font=('Consolas', 9), height=8, wrap=tk.WORD)
            txt.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)
            txt.insert('1.0', json.dumps(data, indent=2))
            txt.config(state='disabled')

    def _delete(self):
        tag = self._get_selected()
        if not tag:
            return
        if not messagebox.askyesno('Conferma', f'Eliminare il fit "{tag}"?'):
            return
        fd = self.get_fits_dir()
        dm.delete_fit_entry(fd, tag)
        self.refresh()
