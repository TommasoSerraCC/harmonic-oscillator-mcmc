import tkinter as tk
from tkinter import ttk, messagebox
import threading
import subprocess
import sys
import os
import re

from ui.core import data_manager as dm


class AnalysisTab(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        self.root_dir = dm.ROOT

        paned = tk.PanedWindow(self, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        # ===== Left panel: actions =====
        left = tk.Frame(paned)
        paned.add(left, width=360)

        # -- Acquisizione dati --
        acq = tk.LabelFrame(left, text='Acquisizione dati (Fortran)')
        acq.pack(fill=tk.X, padx=4, pady=4)

        r = tk.Frame(acq)
        r.pack(fill=tk.X, padx=4, pady=2)
        tk.Label(r, text='βℏω:').pack(side=tk.LEFT)
        self.acq_bhw = tk.StringVar(value='10')
        tk.Entry(r, textvariable=self.acq_bhw, width=6).pack(side=tk.LEFT, padx=2)
        tk.Label(r, text='nstep:').pack(side=tk.LEFT, padx=4)
        self.acq_nstep = tk.StringVar(value='1000000')
        tk.Entry(r, textvariable=self.acq_nstep, width=10).pack(side=tk.LEFT, padx=2)

        r2 = tk.Frame(acq)
        r2.pack(fill=tk.X, padx=4, pady=2)
        tk.Label(r2, text='Nt (virgola):').pack(side=tk.LEFT)
        self.acq_nts = tk.StringVar(value='4,12,24,30,36,42,50,75,100,150,200')
        tk.Entry(r2, textvariable=self.acq_nts, width=32).pack(side=tk.LEFT, padx=2)

        self.acq_btn = tk.Button(acq, text='Avvia acquisizione', command=self._run_fortran)
        self.acq_btn.pack(pady=4)

        self.energy_only_var = tk.BooleanVar(value=False)
        tk.Checkbutton(acq, text='Solo energia (no correlatori)',
                       variable=self.energy_only_var).pack(pady=2)

        # -- Analisi primaria --
        ana = tk.LabelFrame(left, text='Analisi primaria (Python)')
        ana.pack(fill=tk.X, padx=4, pady=4)

        r3 = tk.Frame(ana)
        r3.pack(fill=tk.X, padx=4, pady=2)
        tk.Label(r3, text='βℏω:').pack(side=tk.LEFT)
        self.ana_bhw_var = tk.StringVar()
        self.ana_bhw_cb = tk.OptionMenu(r3, self.ana_bhw_var, '')
        self.ana_bhw_cb.config(width=6)
        self.ana_bhw_cb.pack(side=tk.LEFT, padx=2)
        tk.Label(r3, text='skip:').pack(side=tk.LEFT, padx=4)
        self.ana_skip = tk.StringVar(value='50000')
        tk.Entry(r3, textvariable=self.ana_skip, width=8).pack(side=tk.LEFT, padx=2)

        r4 = tk.Frame(ana)
        r4.pack(fill=tk.X, padx=4, pady=2)
        tk.Label(r4, text='Nt da analizzare:').pack(side=tk.LEFT)
        self.ana_list = tk.Listbox(r4, height=5, selectmode=tk.EXTENDED, width=28)
        self.ana_list.pack(side=tk.LEFT, padx=2, fill=tk.X, expand=True)

        r5 = tk.Frame(ana)
        r5.pack(fill=tk.X, padx=4, pady=2)
        tk.Button(r5, text='Analizza selezionati', command=self._run_analysis_selected).pack(side=tk.LEFT, padx=4)
        tk.Button(r5, text='Analizza tutti mancanti', command=self._run_analysis_all).pack(side=tk.LEFT, padx=4)

        self.ana_energy_only_var = tk.BooleanVar(value=False)
        tk.Checkbutton(r5, text='Solo energia',
                       variable=self.ana_energy_only_var,
                       command=lambda: self._refresh_unanalyzed()).pack(side=tk.LEFT, padx=6)

        self.ana_bhw_var.trace_add('write', lambda *_: self._refresh_unanalyzed())

        self.status_var = tk.StringVar()
        tk.Label(left, textvariable=self.status_var, fg='blue', wraplength=340,
                 anchor='w', justify='left').pack(fill=tk.X, padx=8, pady=4)

        # ===== Right panel: file explorers =====
        right = tk.Frame(paned)
        paned.add(right, stretch='always')

        nb = ttk.Notebook(right)
        nb.pack(fill=tk.BOTH, expand=True)

        data_frame = tk.Frame(nb)
        nb.add(data_frame, text='data/')
        self.data_tree = ttk.Treeview(data_frame, show='tree')
        dsb = ttk.Scrollbar(data_frame, orient=tk.VERTICAL, command=self.data_tree.yview)
        self.data_tree.config(yscrollcommand=dsb.set)
        dsb.pack(side=tk.RIGHT, fill=tk.Y)
        self.data_tree.pack(fill=tk.BOTH, expand=True)

        res_frame = tk.Frame(nb)
        nb.add(res_frame, text='results/')
        self.res_tree = ttk.Treeview(res_frame, show='tree')
        rsb = ttk.Scrollbar(res_frame, orient=tk.VERTICAL, command=self.res_tree.yview)
        self.res_tree.config(yscrollcommand=rsb.set)
        rsb.pack(side=tk.RIGHT, fill=tk.Y)
        self.res_tree.pack(fill=tk.BOTH, expand=True)

        tk.Button(right, text='Aggiorna', command=self._refresh_all).pack(pady=4)

        self._refresh_all()

    def _refresh_all(self):
        self._build_dir_tree(self.data_tree, dm.DATA_DIR)
        self._build_dir_tree(self.res_tree, dm.RESULTS_DIR)
        self._refresh_bhw_menu()
        self._refresh_unanalyzed()

    def _build_dir_tree(self, tree, root_path):
        tree.delete(*tree.get_children())
        if not os.path.isdir(root_path):
            return
        self._insert_dir(tree, '', root_path)

    def _insert_dir(self, tree, parent, dirpath):
        try:
            entries = sorted(os.listdir(dirpath))
        except OSError:
            return
        for entry in entries:
            full = os.path.join(dirpath, entry)
            if os.path.isdir(full):
                node = tree.insert(parent, 'end', text=entry + '/')
                self._insert_dir(tree, node, full)
            else:
                size = os.path.getsize(full)
                if size > 1024*1024:
                    sz = f'{size/(1024*1024):.1f} MB'
                elif size > 1024:
                    sz = f'{size/1024:.1f} KB'
                else:
                    sz = f'{size} B'
                tree.insert(parent, 'end', text=f'{entry}  ({sz})')

    def _refresh_bhw_menu(self):
        sets = dm.scan_data_sets()
        menu = self.ana_bhw_cb['menu']
        menu.delete(0, 'end')
        for bhw, nstep in sets:
            label = f'{bhw} (nstep={nstep})'
            menu.add_command(label=label, command=tk._setit(self.ana_bhw_var, str(bhw)))
        if sets:
            self.ana_bhw_var.set(str(sets[0][0]))

    def _refresh_unanalyzed(self):
        self.ana_list.delete(0, tk.END)
        bhw_s = self.ana_bhw_var.get()
        if not bhw_s:
            return
        bhw = int(bhw_s)
        nstep_map = {b: n for b, n in dm.scan_data_sets()}
        nstep = nstep_map.get(bhw)
        if not nstep:
            return
        skip_s = self.ana_skip.get().strip()
        skip = int(skip_s) if skip_s else 0
        if self.ana_energy_only_var.get():
            unanalyzed = dm.get_unanalyzed_energy_only(bhw, nstep, skip)
        else:
            unanalyzed = dm.get_unanalyzed_raw(bhw, nstep, skip)
        for nt in unanalyzed:
            self.ana_list.insert(tk.END, str(nt))
        if not unanalyzed:
            self.ana_list.insert(tk.END, '(tutti analizzati)')

    def _get_ana_nstep(self):
        bhw = int(self.ana_bhw_var.get())
        return dict(dm.scan_data_sets()).get(bhw)

    def _run_analysis_selected(self):
        sel = self.ana_list.curselection()
        if not sel:
            messagebox.showinfo('Info', 'Selezionare almeno un Nt')
            return
        nts = []
        for i in sel:
            val = self.ana_list.get(i)
            if val.startswith('('):
                return
            nts.append(int(val))
        self._run_analysis_batch(nts)

    def _run_analysis_all(self):
        bhw_s = self.ana_bhw_var.get()
        if not bhw_s:
            return
        bhw = int(bhw_s)
        nstep = self._get_ana_nstep()
        skip = int(self.ana_skip.get()) if self.ana_skip.get().strip() else 0
        nts = dm.get_unanalyzed_raw(bhw, nstep, skip)
        if not nts:
            messagebox.showinfo('Info', 'Tutti i dati sono già analizzati')
            return
        self._run_analysis_batch(nts)

    def _run_analysis_batch(self, nts):
        bhw = int(self.ana_bhw_var.get())
        nstep = self._get_ana_nstep()
        skip = int(self.ana_skip.get()) if self.ana_skip.get().strip() else 0
        is_energy_only = self.ana_energy_only_var.get()
        self.status_var.set(f'Analisi in corso: {len(nts)} file...')
        threading.Thread(target=self._run_analysis_thread,
                         args=(bhw, nstep, skip, nts, is_energy_only), daemon=True).start()

    def _run_analysis_thread(self, bhw, nstep, skip, nts, energy_only=False):
        if energy_only:
            for i, nt in enumerate(nts):
                self.after(0, lambda i=i, nt=nt: self.status_var.set(
                    f'Analisi energia {i+1}/{len(nts)}: nt={nt}...'))
                try:
                    dm.analyze_energy_only(bhw, nstep, nt, skip)
                except Exception as e:
                    self.after(0, lambda e=e: messagebox.showerror(
                        'Errore', f'nt={nt}: {str(e)[:400]}'))
                    self.after(0, lambda: self.status_var.set('Errore'))
                    return
            self.after(0, lambda: self.status_var.set(f'Completato ({len(nts)} file energia)'))
            self.after(0, self._refresh_all)
            return
        script = os.path.join(self.root_dir, 'python_scripts', 'analyze_and_save.py')
        for i, nt in enumerate(nts):
            self.after(0, lambda i=i, nt=nt: self.status_var.set(
                f'Analisi {i+1}/{len(nts)}: nt={nt}...'))
            cmd = [sys.executable, script,
                   '--bhw', str(bhw), '--nt', str(nt),
                   '--nstep', str(nstep), '--skip', str(skip)]
            try:
                subprocess.run(cmd, cwd=self.root_dir, check=True,
                               capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                self.after(0, lambda e=e: messagebox.showerror(
                    'Errore', f'nt={nt}: {e.stderr[:400]}'))
                self.after(0, lambda: self.status_var.set('Errore'))
                return
        self.after(0, lambda: self.status_var.set(f'Completato ({len(nts)} file)'))
        self.after(0, self._refresh_all)

    def _run_fortran(self):
        bhw_s = self.acq_bhw.get().strip()
        nstep_s = self.acq_nstep.get().strip()
        nts_s = self.acq_nts.get().strip()
        if not bhw_s or not nstep_s or not nts_s:
            messagebox.showerror('Errore', 'Compilare tutti i campi')
            return
        nts = [x.strip() for x in nts_s.split(',') if x.strip()]
        self.acq_btn.config(state='disabled')
        self.status_var.set('Acquisizione Fortran in corso...')
        threading.Thread(target=self._run_fortran_thread,
                         args=(bhw_s, nstep_s, nts), daemon=True).start()

    def _run_fortran_thread(self, bhw_s, nstep_s, nts):
        exe = os.path.join(self.root_dir, 'main.exe')
        if not os.path.isfile(exe):
            exe = os.path.join(self.root_dir, 'main')
            if not os.path.isfile(exe):
                self.after(0, lambda: messagebox.showerror(
                    'Errore', 'Eseguibile main non trovato. Compilare main.f'))
                self.after(0, lambda: self.acq_btn.config(state='normal'))
                self.after(0, lambda: self.status_var.set(''))
                return
        energy_only = 1 if self.energy_only_var.get() else 0
        n_nt = len(nts)
        nt_lines = ', '.join(f'nt_vals({i+1})={v}' for i, v in enumerate(nts))
        namelist = f'&params bhw={bhw_s}, nsteps={nstep_s}, n_nt={n_nt}, {nt_lines}, energy_only={energy_only} /\n'
        try:
            subprocess.run(exe, input=namelist, cwd=self.root_dir, check=True,
                           capture_output=True, text=True, timeout=7200)
            self.after(0, lambda: self.status_var.set('Acquisizione completata'))
        except subprocess.CalledProcessError as e:
            self.after(0, lambda: messagebox.showerror('Errore', e.stderr[:500]))
            self.after(0, lambda: self.status_var.set('Errore'))
        except subprocess.TimeoutExpired:
            self.after(0, lambda: self.status_var.set('Timeout'))
        finally:
            self.after(0, lambda: self.acq_btn.config(state='normal'))
            self.after(0, self._refresh_all)
