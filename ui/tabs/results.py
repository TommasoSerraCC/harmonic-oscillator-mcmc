import tkinter as tk
from tkinter import ttk
import os

from ui.core import data_manager as dm

try:
    from PIL import Image, ImageTk
    HAS_PIL = True
except ImportError:
    HAS_PIL = False


class ResultsTab(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent)

        paned = tk.PanedWindow(self, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        left = tk.Frame(paned)
        paned.add(left, width=300)

        self.tree = ttk.Treeview(left, show='tree')
        sb = ttk.Scrollbar(left, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.config(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree.pack(fill=tk.BOTH, expand=True)
        self.tree.bind('<<TreeviewSelect>>', self._on_select)

        tk.Button(left, text='Aggiorna', command=self._build_tree).pack(pady=4)

        right = tk.Frame(paned)
        paned.add(right, stretch='always')

        self.content = tk.Text(right, wrap=tk.WORD, font=('Consolas', 9))
        sb2 = ttk.Scrollbar(right, orient=tk.VERTICAL, command=self.content.yview)
        self.content.config(yscrollcommand=sb2.set)
        sb2.pack(side=tk.RIGHT, fill=tk.Y)
        self.content.pack(fill=tk.BOTH, expand=True)

        self._img_ref = None
        self._build_tree()

    def _build_tree(self):
        self.tree.delete(*self.tree.get_children())
        self._paths = {}

        for base_dir, label in [(dm.RESULTS_DIR, 'results'), (dm.PLOTS_DIR, 'plots')]:
            if not os.path.isdir(base_dir):
                continue
            root_id = self.tree.insert('', 'end', text=label, open=True)
            self._insert_dir(root_id, base_dir)

    def _insert_dir(self, parent_id, dirpath):
        try:
            entries = sorted(os.listdir(dirpath))
        except OSError:
            return
        for entry in entries:
            full = os.path.join(dirpath, entry)
            if os.path.isdir(full):
                node = self.tree.insert(parent_id, 'end', text=entry)
                self._paths[node] = full
                self._insert_dir(node, full)
            else:
                node = self.tree.insert(parent_id, 'end', text=entry)
                self._paths[node] = full

    def _on_select(self, event):
        sel = self.tree.selection()
        if not sel:
            return
        path = self._paths.get(sel[0])
        if not path or os.path.isdir(path):
            return
        self.content.delete('1.0', tk.END)
        self._img_ref = None

        if path.lower().endswith('.png'):
            self._show_image(path)
        elif path.lower().endswith('.dat') or path.lower().endswith('.txt'):
            self._show_text(path)

    def _show_text(self, path):
        try:
            with open(path, encoding='utf-8', errors='replace') as f:
                text = f.read()
            self.content.insert('1.0', text)
        except Exception as e:
            self.content.insert('1.0', f'Errore: {e}')

    def _show_image(self, path):
        if not HAS_PIL:
            self.content.insert('1.0', f'Pillow non installato.\n{path}')
            return
        try:
            img = Image.open(path)
            max_w, max_h = 700, 500
            w, h = img.size
            ratio = min(max_w / w, max_h / h, 1.0)
            img = img.resize((int(w * ratio), int(h * ratio)), Image.LANCZOS)
            self._img_ref = ImageTk.PhotoImage(img)
            self.content.image_create('1.0', image=self._img_ref)
        except Exception as e:
            self.content.insert('1.0', f'Errore: {e}')
