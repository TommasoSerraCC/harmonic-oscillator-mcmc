import tkinter as tk
from tkinter import ttk

from ui.tabs.analysis import AnalysisTab
from ui.tabs.observables import ObservablesTab
from ui.tabs.correlators import CorrelatorsTab
from ui.tabs.energy_gaps import EnergyGapsTab
from ui.tabs.energy_vs_temp import EnergyVsTempTab
from ui.tabs.results import ResultsTab


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('Oscillatore Armonico — Analisi MCMC')
        self.geometry('1100x700')

        notebook = ttk.Notebook(self)
        notebook.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        notebook.add(AnalysisTab(notebook), text='Analisi')
        notebook.add(ObservablesTab(notebook), text='Osservabili')
        notebook.add(CorrelatorsTab(notebook), text='Correlatori')
        notebook.add(EnergyGapsTab(notebook), text='Gap Energetici')
        notebook.add(EnergyVsTempTab(notebook), text='E vs T')
        notebook.add(ResultsTab(notebook), text='Risultati')


def main():
    app = App()
    app.mainloop()


if __name__ == '__main__':
    main()
