"""Sphinx configuration for the harmonic oscillator MCMC documentation."""

import os
import sys

sys.path.insert(0, os.path.abspath('..'))

# -- Project information ---------------------------------------------

project = 'Harmonic Oscillator MCMC'
author = 'Tommaso Serra'
copyright = '2026, Tommaso Serra'
release = '1.0'

# -- General configuration -------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'myst_parser',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Autodoc ----------------------------------------------------------

# The GUI modules import tkinter and the Tk matplotlib backend, neither
# of which is available on a headless documentation builder. Mocking them
# lets autodoc read the docstrings without a display.
autodoc_mock_imports = [
    'tkinter',
    'PIL',
    'matplotlib.backends.backend_tkagg',
]

autodoc_member_order = 'bysource'
autodoc_default_options = {
    'members': True,
    'show-inheritance': True,
}

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = False

# -- MyST -------------------------------------------------------------

myst_enable_extensions = [
    'amsmath',
    'colon_fence',
    'dollarmath',
]

myst_heading_anchors = 3

# -- HTML output ------------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = []
html_title = 'Harmonic Oscillator MCMC'
