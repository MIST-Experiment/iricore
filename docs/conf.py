import os
import sys

sys.path.insert(0, os.path.abspath('../src'))

# -- Project information -----------------------------------------------------

project = 'iricore'
copyright = '2023, Vadym Bidula'
author = 'Vadym Bidula'
email = 'vadym.bidula@gmail.com'

# The full version, including alpha/beta/rc tags
release = '1.8.3'

# -- General configuration ---------------------------------------------------

html_static_path = []
typehints_defaults = 'comma'
autodoc_typehints = 'description'
simplify_optional_unions = False
html_theme = "sphinx_rtd_theme"

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx_autodoc_typehints',
    'sphinxemoji.sphinxemoji',
    # 'sphinxcontrib.bibtex',
    'sphinx_rtd_theme',
]

