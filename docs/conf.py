import os
import sys

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'opt-syn'
copyright = '2026, Jared Miller'
author = 'Jared Miller'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.viewcode', 
    'sphinxcontrib.matlab', 
    'sphinx.ext.autodoc',  
    'sphinx_favicon', 
    'myst_parser'
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# matlab source code
matlab_src_dir = '../src'

primary_domain = 'mat'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#httpsdocumentation for theme ://sphinx-book-theme.readthedocs.io/en/latest/
html_theme = 'sphinx_book_theme'
html_static_path = ['_static']

favicons = [    
    "favicon/sun_part_16.svg"
]
#"img/logo/sun_part.icofavicon-16x16.png",
    #"favicon-32x32.png",