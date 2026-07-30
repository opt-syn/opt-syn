import os
import sys
from datetime import date, datetime, timezone

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'opt-syn'
author = "opt-syn team"
copyright = f"{date.today().year}, opt-syn team"

release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

#configuration options are based on https://github.com/AutoLyap/AutoLyap/blob/main/docs/source/conf.py

extensions = [
    'sphinx.ext.viewcode', 
    'sphinxcontrib.matlab', 
    'sphinx.ext.autodoc',  
    'sphinx.ext.napoleon',  
    'sphinx.ext.mathjax',
    'sphinx_favicon',     
    'sphinxcontrib.bibtex',    
    'myst_parser'
]

html_context = {
    "seo_site_name":"opt-syn",
    "seo_site_description": ("opt-syn is a MATLAB package for the analysis and design of first-order methods"),
    "seo_default_keywords": [
        "opt-syn",
        "optimization",
        "first-order",
        "design",
        "synthesis",
        "convergence",
        "LMI",
        "semidefinite programming"
    ]
}

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# matlab source code
matlab_src_dir = '../src'

primary_domain = 'mat'
matlab_auto_link ="basic"

suppress_warnings = ["myst.xref_missing", "docutils"]

napoleon_custom_sections = [('Returns', 'params_style'), ('Return', 'params_style')]

autodoc_default_options = {
    # other options
    'show-inheritance': True
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#documentation for theme https://sphinx-book-theme.readthedocs.io/en/latest/
html_theme = 'sphinx_book_theme'
html_static_path = ['_static']

favicons = [    
    "favicon/sun_part_16.svg"
]
 
 
html_theme_options = {
   "logo": {
      "image_light": "opt_syn_light.svg",
      "image_dark": "opt_syn_dark.svg",
   },
   "repository_url": "https://github.com/jarmill/opt-syn",
                         "use_repository_button": True
}

html_title = "opt-syn"

html_sidebars={}


#myst options
myst_enable_extensions = [
    "alert",
    "amsmath",
    "attrs_inline",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "gfm_autolink",
    "html_admonition",
    "html_image",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
    "substitution"
]

myst_substitutions = {
    "osyn" : '<span style="color:#5D93BF">opt</span>-syn'
}

#math options
mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-mml-chtml.js"

mathjax3_config = {
    "loader": {
        "load": ["[tex]/html"],
    },
    "tex" : {
        "packages": {"[+]": ["html", "arydshln", "xcolor"]},
        "macros":{
            "abs": [r"\left\lvert #1 \right\rvert", 1],
            "Bignorm": [r"\left\lVert #1 \right\rVert", 1],
            "norm": [r"\lVert #1 \rVert", 1],
            "Biginner": [r"\left\langle #1, #2 \right\rangle", 2],
            "Biginner": [r"\left\langle #1, #2 \right\rangle", 2],
            "mas": [r"\left[\begin{array}{#1}#2\end{array}\right]", 2],
            "mat": [r"\left(\begin{array}{#1}#2\end{array}\right)", 2],
            "mav": [r"\left\lVert\begin{array}{#1}#2\end{array}\right\rVert", 2],
            "argmax": r"\operatorname*{argmax}",
            "argmin": r"\operatorname*{argmin}",
            "R" : r"\mathbb{R}",
            "Rbar": r"\overline{\mathbb{R}}",
            "rhoi" : r"\rho^{-1}",
            "C" : r"\mathbb{C}",
            "sC" : r"\mathcal{C}",
            "N" : r"\mathbb{N}",
            "Z" : r"\mathbb{Z}",
            "F" : r"\mathcal{F}",
            "1" : r"\mathbb{1}",
            "0" : r"\mathbb{0}",
            "hdots" : r"\cdot \cdot \cdot ",
            "hl" : r"\\ \hline",                    
            "hdl" : r"\\ \hdashline", 
            "nulls" : r"\textrm{null}",
            "ov": [r"\overline{#1}", 1],
            "Acl" : r"\mathcal{A}",
            "Bcl" : r"\mathcal{B}",
            "Ccl" : r"\mathcal{C}",
            "Dcl" : r"\mathcal{D}",
            "Ac" : r"A_c",
            "Bc" : r"B_c", #change to make colors
            "Cc" : r"C_c",
            "Dc" : r"D_c", 
            "Kc" : r"K_c", 
            "blam" : r"{\lambda}", 
        }
    }
}

#bibtex
bibtex_bibfiles = ['references.bib']