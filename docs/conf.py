# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Add path to module
#sys.path.insert(0, os.path.abspath('/home/mayra/source/Nichesphere/nichesphere/nichesphere'))
#sys.path.insert(0, os.path.abspath('../'))
#sys.path.insert(0, os.path.abspath('../../'))
sys.path.insert(0, os.path.abspath('../src'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'NicheSphere'
copyright = '2025, Mayra Ruiz, James Nagai'
author = 'Mayra Ruiz, James Nagai'
release = '1.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx_rtd_theme', 'nbsphinx', 'sphinx.ext.autosummary', 'sphinx.ext.intersphinx']

templates_path = ['_templates']
exclude_patterns = []

language = 'English'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#import sphinx_rtd_theme

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Path to the logo relative to the configuration directory
html_logo = "_static/logo.png"

# Ensure the static path is included so Sphinx finds the file
html_static_path = ['_static']

# (Optional) If you want ONLY the logo to show without the project name text
html_theme_options = {
    'logo_only': True,
    'display_version': False,
}


