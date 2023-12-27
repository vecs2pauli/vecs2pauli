# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import subprocess, os
# Doxygen
subprocess.call('doxygen doxyfile.in', shell=True)

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'vecs2pauli'
copyright = '2023, vecs2pauli'
author = 'vecs2pauli'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'breathe'
              ]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


highlight_language = 'c++'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
# html_static_path = ['_static']

# -- Breathe configuration -------------------------------------------------

breathe_projects = {
	"vecs2pauli": "_build/xml/"
}
breathe_domain_by_extension = {"h" : "cpp"}
breathe_default_project = "vecs2pauli"
