# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------

project = u"pytc2"
copyright = u"2023, Mariano Llamedo Soria"
author = u"Mariano Llamedo Soria"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
Sphinx = "^5.3.0"
docutils = "^0.18.1"


extensions = [
    "myst_nb",
    "autoapi.extension",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]
autoapi_dirs = ["../src"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

myst_enable_extensions = [
    "amsmath",
    "html_image",
    "dollarmath"
]

# myst_heading_anchors = 1

suppress_warnings = ["myst.xref_missing","myst.header"]

myst_url_schemes = ("http", "https", "mailto")

#html_extra_path =  ['img']


