# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

import sphinx_autosummary_accessors

import new_majiq

# -- Project information -----------------------------------------------------

project = "MAJIQ"
copyright = "2020, University of Pennsylvania"

# the short x.y version
version = new_majiq.__version__.split("+", 1)[0]
release = new_majiq.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "IPython.sphinxext.ipython_directive",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinx_autosummary_accessors",
    "nbsphinx",
    "sphinxcontrib.programoutput",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates", sphinx_autosummary_accessors.templates_path]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["**.ipynb_checkpoints"]


autosummary_generate = True
autodoc_typehints = "none"


# Napoleon configurations
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = False
napoleon_use_rtype = False
napoleon_preprocess_types = True
napoleon_type_aliases = {
    # without namespace
    "GeneJunctions": "new_majiq.GeneJunctions",
    "GeneJunctionsAccumulator": "new_majiq.GeneJunctionsAccumulator",
    "DPsiPrior": "new_majiq.DPsiPrior",
    "PsiCoverage": "new_majiq.PsiCoverage",
    "SJIntronsBins": "new_majiq.SJIntronsBins",
    "SJJunctionsBins": "new_majiq.SJJunctionsBins",
    # internals
    "_SpliceGraph": "new_majiq.internals.SpliceGraph",
    "_Events": "new_majiq.internals.Events",
    # abbreviated namespace
    "xr.DataArray": "xarray.DataArray",
    "xr.Dataset": "xarray.Dataset",
}


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
html_title = "MAJIQ documentation"
html_static_path = ["_static"]
html_logo = "_static/logo1.png"
html_favicon = "_static/favicon.ico"
html_css_files = ["style.css"]
html_theme_options = {
    "logo_only": True,
    "display_version": True,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".


# options for LaTeX output
latex_elements = {"preamble": r"\usepackage{enumitem}\setlistdepth{99}"}
