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
import os
import re
import sys

sys.path.insert(0, os.path.abspath("../../src"))

# -- Project information -----------------------------------------------------

project = "GEOUNED"
copyright = "2024, UNED"
author = "Juan-Pablo Catalan and Patrick Sauvan"

# The full version, including alpha/beta/rc tags
import geouned

version = geouned.__version__
release = geouned.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx_autodoc_typehints",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.viewcode",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.
# this theme supports versions https://github.com/pydata/pydata-sphinx-theme
html_theme = "pydata_sphinx_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# TODO add logo
# html_favicon = "favicon.ico"

# Version match must match the 'version' key in version_switcher.json
pattern = re.compile(r"^[0-9]+\.[0-9]+")
version_match = pattern.search(version)
if version_match:
    version_match = version_match.group()
elif "dev" in version:
    version_match = "dev"
else:
    version_match = version

html_theme_options = {
    "github_url": "https://github.com/GEOUNED-org/GEOUNED",
    "switcher": {
        "json_url": "https://raw.githubusercontent.com/fusion-neutronics/GEOUNED/adding_version_support_to_docs/docs/version_switcher.json",
        "version_match": version_match,
    },
    "nav_title": "Geouned",
    "navbar_start": ["version-switcher", "navbar-icon-links"],
}
