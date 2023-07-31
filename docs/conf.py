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


# -- Project information -----------------------------------------------------

project = 'easyunfold'
copyright = '2023, Bonan Zhu'
author = 'Bonan Zhu'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc', "myst_parser", "sphinx_click", "autodoc2", 'sphinx_copybutton'
]

autodoc2_packages = [
    "../easyunfold",
]

myst_enable_extensions = [
    "amsmath",
    "attrs_inline",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "replacements",
    "smartquotes",
    "linkify",
    "strikethrough",
    "substitution",
    "tasklist",
]

autodoc2_render_plugin = "myst"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_book_theme'
html_theme_options = {
  'navigation_depth': 2,
  'repository_url': 'https://github.com/SMTG-UCL/easyunfold',
  'use_repository_button': True,
  "show_navbar_depth": 2,
  'home_page_in_toc': True,
  "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/SMTG-UCL/easyunfold",  # required
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/easyunfold/",
            "icon": "https://img.shields.io/pypi/v/easyunfold",
            "type": "url",
        }
    ]
}

html_logo = 'img/logo.svg'
html_title = "easyunfold"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']