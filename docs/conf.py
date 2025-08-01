# Configuration file for the Sphinx documentation builder.

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('..'))

import shutil
import subprocess
from typing import Any, Dict

# Force to adopt the EN_US locale settings to make the data treatment uniform
import locale
us_locale: str = 'en_US.utf8'
try:
    locale.setlocale(locale.LC_ALL, us_locale)
except Exception as e:
    raise ValueError(f"Locale '{us_locale}' not supported!\n" \
                     "Please, install it to generate the doc!")

# -- Source files directory name
SOURCE_DIR_NAME = "source"


# -----------------------------------
# Documentation configuration options
# -----------------------------------
project = "GLOW"
authors = "Davide Manzione"
version = "1.0"
date = "See signature"
doc_title = "GLOW: Generator of Unstructured, Non-Native Geometries for \
Tracking Analyses in DRAGON5"
id_no = "0"
reference_code = "XXX-YYY-ZZZ-???"
revision_no = "1"
modified_pages = "All"
modif_descr = "First Release"
abstract = "This document is the reference manual for the \
\\sphinxstyleemphasis{GLOW} \
(\\sphinxstylestrong{G}eometry \\sphinxstylestrong{L}ayout \
\\sphinxstylestrong{O}riented \\sphinxstylestrong{W}orkflows) Python package, \
which implements a generator of unstructured, non-native geometries for the \
\\sphinxstyleemphasis{Method} \\sphinxstyleemphasis{of} \
\\sphinxstyleemphasis{Characteristics} (MoC) in \
\\sphinxstyleemphasis{DRAGON5}. It uses the APIs of \
\\sphinxstyleemphasis{SALOME} to provide the functionalities to \
build and visualize the geometry, and to generate the corresponding surface \
geometry representation according to the \\sphinxstyleemphasis{APOLLO2} TDT \
format which \\sphinxstyleemphasis{DRAGON5} uses in tracking analyses."
reviewers = "Gabriele Ottino"
approvers = "Daniele Tomatis"
doc_purpose = "1"
business_mark = "4"
exprt_ctrl = "3"
national_sec = {"country": "4", "sec_level": "1", }
itns = "1"
distribution_list = "@newcleo all"
bibtex_bibfiles = [os.path.join(SOURCE_DIR_NAME, 'glow.bib')]
latex_theme_to_use = "manual"


# -- Global substitutions ----------------------------------------------------
rst_prolog = fr"""
.. |TOOL| replace:: **GLOW**
.. |newcleo| replace:: *new*\cleo
.. |LICENSE| replace:: **LGPL-2.1**
"""

# -- General configuration ---------------------------------------------------

# Extensions
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'myst_parser',
    'sphinx.ext.autosectionlabel',
    'sphinxcontrib.bibtex',
    "sphinx.ext.extlinks",
]

myst_enable_extensions = [
    "dollarmath",    # To enable parsing of inline and block math using dollar signs
    "amsmath",       # To include more complex mathematical expressions
    "deflist",       # To include definition lists for creating glossaries
    "colon_fence",   # To make code blocks more readable in the documentation
]

# External links configuration
extlinks = {}

# Labels automations
numfig = True
math_numfig = True

# Code documentation settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
add_module_names = False

autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'private-members': False,
    'inherited-members': False,
    'show-inheritance': True,
}
autodoc_member_order = 'groupwise'
autodoc_mock_imports = ["salome", "SALOMEDS", "GEOM"]
autosectionlabel_prefix_document = True

autodoc_preserve_defaults = True

# Sources definition
source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

# The master toctree document
master_doc = "index"

language = "en"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# -- Options of MyST-NB ------------------------------------------------------
# Execute notebooks only if output files not present.
# nb_execution_mode = "auto"

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
    "navigation_depth": 5,
}

# Add any path that contains custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = [os.path.join(SOURCE_DIR_NAME, "_static")]

html_css_files = [
    os.path.join("css", "theme.css"),
    os.path.join("css", "eqno.css"),
]

html_logo = os.path.join(SOURCE_DIR_NAME, "images", "newcleo_logo.png")


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = project + " documentation"


# -- Options for LaTeX output ------------------------------------------------

# List of folders where latex templates are placed
templates_path = ["_templates"]

# Mapping of the available latex docclasses
latex_docclasses = {
    "manual": "report",
    "nwcldocs": "nwcldocs"
}

# Extract the list of authors
authors_list = [str(item) + 3*r"\break " for item in authors.split(",")]

#############################
# Setup for "manual" template
latex_theme = "manual"
latex_elements: Dict[str, Any] = {
    "papersize": "letterpaper",
    "pointsize": "11pt",
    "extrapackages": r"\input{../../_templates/extra_manual.texsty}",
    "makeindex": "\\usepackage[columns=1]{idxlayout}\\makeindex",
    "figure_align": "H",
}
latex_authors = "".join(str(item) for item in authors_list)

##########################################
# Settings to use with "nwcldocs" template
extrapackages_nwcldocs = r"\input{../../_templates/extra_nwcl.texsty}"
preamble_nwcldocs = r"\input{../../_templates/preamble_nwcl.texsty}"
makeindex_nwcldocs = "\\usepackage[columns=1,totoc]{idxlayout}\\makeindex"
date_nwcldocs = date

# Extract the list of reviewers and approvers
reviewers_list = [str(item) + 3*r"\break " for item in reviewers.split(",")]
approvers_list = [str(item) + 3*r"\break " for item in approvers.split(",")]

# date_nwcldocs = date_object.strftime("%d/%m/%Y")
atendofbody_nwcldocs = {
    "id": id_no,
    "ref_code": reference_code,
    "rev": revision_no,
    "abstract": abstract,
    "date": date_nwcldocs,
    "pages": modified_pages,
    "desc": modif_descr,
    "reviewers": "".join(str(item) for item in reviewers_list),
    "approvers": "".join(str(item) for item in approvers_list),
    "doc_for": doc_purpose,
    "bus_mark": business_mark,
    "exp_ctrl": exprt_ctrl,
    "national_sec": [item for item in national_sec.values()],
    "itns": itns,
    "dist_list": distribution_list
}

# Check the correspondence of the indicated theme with the currently available
# ones
if latex_theme_to_use not in latex_docclasses:
    raise RuntimeError("Only 'manual' and 'nwcldocs' document classes are "
                       f"allowed, while '{latex_theme_to_use}' has been "
                       "provided!\n Please change document class to use "
                       "and try again!")

# Check the availability of the 'nwcldocs' theme and set it accordingly
if latex_theme_to_use == "nwcldocs":
    if shutil.which("kpsewhich"):
        check = subprocess.run(["kpsewhich", "nwcldocs.cls"],
                               stdout=subprocess.PIPE, text=True)
        if not check.stdout:
            raise RuntimeError("'nwcldocs.cls' tex class not found!\n"
                               "PDF file cannot be generated!\n Please "
                               "change document class to use and try again!")
        latex_theme = "nwcldocs"
        latex_toplevel_sectioning = "section"
        latex_elements.pop("papersize", None)
        latex_elements["classoptions"] = "techdoc"
        latex_elements["extrapackages"] = extrapackages_nwcldocs
        latex_elements["preamble"] = preamble_nwcldocs
        latex_elements["makeindex"] = makeindex_nwcldocs
        latex_elements["atendofbody"] = atendofbody_nwcldocs
        latex_elements["sphinxsetup"] = "TitleColor={named}{black}"
        latex_authors = "".join(str(item) for item in authors_list)
    else:
        raise RuntimeError(
            "No tex environment found!\n"
            "PDF file cannot be generated!\n Please change "
            "document class to use and try again!"
        )

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
# author, documentclass [howto, manual, or own class]).
latex_documents = [
    (
        master_doc,
        "glow.tex",
        doc_title,
        latex_authors,
        latex_theme,
    ),
]

latex_logo = html_logo
