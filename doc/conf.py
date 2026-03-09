import os
import sys
from importlib.metadata import version, PackageNotFoundError

# -- Path setup --------------------------------------------------------------

sys.path.insert(0, os.path.abspath(".."))

import sphinx_gallery

project = "bioconvert"
author = "Bioconvert developers"

try:
    release = version("bioconvert")
except PackageNotFoundError:
    release = "0.0.0"

version = release

# -- ReadTheDocs detection ---------------------------------------------------

on_rtd = os.environ.get("READTHEDOCS") == "True"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "numpydoc",
    "matplotlib.sphinxext.plot_directive",
    "sphinx_gallery.gen_gallery",
]

autosummary_generate = True
autoclass_content = "both"

templates_path = ["_templates"]
exclude_patterns = ["_build"]

source_suffix = ".rst"
master_doc = "index"

pygments_style = "sphinx"

# -- Mock heavy dependencies for RTD ----------------------------------------

autodoc_mock_imports = [
    "pysam",
    "pyBigWig",
    "sambamba",
]

# -- Matplotlib --------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")

# -- Sphinx Gallery ----------------------------------------------------------

plot_gallery = not on_rtd

sphinx_gallery_conf = {
    "doc_module": "bioconvert",
    "backreferences_dir": "modules/generated",
}

# -- HTML output -------------------------------------------------------------

html_theme = "sphinx_rtd_theme"

html_static_path = ["_static"]
html_logo = "_static/logo.png"

html_last_updated_fmt = "%b %d, %Y"

# -- Intersphinx -------------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}
