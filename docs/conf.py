# Configuration file for the Sphinx documentation builder.

# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/page/usage/configuration.html

# -- Path setup --------------------------------------------------------------
import shutil
import sys
from datetime import datetime
from importlib.metadata import metadata
from pathlib import Path

from sphinxcontrib import katex

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE / "extensions"))


# -- Project information -----------------------------------------------------

# NOTE: If you installed your project in editable mode, this might be stale.
#       If this is the case, reinstall it to refresh the metadata
info = metadata("NicheSphere")
project = info["Name"]
author = info["Author"]
copyright = f"{datetime.now():%Y}, {author}."
version = info["Version"]
urls = dict(pu.split(", ") for pu in info.get_all("Project-URL"))
repository_url = urls["Source"]

# The full version, including alpha/beta/rc tags
release = info["Version"]

# -- BibTeX citation settings ------------------------------------------------
bibtex_bibfiles = ["references.bib"]
bibtex_default_style = "alpha"
bibtex_reference_style = "author_year"

templates_path = ["_templates"]
nitpicky = True  # Warn about broken links
needs_sphinx = "4.0"

html_context = {
    "display_github": True,  # Integrate GitHub
    "github_user": "CostaLab",
    "github_repo": project,
    "github_version": "main",  # change if your default branch is not `main`
    "conf_py_path": "/docs/",
}

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings.
# They can be extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "myst_nb",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinxcontrib.bibtex",
    "sphinxcontrib.katex",
    "sphinx_autodoc_typehints",
    "sphinx_design",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinxext.opengraph",
    "scverse_misc.sphinx_ext",
    *[p.stem for p in (HERE / "extensions").glob("*.py")],
]

autosummary_generate = True
autodoc_member_order = "groupwise"
default_role = "literal"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_preprocess_types = True   # needed for napoleon_type_aliases to apply
myst_heading_anchors = 6  # create anchors for h1-h6
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
    "html_admonition",
]
myst_url_schemes = ("http", "https", "mailto")
nb_output_stderr = "remove"
nb_execution_mode = "off"  # tutorials need external data: use the stored outputs
nb_merge_streams = True
typehints_defaults = "braces"
always_use_bars_union = True  # use `|` instead of `Union` in types even when building with Python ≤3.14

source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
    ".myst": "myst-nb",
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "anndata": ("https://anndata.scverse.org/en/stable/", None),
    "scanpy": ("https://scanpy.scverse.org/en/stable/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "networkx": ("https://networkx.org/documentation/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

# Map the shorthand used in your docstrings to real, linkable names
napoleon_type_aliases = {
    "pd.DataFrame": "pandas.DataFrame",
    "pd.Dataframe": "pandas.DataFrame",
    "pd.Series": "pandas.Series",
    "pd.Index": "pandas.Index",
    "nx.Graph": "networkx.Graph",
    "AnnData": "anndata.AnnData",
    "np.array": "numpy.ndarray",
    "string": "str",
    "dictionary": "dict",
    "sets": "set",
}

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_logo = "_static/logo.png"

html_title = project

html_theme_options = {
    "repository_url": repository_url,
    "use_repository_button": True,
    "use_download_button": True,
    "use_fullscreen_button": True,
    "path_to_docs": "docs/",
    "navigation_with_keys": False,
    # Maximum heading level shown in the right-hand page TOC
    "show_toc_level": 4,
    # Navigation levels expanded by default in the left sidebar
    "show_navbar_depth": 2,
}

pygments_style = "default"
katex_prerender = shutil.which(katex.NODEJS_BINARY) is not None

nitpick_ignore = [
    # If building the documentation fails because of a missing link that is outside your control,
    # you can add an exception to this list.
    #     ("py:class", "igraph.Graph"),
]

nitpick_ignore_regex = [
    ("py:class", r"default.*"),             # "default: 0.05", "default 'tab20'", ...
    ("py:class", r"[0-9.eE+-]+"),           # bare numbers left over from defaults
    ("py:class", r"flag"),
    ("py:class", r"pandas\.core\..*"),      # from type annotations
    ("py:class", r"networkx\.classes\..*"),
    ("py:class", r"sns\..*|seaborn\..*"),
    # shorthand names used in docstrings that can't be linked
    ("py:class", r"pd\..*"),
    ("py:class", r"np\..*"),
    ("py:class", r"nx\..*"),
    ("py:class", r"AnnData"),
    ("py:class", r"string|sets|dictionary"),
    # fragments created when napoleon splits "type (default value)"
    ("py:class", r".*\(default"),      # "float (default", "list (default", ...
    ("py:class", r".*\)"),             # "None)", "0.05)", "0.2])", ")"
    ("py:class", r"\[[0-9.]+"),        # "[0.7"
]

suppress_warnings = ["mystnb.unknown_mime_type"]  # see item 5
