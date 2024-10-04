# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "pygscatalog"
copyright = "2024, PGS Catalog"
author = "PGS Catalog"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "autoapi.extension",
    "sphinx_toolbox.sidebar_links",
    "sphinx_toolbox.github",
    "sphinx_toolbox.shields",
]
extensions.append("autoapi.extension")

autoapi_dirs = [
    "../pgscatalog.core/src/",
    "../pgscatalog.calc/src/",
    "../pgscatalog.match/src/",
]
autoapi_keep_files = True
autoapi_options = ["members", "undoc-members", "show-module-summary"]
autoapi_member_order = "groupwise"
# don't document CLI programs - they have separate documentation
autoapi_ignore = ["*cli*"]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", ".venv"]

github_username = "pgscatalog"
github_repository = "pygscatalog"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "alabaster"
# html_static_path = ["_static"]
