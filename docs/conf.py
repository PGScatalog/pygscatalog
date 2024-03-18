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

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

github_username = "pgscatalog"
github_repository = "pygscatalog"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "alabaster"
# html_static_path = ["_static"]

# use autoapi for packages that provide APIs (libraries)
autoapi_dirs = [
    "../pgscatalog.corelib/src/pgscatalog",
    "../pgscatalog.matchlib/src/pgscatalog",
    "../pgscatalog.calclib/src/pgscatalog",
]
# see _templates/autoapi/index.rst for autoapi fix
autoapi_template_dir = "_templates/autoapi"
autoapi_python_use_implicit_namespaces = True
autoapi_keep_files = True

# hide private members
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    "imported-members",
]
autoapi_member_order = "groupwise"


def skip_submodules(app, what, name, obj, skip, options):
    if what == "module":
        skip = True
    return skip


def setup(sphinx):
    sphinx.connect("autoapi-skip-member", skip_submodules)
