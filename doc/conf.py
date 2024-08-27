# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import ch_ephem

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "ch_ephem"
copyright = "2013–2024, The CHIME Collaboration"
author = "The CHIME Collaboration"

version = ".".join(ch_ephem.__version__.split(".")[:3])
# The full version, including alpha/beta/rc tags.
release = ch_ephem.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinxcontrib.katex",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "*__*"]

source_suffix = ".rst"

master_doc = "index"

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = False
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_type_aliases = None
napoleon_attr_annotations = True
napoleon_custom_sections = [
    "Searching the index",
    "Getting results",
    "Exceptions",
    "Constants",
    "Classes",
    "Functions",
]

autoclass_content = "both"  # include both class docstring and __init__
autodoc_default_options = {
    # Make sure that autodoc declarations show the right members
    "members": True,
    "show-inheritance": True,
}
autosummary_generate = True  # Make _autosummary files and include them
autosummary_imported_members = False

numpydoc_show_class_members = False


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

pygments_style = "sphinx"

html_theme = "sphinx_rtd_theme"
# html_static_path = ["_static"]
