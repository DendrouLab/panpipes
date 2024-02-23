# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'Panpipes'
copyright = '2023, Rich-Griffin & Curion'
author = 'Charlotte Rich-Griffin & Fabiola Curion'

release = '0.4'
version = '0.4.1'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
        'myst_parser'
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

source_suffix = ['.rst', '.md'] 
# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = "img/panpipeslogo2.png"
html_theme_options = {
    'logo_only': True,
    'display_version': False
}


# -- Options for EPUB output
epub_show_urls = 'footnote'
