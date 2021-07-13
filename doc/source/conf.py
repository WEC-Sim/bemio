import sys

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.todo',
    'sphinx.ext.pngmath',
    'sphinx.ext.viewcode',
    'matplotlib.sphinxext.plot_directive',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = 'bemio'
copyright = '2015, National Renewable Energy Laboratory and Sandia National Laboratories'
version = 'v1.0'

# The full version, including alpha/beta/rc tags.
release = 'v1.0a0'

exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'


# -- Options for HTML output ----------------------------------------------
# html_theme = 'classic'
# html_theme = 'basicstrap'
html_theme = 'sphinx_rtd_theme'

html_static_path = ['_static']


htmlhelp_basename = 'bemiodoc'

html_theme_options = {}

# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('bemio_doc', 'bemio', 'bemio Documentation',
     ['Michael Lawson, National Renewable Energy Laboratory'], 1)
]

# Custom paths
sys.path.insert(0,'/Users/mlawson/Applications/bemio_lawsonro3')
