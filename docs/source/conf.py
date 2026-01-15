# Configuration file for the Sphinx documentation builder.

# -- Project information

import sys
import os
import subprocess
import memilio
project = 'MEmilio'
copyright = '2020-2026 MEmilio'
author = ''

release = ''
version = '1.3.0'


read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

if read_the_docs_build:

    subprocess.call('cd ..; doxygen', shell=True)
    subprocess.call('doxysphinx build . $READTHEDOCS_OUTPUT/html ../Doxyfile')

# sys.path.insert(0, os.path.abspath('../../pycode'))

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinx_copybutton',
    'sphinx_toolbox.collapse',
    'sphinx_design',
    'hoverxref.extension',
    'sphinxcontrib.doxylink',
    #    'sphinx_remove_toctrees'
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

hoverxref_auto_ref = True
hoverxref_roles = ["term"]
hoverxref_domains = ["py"]
hoverxref_role_types = {
    "hoverxref": "tooltip",
    "ref": "tooltip",
    "term": "tooltip",
    "obj": "tooltip",
    "func": "tooltip",
    "mod": "tooltip",
    "meth": "tooltip",
    "class": "tooltip",
}

doxylink = {
    "C++ API": ("source/cppapi/html/tagfile.xml", 
                "source/cppapi/html")
}

# remove_from_toctrees = ["api/*"]

templates_path = ['_templates']
# -- Options for HTML output

html_static_path = ['_static']
html_css_files = [
    'custom.css',
]

html_favicon = "../memilio.ico"

maximum_signature_line_length = 40

html_theme = 'sphinx_rtd_theme'
html_logo = "../memilio-small.png"
html_theme_options = {
    "collapse_navigation": True,
    "logo_only": True,
    "style_nav_header_background": "#f8f9fb",
}

# Mock heavy dependencies to speed up build
autodoc_mock_imports = [
    "numpy",
    "scipy",
    "pandas",
    "matplotlib",
    "tensorflow",
    "scikit-learn",
    "h5py",
    "tables",
    "geopandas",
    "pyarrow",
    "PyQt6",
    "wget",
    "twill",
    "folium",
    "mapclassify",
    "imageio",
    # Mock C++ extension modules that require compilation
    "memilio.simulation",
    "memilio.generation",
]

# -- Options for EPUB output
epub_show_urls = 'footnote'
