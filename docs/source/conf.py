# Configuration file for the Sphinx documentation builder.

# -- Project information

import sys
import os
import subprocess
import memilio
project = 'MEmilio'
copyright = '2020-2025 MEmilio'
author = ''

release = ''
version = '1.3.0'


read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

if read_the_docs_build:

    subprocess.call('cd ..; doxygen', shell=True)

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
    'sphinx_design',
    'breathe',
    'exhale',
    'hoverxref.extension',
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

exhale_args = {
    "containmentFolder":   "./api",
    "rootFileName":        "library_root.rst",
    "doxygenStripFromPath":    "..",
    "rootFileTitle":       "C++ API",
    "createTreeView":      True,
    "treeViewIsBootstrap": False,
    "contentsDirectives":    False,
}

breathe_projects = {"MEmilio": "../xml"}
breathe_default_project = "MEmilio"

# remove_from_toctrees = ["api/*"]

templates_path = ['_templates']
# -- Options for HTML output

html_static_path = ['_static']
html_css_files = [
    'custom.css',
]
# html_js_files = [
#     'custom.js',
# ]

maximum_signature_line_length = 40

html_theme = 'sphinx_rtd_theme'
html_logo = "../memilio-small.png"
html_theme_options = {
    # "sidebar_hide_name": True,
    # "footer_icons": [
    #     {
    #         "name": "GitHub",
    #         "url": "https://github.com/SciCompMod/memilio",
    #         "html": """
    #             <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
    #                 <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
    #             </svg>
    #         """,
    #         "class": "",
    #     },
    # ],
    "collapse_navigation": False,
    "logo_only": True,
    "style_nav_header_background": "#f8f9fb",
}
# def setup(app):
#     app.add_css_file('custom.css')  # Use app.add_stylesheet in older Sphinx versions
#     app.add_js_file('custom.js')

# -- Options for EPUB output
epub_show_urls = 'footnote'
