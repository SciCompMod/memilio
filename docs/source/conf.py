# Configuration file for the Sphinx documentation builder.

import sys
import os
import subprocess
import types
import pathlib
import re
import importlib
import pkgutil
import inspect

# -- Project information
project = 'MEmilio'
copyright = '2020-2026 MEmilio'
author = ''
release = ''
version = '1.3.0'

# =============================================================================
# C++ extension stub setup (must run before importing memilio)
# =============================================================================

# All C++ extension modules that need mocking/stubbing
_CPP_EXTENSIONS = [
    "memilio.simulation._simulation",
    "memilio.simulation._simulation_abm",
    "memilio.simulation._simulation_omseirs4",
    "memilio.simulation._simulation_osecir",
    "memilio.simulation._simulation_osecirvvs",
    "memilio.simulation._simulation_oseirdb",
    "memilio.simulation._simulation_oseir",
    "memilio.simulation._simulation_osir",
    "memilio.simulation._simulation_ssir",
    "memilio.simulation._simulation_ssirs",
]

# Mapping from stub files to extension module names
_STUB_MAPPING = {
    '__init__.pyi': 'memilio.simulation._simulation',
    'abm.pyi': 'memilio.simulation._simulation_abm',
    'osecir.pyi': 'memilio.simulation._simulation_osecir',
    'osecirvvs.pyi': 'memilio.simulation._simulation_osecirvvs',
    'oseir.pyi': 'memilio.simulation._simulation_oseir',
    'osir.pyi': 'memilio.simulation._simulation_osir',
}

# Wrapper modules that re-export from extension modules
_WRAPPER_MODULES = {
    'memilio.simulation.osecir': 'memilio.simulation._simulation_osecir',
    'memilio.simulation.osecirvvs': 'memilio.simulation._simulation_osecirvvs',
    'memilio.simulation.oseir': 'memilio.simulation._simulation_oseir',
    'memilio.simulation.osir': 'memilio.simulation._simulation_osir',
    'memilio.simulation.abm': 'memilio.simulation._simulation_abm',
}


def _create_module_from_stub(stub_path, module_name):
    """Create a module with classes/functions extracted from a .pyi stub file."""
    content = stub_path.read_text()
    module = types.ModuleType(module_name)
    module.__file__ = str(stub_path)
    module.__doc__ = f"C++ bindings (from stub: {stub_path.name})"
    exported_names = []

    # Extract class definitions with methods
    class_pattern = re.compile(
        r'^class\s+(\w+)(?:\([^)]*\))?:\s*\n((?:[ \t]+.*\n)*)',
        re.MULTILINE
    )
    for match in class_pattern.finditer(content):
        class_name = match.group(1)
        class_body = match.group(2) or ''

        # Extract docstring
        doc_match = re.search(r'^\s+"""(.*?)"""', class_body, re.DOTALL)
        docstring = doc_match.group(1).strip() if doc_match else f'{class_name} - C++ binding class'

        # Extract methods
        methods = {}
        for m in re.finditer(r'^\s+def\s+(\w+)\s*\([^)]*\)[^:]*:', class_body, re.MULTILINE):
            name = m.group(1)
            methods[name] = lambda self, *a, **k: None
            methods[name].__name__ = name
            methods[name].__doc__ = f'{name} method'

        placeholder = type(class_name, (), {'__doc__': docstring, '__module__': module_name, **methods})
        setattr(module, class_name, placeholder)
        if not class_name.startswith('_'):
            exported_names.append(class_name)

    # Extract top-level functions
    for match in re.finditer(r'^def\s+(\w+)\s*\([^)]*\)[^:]*:', content, re.MULTILINE):
        func_name = match.group(1)
        if not func_name.startswith('_'):
            def make_func(name):
                def f(*args, **kwargs): pass
                f.__name__ = name
                f.__doc__ = f'{name} - C++ binding function'
                return f
            setattr(module, func_name, make_func(func_name))
            exported_names.append(func_name)

    module.__all__ = exported_names
    return module


def _setup_cpp_stubs():
    """Set up stub modules for C++ extensions when real bindings unavailable."""
    try:
        import memilio.simulation._simulation
        return  # Real bindings available
    except ImportError:
        pass

    stubs_dir = pathlib.Path(__file__).parent.parent.parent / 'pycode' / 'memilio-simulation-stubs' / 'memilio-stubs' / 'simulation'

    for ext_module in _CPP_EXTENSIONS:
        # Find matching stub file
        stub_file = None
        for stub_name, mod_name in _STUB_MAPPING.items():
            if mod_name == ext_module:
                stub_file = stubs_dir / stub_name
                break

        if stub_file and stub_file.exists():
            try:
                sys.modules[ext_module] = _create_module_from_stub(stub_file, ext_module)
            except Exception:
                sys.modules[ext_module] = types.ModuleType(ext_module)
                sys.modules[ext_module].__all__ = []
        else:
            sys.modules[ext_module] = types.ModuleType(ext_module)
            sys.modules[ext_module].__all__ = []


def _register_submodules():
    """Import and register submodules on their parent packages for autodoc."""
    for pkg_name in ['memilio.epidata', 'memilio.plot', 'memilio.simulation', 'memilio.generation']:
        try:
            pkg = importlib.import_module(pkg_name)
            pkg_path = getattr(pkg, '__path__', None)
            if not pkg_path:
                continue
            for _, modname, _ in pkgutil.iter_modules(pkg_path):
                if not modname.startswith('_'):
                    try:
                        setattr(pkg, modname, importlib.import_module(f'{pkg_name}.{modname}'))
                    except Exception:
                        pass
        except Exception:
            pass


def _fix_module_attributes():
    """Fix __module__ attributes so autodoc includes classes from wrapper modules."""
    for wrapper_name, ext_name in _WRAPPER_MODULES.items():
        try:
            wrapper_mod = importlib.import_module(wrapper_name)
        except Exception:
            continue

        for name in dir(wrapper_mod):
            if name.startswith('_'):
                continue
            obj = getattr(wrapper_mod, name, None)
            if obj and (inspect.isclass(obj) or inspect.isfunction(obj)):
                if getattr(obj, '__module__', None) == ext_name:
                    try:
                        obj.__module__ = wrapper_name
                    except (AttributeError, TypeError):
                        pass


# Run setup functions
_setup_cpp_stubs()
import memilio
_register_submodules()
_fix_module_attributes()

# =============================================================================
# Sphinx configuration
# =============================================================================

if os.environ.get('READTHEDOCS', None) == 'True':
    subprocess.call('cd ..; doxygen', shell=True)

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
    'breathe',
    'exhale',
    'hoverxref.extension',
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
    "doxygenStripFromPath":"..",
    "rootFileTitle":       "C++ API",
    "createTreeView":      True,
    "treeViewIsBootstrap": False,
    "contentsDirectives":  False,
}
breathe_projects = {"MEmilio": "../xml"}
breathe_default_project = "MEmilio"

# remove_from_toctrees = ["api/*"]

templates_path = ['_templates']

html_static_path = ['_static']
html_css_files = ['custom.css']
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
    "magic",
    "requests",
    "urllib3",
    "progress",
    "clang",
    "graphviz",
    "dataclasses_json",
]

# EPUB output
epub_show_urls = 'footnote'
