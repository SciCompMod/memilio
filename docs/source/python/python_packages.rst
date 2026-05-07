Overview
=========

MEmilio contains several python modules offering an easy-to-use interface to the efficiently implemented C++
models and to complement these by data preparation, code generation, machine learning, or plotting functionality.
Please see the individual package documentation for more details on the functionality and usage.

.. grid:: 1 1 2 2
    :gutter: 2 3 4 4

    .. grid-item-card::
        :img-top: http://martinkuehn.eu/research/images/pybind.png
        :text-align: center

        MEmilio Python Interface
        ^^^

        This package provides a python interface for parts of the C++ main library,
        with the goal of exposing fast mathematical-epidemiological models to
        a bigger user base.

        +++

        .. button-ref:: m-simulation
            :expand:
            :color: secondary
            :click-parent:

            To the python bindings

    .. grid-item-card::
        :img-top: http://martinkuehn.eu/research/images/epidata.png
        :text-align: center

        MEmilio EpiData Package
        ^^^

        This package provides tools to download and structure important 
        data such as infection or derived mobility data.

        +++

        .. button-ref:: m-epidata
            :expand:
            :color: secondary
            :click-parent:

            To the data package


.. grid:: 1 1 1 1
    :gutter: 2 3 4 4

    .. grid-item-card::
        :img-top: http://martinkuehn.eu/research/images/surrogate.png
        :text-align: center

        Machine Learning Surrogate Models
        ^^^

        This package contains machine learning-based surrogate models that were trained based on the MEmilio simulation outputs. 

        +++

        .. button-ref:: m-surrogate
            :expand:
            :color: secondary
            :click-parent:

            To the intro of surrogate models

    .. grid-item-card::
        :img-top: http://martinkuehn.eu/research/images/plot.png
        :text-align: center

        Visualization
        ^^^

        Generalized visualization functions for MEmilio specific plots.

        +++

        .. button-ref:: m-plot
            :expand:
            :color: secondary
            :click-parent:

            To the visualization
   
    .. grid-item-card::
        :img-top: http://martinkuehn.eu/research/images/pygen.png
        :text-align: center

        Interface Generation
        ^^^

        Easy to use tool for helping with the creation of python bindings or interfaces to (new) C++ models.

        +++

        .. button-ref:: m-generation
            :expand:
            :color: secondary
            :click-parent:

            To the generation package


.. _Python_Installation:



Installation
------------

The Python packages can be installed in two ways depending on your use case.

**Option 1: Install from PyPI (Recommended - no C++ compiler required; currently only supported for memilio-simulation)**

Pre-built wheels are provided for Linux and Windows on Python 3.8 to 3.13.

.. code-block:: console

   pip install memilio-simulation

This is the easiest way to get started. No C++ compiler or CMake is needed.

**Option 2: Install from source (latest development version, or contributing)**

Each package provides a ``pyproject.toml`` that installs the package and its dependencies with pip.

If you want to install the memilio-simulation package: 

* The pyproject.toml is in the **root of the MEmilio repository**.

If you want to install any of the other Python packages:

* The pyproject.toml is in the respective folder ``pycode/memilio-*``

The dependencies of the individual packages are denoted in their documentation.
The installation can be run with the following command from the directory containing the ``pyproject.toml`` file

.. code-block:: console 
    
    python -m pip install .

This copies the package and the required dependencies to your site-packages.

For development of code use this command instead

.. code-block:: console 
    
    python -m pip install -e .[dev]

This command allows you to work on the code without having . It also installs 
all additional dependencies required for development and maintenance.


The ``-e`` flag links the installation to your local source code so Python changes are reflected immediately. Hence, 
you do not need to reinstall the package after a changes in Python. C++ changes require re-running this command to 
recompile.

.. dropdown:: :fa:`gears` Build files for skbuild

    The simulaion and generation packages use skbuild to compile python bindings or parts of the C++ library.
    By default, the cmake build files are put into ``pycode/build/memilio-{package_name}`` to save on time during
    package development. If you get unexpected cmake errors, you can try and delete the respective build directory. If
    you do not want to store the build files at all, you can remove the ``build_dir`` entry from the section
    ``[tool.scikit-build]`` in the ``pyproject.toml``. Then skbuild will use a temporary directory instead.

Testing
-------

Each package provides a test suite under ``pycode/memilio-{package_name}/tests``. 
To run the tests, simply use the following command inside the package folder after installation:

.. code-block:: console 

    cd tests
    python -m unittest

Coverage Report
----------------

Dependencies for coverage report:

* coverage

To get the coverage report do in the package folder

.. code-block:: console

    python -m coverage run -m unittest
    python -m coverage report
    python -m coverage xml -o coverage_python.xml
    python -m coverage html -d coverage_python

Coverage report for actual master:

`Coverage Report <https://scicompmod.github.io/memilio/coverage/python/>`__

Inspection via pylint
---------------------

The following packages have to be installed to run pylint:

* pylint
* pylint-json2html

Run pylint with the commands in the package folder

.. code-block:: console

    python ../run_pylint.py
    pylint-json2html -f jsonextended -o build_pylint/pylint.html < build_pylint/pylint_extended.json

From the repository root you can also target a package explicitly, for example
``python pycode/run_pylint.py --package-dir memilio-plot``.
