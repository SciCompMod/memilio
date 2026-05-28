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

.. _python-package-installation :

Package installation
--------------------

You can run simulations, download data, or create plots with our Python packages.
First, however, you need a working Python installation.
If you are not sure that you have one, head to :ref:`<general-setup-section>`, then continue here.

.. note::

    We highly recommend using a `virtual environment <https://docs.python.org/3/library/venv.html>`__ for installing
    any python package (including ours)! You can find instructions on how to create one
    `here <https://docs.python.org/3/library/venv.html#creating-virtual-environments>`__.
    You can freely choose name and location of the environment, just make sure that you can find it again later.
    And remember to active the environment before trying to use MEmilio!

The Python packages can be installed in two ways, depending on your use case.
We will use pip here, but you can use the python package manager of your choice as well.

Option 1: Install from PyPI
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. info:: Currently only supported for memilio-simulation.

If you just want to run simulations with the latest released version, install the pre-built wheel directly from PyPI:

.. code-block:: console

    python -m pip install memilio-simulation

This is the recommended way to install a package as a user. If you want to develop or contribute to the
project, check out the second option (installing from source).

Pre-built wheels are provided for Linux, Windows and MacOS(ARM).

Option 2: Install from source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you need the latest (unreleased) code, or want to contribute to the package, you need to build from source.
For the :doc:`memilio-simulation <m-simulation>`, :doc:`memilio-surrogatemodel <m-surrogate>` and
:doc:`memilio-generation <m-generation>` packages,
this requires a C++ compiler and CMake, see section :ref:`<general-setup-section>` for setting these up.

Each package provides a ``pyproject.toml`` in the respective directory ``pycode/memilio-{package_name}``, that installs
the package and its python dependencies with pip.

The dependencies of the individual packages are denoted in their documentation.
The installation can be run with the following command from the directory containing the ``pyproject.toml`` file

.. code-block:: console 
    
    cd pycode/memilio-{package_name}
    python -m pip install .

This installs the package and the required dependencies to your site-packages. If you want to avoid changing
directories too often, simply run ``python -m pip install pycode/memilio-{package_name}`` from the project root.

.. tip:: For Contributors: Installing development packages

    The ``-e`` flag installs the package in a mode, which links the installation to your local source code directory.
    This allows working on the python code without having to reinstall the package after every change.
    If you plan to contribute to any package, install all development dependencies by adding ``[dev]``:

    .. code-block:: console

        python -m pip install -e .[dev]

.. warning::
   C++ code changes always require re-running this command for the changes to take effect.

.. dropdown:: :fa:`gears` Expert's knowledge: Build files for skbuild

    The simulaion and generation packages use skbuild to compile python bindings and parts of the C++ library.
    By default, the cmake build files are put into ``pycode/build/memilio-{package_name}`` to save on time during
    package development. If you get unexpected cmake errors, you can try deleting the respective build directory. If
    you do not want to store the build files at all, you can remove the ``build_dir`` entry from the section
    ``[tool.scikit-build]`` in the ``pyproject.toml``. Then skbuild will use a temporary directory instead.

Testing
-------

Each package provides a test suite under ``pycode/memilio-{package_name}/tests``. 
To run the tests, simply use the following command inside the package directory after installation:

.. code-block:: console 

    cd tests
    python -m unittest

This works for both the normal and the editable (with "-e") installations. 

Alternatively, you can start the tests from outside the source directory, but you need to tell the unittest module where
to look. From the MEmilio project directory, run

.. code-block:: console 

    python -m unittest discover -s pycode/memilio-surrogatemodel/tests

Coverage Report
----------------

Dependencies for coverage report:

* coverage

To get the coverage report, run the following in a package directory

.. code-block:: console

    python -m coverage run -m unittest discover -s tests
    python -m coverage report
    python -m coverage xml -o coverage_python.xml
    python -m coverage html -d coverage_python

Inspection via pylint
---------------------

The following packages have to be installed to run pylint:

* pylint
* pylint-json2html

Run pylint with the commands in the package directory

.. code-block:: console

    python ../run_pylint.py
    pylint-json2html -f jsonextended -o build_pylint/pylint.html < build_pylint/pylint_extended.json

From the repository root you can also target a package explicitly, for example
``python pycode/run_pylint.py --package-dir memilio-plot``.
