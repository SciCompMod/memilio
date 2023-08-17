MEmilio Plot Package
=======================

Content
-------

- Introduction
- Installation
- Dependencies
- Running the scripts
- Testing and Coverage
- Pylint
- Additional tools
- Notes for developers (!)
- Troubleshooting

Introduction
------------

This package provides modules and scripts to plot epidemiological or simulation data as returned
by other packages of the MEmilio software.

Installation
------------

Use the provided ``setup.py`` script to install the package and its dependencies.

To install the package, use (from the directory that contains ``setup.py``)

.. code:: sh

    pip install .

This copies everything required to your site-packages.

For developement of code use the command 

.. code:: sh

    pip install -e .[dev]

This command allows you to work on the code without having to reinstall the package after a change. It also installs all additional dependencies required for development and maintenance.

Dependencies
------------

Required python packages:

- pandas>=1.2.2
- matplotlib
- numpy>=1.22,<1.25
- openpyxl
- xlrd
- requests
- pyxlsb
- wget
- folium
- matplotlib
- mapclassify
- geopandas
- h5py

Testing and Coverage
--------------------

The following packages are used by the tests:

- pyfakefs (creates fake directory to test that expected folders are created and data is written)
- coverage

See Installation on how to install all these dependencies automatically.

To run the tests make 

.. code:: sh

    python -m unittest

To get the coverage report do

.. code:: sh

    python -m coverage run -m unittest
    python -m coverage report
    python -m coverage xml -o coverage_python.xml
    python -m coverage html -d coverage_python

Coverage report for actual master:

:Coverage Report: https://dlr-sc.github.io/memilio/coverage/python/

Inspection via pylint
---------------------
The following packages have to be installed to run pylint:

* pylint
* pylint-json2html

See Installation on how to install all these dependencies automatically.

Run pylint with the commands

.. code:: sh

    python setup.py pylint
    pylint-json2html -f jsonextended -o build_pylint/pylint.html < build_pylint/pylint_extended.json

Pylint report for actual master:

:Pylint Report: https://dlr-sc.github.io/memilio/pylint/

Notes for developers
--------------------

If a new functionality shall be added please stick to the following instructions:

When you start creating a new script:

- use doxygen like comments in code as
    - add description in the beginning of the file
        - ## Header
        - # @brief name descr
        - # longer description
    - add description in the beginning of every function directly after the definition
        - start and end with """
        - add a short description to first line
        - afterwards add a longer description
        - # @param name of parameter
        - # @return type description

General
- Always add unittests
- Check test coverage report, if every new feature is covered.
- Check the pylint report just comments with "refactor" are allowed.
