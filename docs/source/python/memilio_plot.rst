MEmilio Plot
=============

MEmilio Plot provides modules and scripts to plot epidemiological or simulation data as returned
by other packages of the MEmilio software.

The package is contained inside the folder `pycode/memilio-plot <https://github.com/SciCompMod/memilio/blob/main/pycode/memilio-plot>`_.


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
- imageio
- datetime

Testing and Coverage
--------------------

The following packages are used by the tests:

* pyfakefs (creates fake directory to test that expected folders are created and data is written)
* coverage

See Installation on how to install all these dependencies automatically.

To run the tests make 

.. code-block:: console

    python -m unittest

To get the coverage report do

.. code-block:: console

    python -m coverage run -m unittest
    python -m coverage report
    python -m coverage xml -o coverage_python.xml
    python -m coverage html -d coverage_python

Coverage report for actual master:

`Coverage Report <https://scicompmod.github.io/memilio/coverage/python/>`_

Inspection via pylint
---------------------

The following packages have to be installed to run pylint:

* pylint
* pylint-json2html

See Installation on how to install all these dependencies automatically.

Run pylint with the commands

.. code-block:: console

    python setup.py pylint
    pylint-json2html -f jsonextended -o build_pylint/pylint.html < build_pylint/pylint_extended.json

Pylint report for actual master:

`Pylint Report <https://dlr-sc.github.io/memilio/pylint/>`_

