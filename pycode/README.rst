Python packages / bindings for the corona project
=================================================

This package collects the python modules for the corona project. Currently, it contains three modules:
 - ``epidemiology.epidata``: Contains scripts to download different kinds of data as RKI, John Hopkins, Spain, Population, DIVI.
To do so, install the package (see below) and than the following executables are available:
     
      - getrkidata
      - getpopuldata
      - getjhdata
      - getspaindata
      - getdividata
      - getalldata

For all executables there are several run options.
Use -h for more information or the :ref:'README <epidata/README>' in the corresponding subdirectory.

 - ``epidemiology.secir``: Contains a python/pybind11 wrapper to the C++ code. An example is provided

 - ``epidemiology.seir``: Contains a python/pybind11 wrapper to the C++ code. An example is provided

More details can be found in the corresponding subdirectories.

Environment
-----------

We recommend to use a virtual environment.
Therefore, do the following.

Create virtiual environment with name "epi_venv" [name can be changed as you want]
.. code:: sh
   # Setup virtual environment

   python3 -m venv epi_venv
   source epi_venv/bin/activate
   pip install --upgrade pip


Installation
------------

For developement of code use

.. code:: sh

    python setup.py develop


To install, just enter

.. code:: sh

    python setup.py install

This builds the C++ extension and copies everything required to your site-packages.


Testing and Coverage
--------------------

The following packages have to be installed to run the tests:

- pyfakefs (creates fake directory to test that expected folders are created and data is written)
- freezegun (freezes the time. Thus, the call today can be changed to a specific date.)

To run the tests make 

.. code:: sh

    python -m unittest

To get the coverage report do
.. code:: sh

    python3 -m coverage report
    python3 -m coverage xml -o coverage_python.xml
    python3 -m coverage html -d coverage_python


Inspection via pylint
---------------------
The following packages have to be installed to run pylint:

- pylint
- pylint-json2html

After installing the package, run

.. code:: sh

    python3.6 setup.py pylint
    pylint-json2html -f jsonextended -o build_pylint/pylint.html < build_pylint/pylint_extended.json

Pylint report for actual master:

:Pylint Report: https://hpc-against-corona.pages.gitlab.dlr.de/epidemiology/master/pylint/pylint.html
:Pylint Report: https://hpc-against-corona.pages.gitlab.dlr.de/epidemiology/269-improve-documentation-of-python-data/pylint/pylint.html
