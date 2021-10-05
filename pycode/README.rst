Python packages / bindings for the corona project
=================================================

This package collects the python modules for the corona project. Currently, it contains three modules:

* ``epidemiology.epidata``: Contains scripts to download different kinds of data as RKI, John Hopkins, Population, DIVI.

To do so, install the package (see below) and than the following executables are available:
     
      * getrkidata
      * getpopuldata
      * getjhdata
      * getdividata
      * getsimdata
      * cleandata
      * epidemiology/epidata/defaultDict.py
      * getcommutermigration

For all executables there are several run options.
Use -h for more information or the `epidata-readme <epidemiology/epidata/README.rst>`_ in the corresponding subdirectory.

* ``epidemiology.secir``: Contains a python/pybind11 wrapper to the C++ code. An example is provided

More details can be found in the corresponding subdirectories.

Environment
-----------

We recommend to use a virtual environment.
Therefore, do the following.

Create virtiual environment with name "epi_venv" [name can be changed as you want] somewhere into the epidemiology folder(!) 

.. code:: sh

   # Setup virtual environment

   python -m venv epi_venv
   source epi_venv/bin/activate
   pip install --upgrade pip

If the virtual env is not created in the epidemiology-folder the default folder for writing and reading data files of the pidata package is in the site-package folder of the env, see output while writing the data. To avoid this use the -o flag, for details see `epidata-readme <epidemiology/epidata/README.rst>`_  or the `Documentation <https://hpc-against-corona.pages.gitlab.dlr.de/epidemiology/master/documentation/index.html>`_.


Installation
------------

For developement of code use

.. code:: sh

    python setup.py develop


To install, just enter

.. code:: sh

    python setup.py install

This builds the C++ extension and copies everything required to your site-packages. 

All the requirements of the `C++ library <../cpp/README.md>`_ must to be met in order to build and use the python bindings.
CMake is executed internally by the `setup.py` script.

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

    python -m coverage run -m unittest
    python -m coverage report
    python -m coverage xml -o coverage_python.xml
    python -m coverage html -d coverage_python

Coverage report for actual master:

:Coverage Report: https://hpc-against-corona.pages.gitlab.dlr.de/epidemiology/master/coverage/python/index.html


Inspection via pylint
---------------------
The following packages have to be installed to run pylint:

* pylint
* pylint-json2html

After installing the package, run

.. code:: sh

    python setup.py pylint
    pylint-json2html -f jsonextended -o build_pylint/pylint.html < build_pylint/pylint_extended.json

Pylint report for actual master:

:Pylint Report: https://hpc-against-corona.pages.gitlab.dlr.de/epidemiology/master/pylint/pylint.html
