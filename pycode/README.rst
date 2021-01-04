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
Use -h for more information or the README in the corresponding subdirectory.

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


Testing
-------

The following packages have to be installed to run the tests:

- pyfakefs (creates fake directory to test that expected folders are created and data is written)
- freezegun (freezes the time. Thus, the call today can be changed to a specific date.)

To run the tests make 

.. code:: sh

    python -m unittest


Inspection via pylint
---------------------
Pylint report for actial master:


:Pylint Report: https://ssa.pages.gitlab.dlr.de/hpc-against-corona/epidemiology/master/pylint/pylint.html
:Pylint Report: https://ssa.pages.gitlab.dlr.de/hpc-against-corona/epidemiology/269-improve-documentation-of-python-data/pylint/pylint.html



