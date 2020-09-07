Python packages / bindings for the corona project
=================================================

This package collects the python modules for the corona project. Currently, it contains two modules:
 - ``epidemiology.epidata``: Contains scripts to download different kinds of data as RKI, John Hopkins, Spain, Population, DIVI. To do so, install the package (see below) and just call the executable
     
      getrkidata
      getpopuldata
      getjhdata
      getspaindata
      getdividata
      getalldata

 - ``epidemiology.secir``: Contains a python/pybind11 wrapper to the C++ code. An example is provided

More details can be found in the corresponding subdirectories.

Installation
------------

For developement of code use

    python setup.py develop


To install, just enter

    python setup.py install

This builds the C++ extension and copies everything required to you site-packages.


Testing
-------

Following packages have to be installed to run the tests:

- pyfakefs (creates fake directory to test, that expected folders are created and data is written)
- freezegun (freezes the time. Thus, the call today can be changed to a specific date.)

To run the tests make 

    python -m unittest




