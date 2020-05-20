Python packages / bindings for the corona project
=================================================

This package collects the python modules for the corona project. Currently, it contains two modules:
 - ``epidemiology.epidata``: Contains scripts to download different kind of data as RKI, John Hopkins, Spain, Population. To do so, install the package (see below) and just call the executable
     
      getrkidata
      getpopuldata
      getjhdata
      getspaindata

 - ``epidemiology.secir``: Contains a python/pybind11 wrapper to the C++ code. An example is provided

More details can be found in the corresponding subdirectories.

Installation
------------

For developement of code use

    python setup.py develop


To install, just enter

    python setup.py install

This builds the C++ extension and copies everything required to you site-packages.



