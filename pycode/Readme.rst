Python packages / bindings for the corona project
=================================================

This package collects the python modules for the corona project. Currently, it contains two modules:
 - ``epidemiology.epidata``: Contains Kathrins scripts to download the RKI data. To do so, just call the executable
     
      getrkidata

 - ``epidemiology.secir``: Contains a python/pybind11 wrapper to the C++ code. An example is provided

More details can be found in the corresponding subdirectories.

Installation
------------

To install, just enter

    python setup.py install

This builds the C++ extension and copies everything required to you site-packages.



