# Epidemiology #

[![pipeline status](https://gitlab.dlr.de/hpc-against-corona/epidemiology/badges/master/pipeline.svg)](https://gitlab.dlr.de/hpc-against-corona/epidemiology/-/commits/master)

This is a common project between the department of Systems Immunology (SIMM) of the he Helmholtz Center for Infection Research (HZI) and the Institute for Software Technology of the German Aerospace Center (DLR). This project will bring cutting edge and compute intensive epidemiological models to a large scale, which enables a precise and high-resolution spatiotemporal pandemic simulation for entire countries.


**Getting started**

This project is divided into multiple building blocks. The implementation of the epidemiological models is to be found in cpp. Data acquisition tools and data is to be found in data. The interactive frontend is to be found under frontend. It is regularly deployed to http://hpcagainstcorona.sc.bs.dlr.de/index.html. In pycode you find python bindings to call the C++ code available in cpp. At the moment, some data tools are still under pycode, too.

The C++ code is divided into two libraries, *epidemiology* (models, math, ensemble runs etc.) and *epidemiology_io* (IO routines for parameter sets and simulation results). By standard, *epidemiology* is built without *epidemiology_io*.


**Requirements**

By standard *epidemiology* library is bundled with
 * spdlog (https://github.com/gabime/spdlog)
 * eigen v3.3 (http://gitlab.com/libeigen/eigen and http://eigen.tuxfamily.org)

See thirdparty/CMakeLists.txt for details.

In order to use IO of parameters and simulation results (*epidemiology_io* library), the tools
  * tixi3 (https://github.com/DLR-SC/tixi) and 
  * hdf5 (https://www.hdfgroup.org/ e.g., via apt install libhdf5-serial-dev)
  * jsoncpp (via apt install libjsoncpp-dev)
  need to be installed.


**Installation** 

*Making and executing C++ code*

* (Create a build folder and) do cmake .. (without IO library) or cmake -DEPI_BUILD_EPI_IO=ON .. (with IO library) in epidemiology-cpp/cpp/build
* Do cmake --build . 
* Run 
  * an example via ./examples/secir_ageres
  * all unit tests via ./tests/runUnitTests

*Steps to execute C++ code via python bindings*

*  Create a python virtual environment via python3 -m venv virtualenv/
*  Activate the environment via source virtualenv/bin/activate
*  Do pip3 install scikit-build
*  In epidemiology-cpp/pycode do
   *  python3 setup.py build
   *  python3 setup.py install
   *  execute some example
*  Run the python tests in the pycode/test folder by typing python -m unittest

**Development**
* [Git workflow and change process](https://gitlab.dlr.de/hpc-against-corona/epidemiology/-/wikis/Git-workflow-and-change-process)
* [C++ Coding Guidelines](https://gitlab.dlr.de/hpc-against-corona/epidemiology/-/wikis/Cpp-Coding-Guidlines)
* [Python Coding Guidelines](https://gitlab.dlr.de/hpc-against-corona/epidemiology/-/wikis/Python%20Coding%20Guidelines)
