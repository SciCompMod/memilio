# HPC against Corona #

[![pipeline status](https://gitlab.dlr.de/hpc-against-corona/epidemiology/badges/master/pipeline.svg)](https://gitlab.dlr.de/hpc-against-corona/epidemiology/-/commits/master)
C++: [![c++ coverage report](https://gitlab.dlr.de/hpc-against-corona/epidemiology/badges/master/coverage.svg?job=test-cpp)](https://gitlab.dlr.de/hpc-against-corona/epidemiology/-/jobs/artifacts/master/file/coverage_report/index.html?job=test-cpp)
Python: [![python coverage report](https://gitlab.dlr.de/hpc-against-corona/epidemiology/badges/master/coverage.svg?job=test-py)](https://gitlab.dlr.de/hpc-against-corona/epidemiology/-/jobs/artifacts/master/file/coverage_python/index.html?job=test-py)

This is a common project between the department of Systems Immunology (SIMM) of the he Helmholtz Center for Infection Research (HZI) and the Institute for Software Technology of the German Aerospace Center (DLR). This project will bring cutting edge and compute intensive epidemiological models to a large scale, which enables a precise and high-resolution spatiotemporal pandemic simulation for entire countries.


**Getting started**

This project is divided into multiple building blocks. The implementation of the epidemiological models is to be found in cpp.
Data acquisition tools and data is to be found in data.
The interactive frontend is to be found under frontend.
It is regularly deployed to http://hpcagainstcorona.sc.bs.dlr.de/index.html.
In pycode the epidemiology python package is defined.
It consists of python bindings (seir, secir) to call the C++ code available in cpp and the epidata subpackage which provides tools to download and structure important data.
More about the python package can be found in [Python README](pycode/README.rst).

The C++ code is divided into two libraries, *epidemiology* (models, math, ensemble runs etc.) and *epidemiology_io* (IO routines for parameter sets and simulation results).
By standard, *epidemiology* is built without *epidemiology_io*.


**Documentation**

In every important part of the project, desribed above, a README can be found which contains further information about the structure of the specific project part, the instructions for the user and very important the instructions for the developer.

Furthermore, the code is documented via doxygen and an instruction how to get it, is provided in the docs folder.
The documentation of the code of the master branch can be found here:
https://hpc-against-corona.pages.gitlab.dlr.de/epidemiology/master/documentation/


**Requirements**

By standard *epidemiology* library is bundled with
 * spdlog (https://github.com/gabime/spdlog)
 * eigen v3.3 (http://gitlab.com/libeigen/eigen and http://eigen.tuxfamily.org)
 * boost outcome and optional (https://www.boost.org/)

See thirdparty/CMakeLists.txt for details.

In order to use IO of parameters and simulation results (*epidemiology_io* library), the tool HDF5 (https://www.hdfgroup.org/ e.g., via apt install libhdf5-serial-dev) needs to be installed.

In addition, *epidemiology_io* is bundled with
 * jsoncpp (https://github.com/open-source-parsers/jsoncpp)
 * Boost Filesystem (https://www.boost.org/)

 All bundled libraries can be built as part of this project and don't need to be installed.

 See here for more information on 3rdparty-dependencies: [cpp/thirdparty/Readme.rst](cpp/thirdparty/Readme.rst)

**Installation** 

*Making and executing C++ code*

* (Create a build folder and) do cmake .. (without IO library) or cmake -DEPI_BUILD_EPI_IO=ON .. (with IO library) in epidemiology-cpp/cpp/build
* Do cmake --build . 
* Run 
  * an example via ./examples/secir_ageres
  * all unit tests via ./tests/runUnitTests

*Running code coverage analysis*

Code coverage analysis currently only works on linux with the "Makefile" generator and in Debug mode. To perform
the analysis, configure cmake as follows

    cmake -DCMAKE_BUILD_TYPE=Debug -DEPI_TEST_COVERAGE=ON ..

This step needs to have the tool lcov installed. To execute the coverage, run

    cmake --build . --target coverage

It will generate a html report inside the `coverage` directory.

*Steps to execute C++ code via python bindings*

*  Create a python virtual environment via python3 -m venv virtualenv/ somewhere under the main epidemiology/ folder
*  Activate the environment via source virtualenv/bin/activate
*  Do pip3 install scikit-build
*  In epidemiology/pycode do
   *  python3 setup.py build
   *  python3 setup.py install
   *  execute some example
*  Run the python tests in the pycode/test folder by typing python -m unittest

**Development**
* [Git workflow and change process](https://gitlab.dlr.de/hpc-against-corona/epidemiology/-/wikis/Git-workflow-and-change-process)
* [C++ Coding Guidelines](https://gitlab.dlr.de/hpc-against-corona/epidemiology/-/wikis/Cpp-Coding-Guidlines)
* [Python Coding Guidelines](https://gitlab.dlr.de/hpc-against-corona/epidemiology/-/wikis/Python%20Coding%20Guidelines)
