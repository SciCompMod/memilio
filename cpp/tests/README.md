# MEmilio Unit Tests #

Unit tests for the MEmilio C++ library using the GoogleTest framework.

*Data handling*

Inputdata that cannot be included in the source code itself is collected in the `data` directory. The path to the directory is discovered by CMake during configuration. That means the compiled unit tests can only be run succesfully if the source directory has not been deleted after building the tests and is not otherwise unreachable. Use [load_test_data.h](load_test_data.h) to access the data.

Some tests write data to the filesystem. The tests try to use the temp directory provided by the operating system. The current directory is used if no temp directory can be found. Files created during tests are deleted afterwards if possible. See [temp_file_register.h](temp_file_register.h)

*Running code coverage analysis*

Code coverage analysis currently only works on linux with the "Makefile" generator and in Debug mode. To perform
the analysis, configure cmake as follows

    cmake -DCMAKE_BUILD_TYPE=Debug -DMEMILIO_TEST_COVERAGE=ON ..

This step needs to have the tool lcov installed. To execute the coverage, run

    cmake --build . --target coverage

It will generate a html report inside the `coverage` directory.

- MEMILIO_TEST_COVERAGE: compile libraries and unit tests for coverage analysis, only available in Debug builds, ON or OFF, default OFF.