Build instructions
==================

The MEmilio core library (MEmilio C++) is written C++ and uses `CMake <https://cmake.org/>`_ as build system. Before
installing MEmilio C++, make sure a C++20 compiler, CMake and a build tool (like GNU Make or Ninja) is installed on your
device. The following guide will make use of the command line, but you can use graphical build tools from an IDE as
well.

Quickstart
^^^^^^^^^^

These are minimal build instructions. More detailed steps and explanations are given below.

.. code:: bash

    git clone https://github.com/SciCompMod/memilio  # download the project
    cd memilio                                       # go into the project directory
    cmake -S cpp -B cpp/build                        # *configure* the project, creating a "build" directory under cpp/
    cmake --build cpp/build -j 2                     # *build* all targets from the project with 2 threads

After the build process is done, you can run files from "cpp/build/bin", for example our test suite:

.. code:: bash

    ./cpp/build/bin/memilio-test

This will run several tests and should write out ``[  PASSED  ]`` in the end.

Requirements
^^^^^^^^^^^^
MEmilio C++ is regularly tested with the following compilers (list will be extended over time):

- GCC, versions 11 and 13
- Clang, version 14 and 17
- MSVC, version 19.43 (Visual Studio 2022)

MEmilio C++ is regularly tested on GitHub runners using Ubuntu 22.04 and 24.04 and Windows Server 2022 and 2025. It is
expected to run on any comparable Linux or Windows system. It is currently not tested on macOS.

The following table lists the dependencies that are used. Most of them are required, but some are optional. The library
can be used without them but with slightly reduced features. CMake will warn about them during configuration. Most of
them are bundled with this library and do not need to be installed manually. Bundled libraries are either included with
this project or loaded from the web on demand. For each dependency, there is a CMake option to use an installed version
instead. Version compatibility needs to be ensured by the user, the version we currently use is included in the table.

.. list-table::
    :header-rows: 1

    * - Library 
      - Version  
      - Required 
      - Bundled               
      - Notes
    * - spdlog  
      - 1.15.0   
      - Yes      
      - Yes (git repo)        
      - https://github.com/gabime/spdlog
    * - Eigen   
      - 3.4.0    
      - Yes      
      - Yes (git repo)        
      - http://gitlab.com/libeigen/eigen
    * - Boost   
      - 1.84.0   
      - Yes      
      - Yes (git repo)        
      - https://github.com/boostorg/boost
    * - JsonCpp 
      - 1.9.6    
      - No       
      - Yes (git repo)        
      - https://github.com/open-source-parsers/jsoncpp
    * - HDF5    
      - 1.12.0   
      - No       
      - No                    
      - https://www.hdfgroup.org/, package libhdf5-dev on apt (Ubuntu)
    * - GoogleTest 
      - 1.10  
      - For Tests only 
      - Yes (git repo)  
      - https://github.com/google/googletest
    * - LibSBML 
      - 5.20.2 
      - No 
      - No 
      - https://sbml.org/software/libsbml/ (For SBML integration only)

See the `thirdparty directory <https://github.com/SciCompMod/memilio/blob/main/cpp/thirdparty/README.md>`_ for more details.

Step-by-step instructions
^^^^^^^^^^^^^^^^^^^^^^^^^

Start by download the newest version (or a specific release) of MEmilio from our
`github repository <https://github.com/SciCompMod/memilio>`_ to a directory of your choice, or clone it directly using
git by first opening a terminal in that directory and then running

.. code:: bash

    git clone https://github.com/SciCompMod/memilio

.. dropdown:: :fa:`gears` Note for developers

    If you need to push changes to the main repo, register an ssh key with GitHub and clone from 
    *git@github.com:SciCompMod/memilio.git* instead. You can also change to the ssh address later using

    .. code:: bash

        git remote set-url origin git@github.com:SciCompMod/memilio.git

This will create a new directory called "memilio". Change into this directory.

.. code:: bash

    cd memilio
    
Before we can *build* anything, we need to *configure* the project first. If you want to use its default options,
simply run

.. code:: bash

    cmake -S cpp -B cpp/build

Additional options can be specified by appending one or more ``-D<OPTION>=<VALUE>``, or by editing the file
``cpp/build/CMakeCache.txt`` after a successful configuration. The following options are known to the library:

.. list-table::
    :header-rows: 1

    * - Option
      - Description
    * - ``MEMILIO_BUILD_TESTS``
      - Build unit tests in the test directory, ON or OFF, default ON.
    * - ``MEMILIO_BUILD_EXAMPLES``
      - Build the example applications in the examples directory, ON or OFF, default ON.
    * - ``MEMILIO_BUILD_MODELS``
      - Build the separate model libraries in the models directory, ON or OFF, default ON.
    * - ``MEMILIO_BUILD_SIMULATIONS``
      - Build the simulation applications in the simulations directory, ON or OFF, default ON.
    * - ``MEMILIO_BUILD_SBML_MODELS``
      - Build the SBML importer and imported models, i.e. everything in the folder ``sbml_model_generation``, ON or OFF, default ON. You may need to set ``sbml_DIR``
    * - ``MEMILIO_USE_BUNDLED_SPDLOG/_BOOST/_EIGEN/_JSONCPP``:
      - Use the corresponding dependency bundled with this project, ON or OFF, default ON.
    * - ``MEMILIO_BUILD_BENCHMARKS``
      - Build the benchmarks for this project, ON or OFF, default OFF.
    * - ``MEMILIO_SANITIZE_ADDRESS/_UNDEFINED``
      - Compile with specified sanitizers to check correctness, ON or OFF, default OFF.
    * - ``MEMILIO_ENABLE_OPENMP``
      - Compile MEmilio with multithreading using OpenMP, ON or OFF, default OFF.
    * - ``MEMILIO_ENABLE_MPI``
      - Compile MEmilio with distributed memory parallelization using MPI. ON or OFF, default OFF. Requires an MPI implementation to be installed on the system. 
    * - ``MEMILIO_ENABLE_WARNINGS``
      - Enable compilation warnings (beyond those enabled in the compiler by default). ON or OFF, default ON.
    * - ``MEMILIO_ENABLE_WARNINGS_AS_ERRORS``
      - Compilation warnings are treated as compilation errors. ON or OFF, default ON.
    * - ``MEMILIO_ENABLE_PROFILING``
      - Compile with runtime profiling support. ON or OFF, default OFF. See `here <https://github.com/SciCompMod/memilio/blob/main/cpp/benchmarks/profiling.md>`_ for information.
    * - ``MEMILIO_ENABLE_LIKWID_MARKER``
      - Compile MEmilio with likwid markers. ON or OFF, default OFF.

Other important options may need:

.. list-table::
    :header-rows: 1

    * - Option
      - Description
    * - ``CMAKE_BUILD_TYPE``
      - Controls compiler optimizations and diagnostics, **Debug**, **Release**, or **RelWithDebInfo**; not available for Multi-Config CMake Generators like Visual Studio, set the build type in the IDE or when running the compiler.
    * - ``CMAKE_INSTALL_PREFIX``
      - Controls the location where the project will be installed
    * - ``HDF5_DIR``
      - If you have HDF5 installed, but it is not found by CMake (usually on the Windows OS), you may have to set this option to the directory in your installation that contains the ``hdf5-config.cmake`` file.
    * - ``sbml_DIR``
      - If you have the SBML library installed, but it is not found by CMake, you may have to set this option to the directory in your installation that contains the ``sbml-config.cmake`` file.

.. note::

    **Example**: You can disable unit tests and enable building benchmarks using 
   
    .. code:: bash
    
        cmake -S cpp -B cpp/build -DMEMILIO_BUILD_TESTS=OFF -DMEMILIO_BUILD_BENCHMARKS=ON

Finally, you can *build* the project by running

.. code:: bash

    cmake --build cpp/build -j <N>

Here, ``<N>`` must be set to the number of jobs used for building MEmilio C++, e.g. the number of available CPU threads
on your system minus two. The argument ``-j <N>`` is optional, but will significantly speed up the compilation.

Once the build command has finished successfully, you can find the compiled binaries in the directory
``cpp/build/bin/``.

You can check that everything is working as intended by running the test suite (if you did not disable it during
configuration):

.. code:: bash

    ./cpp/build/bin/memilio-test

Also try out the example binaries (ending in ``_example``)!

Integration into other projects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using CMake, integration is simple. 

If you installed the project, there is a `memilio-config.cmake` file included with your installation. This config file will tell CMake which libraries and directories have to be included. Look up the config using the command `find_package(memilio)` in your own `CMakeLists.txt`. On Linux, the file should be found automatically if you installed in the normal GNU directories. Otherwise, or if you are working on Windows, you have to specify the `memilio_DIR` variable when running CMake to point it to the `memilio-config.cmake` file. Add the main framework as a dependency with the command `target_link_libraries(<your target> PRIVATE memilio::memilio)`. Other targets that are exported are `memilio::secir`, `memilio::seir`, and `memilio::abm`. This will set all required include directories and libraries, even transitive ones.

Alternatively, `MEmilio` can be integrated as a subdirectory of your project with `add_subdirectory(memilio/cpp)`, then you can use the same `target_link_libraries` command as above.

Installation
~~~~~~~~~~~~

.. note::
    
    **Warning**: Installing currently is not tested and probably does not work as expected or at all. If you want to
    integrate the project into yours, use the `add_subdirectory` way.

After having built MEmilio C++ as described above, you can install it to the location given in the
`CMAKE_INSTALL_PREFIX` variable by running

.. code:: bash

    cmake --install cpp/build

This will install the libraries, headers, and executables that were built, i.e. where `MEMILIO_BUILD_<PART>=ON`.

Known issues
^^^^^^^^^^^^

- Installing currently is not tested and probably does not work as expected or at all. If you want to integrate the project into yours, use the `add_subdirectory` way.
- On Windows, automatic detection of HDF5 installations does not work reliably. If you get HDF5 related errors during the build, you may have to supply the HDF5_DIR variable during CMake configuration, see above.
