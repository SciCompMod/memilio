Build instructions
==================

The MEmilio core library (MEmilio C++) is written in C++ and uses `CMake <https://cmake.org/>`__ as build system. For
building MEmilio C++, you need a C++20 compiler, CMake and a build tool (like GNU Make or Ninja) installed on your
device. Follow the steps from :doc:`../download_and_setup` to set up the required tools and the project.

The following guide will make use of the command line, but you can use graphical build tools from an IDE as well.

Quick start
-----------

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
------------

MEmilio C++ is regularly tested with the following compilers (list will be extended over time):

- GCC, versions 11 and 13
- Clang, version 14 and 17
- AppleClang, version 21
- MSVC, versions 19.44+ (Visual Studio 2022) and 19.51+ (Visual Studio 2026)

MEmilio C++ is regularly tested on GitHub runners using Ubuntu 22.04 and 24.04, MacOS 26 and Windows Server 2022
and 2025. It is expected to run on any comparable Linux, MacOS or Windows system.

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
      - 1.17.0   
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

See the `thirdparty directory <https://github.com/SciCompMod/memilio/blob/main/cpp/thirdparty/README.md>`__ for more details.

Step-by-step instructions
-------------------------

In case you have not installed the required tools, or downloaded the MEmilio repository, first follow the steps from the
:doc:`../download_and_setup` section. We assume that you have a terminal open and have navigated into the the memilio
repository.

Configuration
~~~~~~~~~~~~~

Before we can *build* anything, we need to *configure* the project first. If you want to use its default options,
simply run cmake in MEmilio's project directory.

.. code:: bash

    cmake -S cpp -B cpp/build

The first time you run this it may take some time to complete, as this downloads all bundled dependencies, if not
disabled by the options below. Depending on your internet connection, this may take a couple minutes.

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
    * - ``MEMILIO_ENABLE_HDF5``
      - Build MEmilio with HDF5 IO support, ON or OFF, default ON. If OFF, HDF5 is not searched for and features like ``save_result``/``read_result`` are disabled.
    * - ``MEMILIO_ENABLE_SBML``
      - Build the SBML importer and imported models, i.e. everything in the folder ``sbml_model_generation``, ON or OFF, default OFF. You may need to set ``sbml_DIR``
    * - ``MEMILIO_ENABLE_WARNINGS``
      - Enable compilation warnings (beyond those enabled in the compiler by default). ON or OFF, default ON.
    * - ``MEMILIO_ENABLE_WARNINGS_AS_ERRORS``
      - Compilation warnings are treated as compilation errors. ON or OFF, default ON.
    * - ``MEMILIO_ENABLE_PROFILING``
      - Compile with runtime profiling support. ON or OFF, default OFF. See `here <https://github.com/SciCompMod/memilio/blob/main/cpp/benchmarks/profiling.md>`__ for information.
    * - ``MEMILIO_ENABLE_LIKWID_MARKER``
      - Compile MEmilio with likwid markers. ON or OFF, default OFF.

Other important options you may need:

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

.. tip::

    **Example**: You can disable unit tests and enable building benchmarks using 
   
    .. code:: bash
    
        cmake -S cpp -B cpp/build -DMEMILIO_BUILD_TESTS=OFF -DMEMILIO_BUILD_BENCHMARKS=ON

Building
~~~~~~~~

Finally, you can *build* the project by running

.. code:: bash

    cmake --build cpp/build -j 4

Here, ``-j 4`` is the number of jobs used for building MEmilio C++. Optimally, instead of using "4", you should set it
to the number of available CPU threads on your system, minus two. But "4" should work well on all modern systems.
The argument ``-j 4`` *is* optional, but will significantly speed up the compilation.

Once the build command has finished successfully, you can find the compiled executables in the directory
``cpp/build/bin/``.

You can check that everything is working as intended by running the test suite (if you did not disable the tests during
configuration):

.. code:: bash

    ./cpp/build/bin/memilio-test

Also try out the example binaries (ending in ``_example``)!

If you want to only build a specific example, you can specify it with the ``--target`` flag:

.. code-block:: console

   cmake --build cpp/build --target <example_name>

If you experience errors, feel free to contact martin.kuehn@dlr.de or open a
`discussion on GitHub <https://github.com/SciCompMod/memilio/discussions>`__!

Integration into other projects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using CMake, integration is simple.

MEmilio can be integrated as a subdirectory of your project with FetchContent:

.. code:: cmake

  include(FetchContent)
  
  FetchContent_Declare(
    memilio
    GIT_REPOSITORY https://github.com/SciCompMod/memilio.git
    GIT_TAG <branch, release or commit hash>
    SOURCE_SUBDIR cpp
  )

  # Set some options for the build.
  set(MEMILIO_BUILD_TESTS OFF)
  set(MEMILIO_BUILD_EXAMPLES OFF)

  # Include MEmilio directly
  FetchContent_MakeAvailable(memilio)

Set the ``GIT_TAG`` to a branch, release, or commit hash. The "main" branch will always work, but we strongly recommend
using a release version, preferably as a hash. 
Then, add the main framework as a dependency with the command ``target_link_libraries(<your target> PRIVATE memilio)``. You may also want to add one or more models as a dependency, like ``ode_secir``, ``ide_seir``, and ``abm``.
This will set all required include directories and libraries, even transitive ones.

Alternatively, if you installed the project (see the next section), there is a `memilio-config.cmake` file included with
your installation. This config file will tell CMake which libraries and directories have to be included. Look up the
config using the command ``find_package(memilio)`` in your own `CMakeLists.txt`. On Linux, the file should be found
automatically if you installed it in the normal GNU directories. Otherwise, or if you are working on Windows, you have
to specify the ``memilio_DIR`` variable when running CMake to point it to the `memilio-config.cmake` file.
Then you can add dependencies as above, but with the added prefix of ``memilio::``. For example, to add the main
framework, use ``target_link_libraries(<your target> PRIVATE memilio::memilio)``.

Installation
~~~~~~~~~~~~

.. warning::
    
    Installing currently is not tested and probably does not work as expected or at all. If you want to
    integrate the project into yours, use the `FetchContent` way.

After having built MEmilio C++ as described above, you can install it to the location given in the
`CMAKE_INSTALL_PREFIX` variable by running

.. code:: bash

    cmake --install cpp/build

This will install the libraries, headers, and executables that were built, i.e. where ``MEMILIO_BUILD_<PART>=ON``.

Next steps
----------

Running simulations
~~~~~~~~~~~~~~~~~~~

You can run simulations either via the C++ interface where they are originally implemented or via the python bindings. 
For the C++ Interface, you can find explanations of the models as well as guides on their usage in the :doc:`C++ model usage <model_usage>` section.
In short, the executables for different model instantiations are built as described above and can be run via 

.. code-block:: console

   ./cpp/build/bin/<example_name>


Out of the box this works for all examples in the ``cpp/examples`` folder of our
`github repository <https://github.com/SciCompMod/memilio/tree/main/cpp/examples>`__,
that do not depend on user-provided external libraries. 
Additional explanations for our models are linked at the corresponding sites of this documentation.

Simulations used in publications
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For simulations used in publications, we maintain a separate repository: 
`memilio-simulations <https://github.com/SciCompMod/memilio-simulations>`__. 
This repository contains simulations organized in separate folders, each with the specific version of MEmilio 
used for the published results. This ensures that simulation results can be easily reproduced.

The repository also includes additional scripts for plotting, data gathering, and pre-/post-processing 
that were used in publications.

Loading data
~~~~~~~~~~~~

The :doc:`memilio-epidata <../python/m-epidata>` package provides tools to download epidemiological relevant datasets.
Some datasets like contact matrices for Germany are also included in the ``data`` folder of the
`github repository <https://github.com/SciCompMod/memilio/tree/main/data>`__ and school holidays (for Germany) are
directly included in the
`C++ code <https://github.com/SciCompMod/memilio/blob/main/cpp/memilio/geography/holiday_data.ipp>`__.  


Creating new models
~~~~~~~~~~~~~~~~~~~

If you want to create new models, you can do so via the C++ interface. For this, we recommend to have a look at 
the :doc:`C++ model creation <model_creation>` section of this documentation.


Visualizations
~~~~~~~~~~~~~~

For visualizations, we provide our :doc:`python package memilio-plot <../python/m-plot>`. Apart from that, we have 
collected some scripts that we used for visualizations in the `tools folder in our github repository <https://github.com/SciCompMod/memilio/tree/main/tools>`__. 
For the latter, no regular testing is conducted. If you encounter errors, please `contact us <mailto:Martin.Kuehn@DLR.de>`__.


Known issues
------------

- Installing currently is not tested and probably does not work as expected or at all. If you want to integrate the project into yours, use the `FetchContent` way.
- On Windows, automatic detection of HDF5 installations does not work reliably. If you get HDF5 related errors during the build, you may have to supply the HDF5_DIR variable during CMake configuration, see above.

Further questions
-----------------

If you have any further questions, please take a look at our :doc:`../faq` and feel free to contact us via
`e-mail <mailto:Martin.Kuehn@DLR.de>`__ or open an issue or discussion on
`GitHub <https://github.com/SciCompMod/memilio/discussions>`__.
