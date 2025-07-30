Configuration Macros
-------------------

MEmilio provides several compile-time configuration options through CMake that control which features and dependencies are available during compilation. These configurations are defined as preprocessor macros and can be used to conditionally compile code sections based on available libraries and enabled features.
The following macros are automatically defined by CMake during the build process based on the availability of dependencies and build options:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Macro Name
     - Description
   * - :code:`MEMILIO_HAS_HDF5`
     - Defined when HDF5 library is available. Enables HDF5 file I/O functionality for reading and writing results.
   * - :code:`MEMILIO_HAS_JSONCPP`
     - Defined when JsonCpp library is available. Enables JSON file I/O functionality for parameter import and configuration files.
   * - :code:`MEMILIO_ENABLE_MPI`
     - Defined when MPI (Message Passing Interface) support is enabled. Allows distributed computing and parallel simulations across multiple processes or compute nodes.
   * - :code:`MEMILIO_ENABLE_OPENMP`
     - Defined when OpenMP support is enabled. Enables shared-memory parallelization for multi-threaded execution within a single process.
   * - :code:`MEMILIO_ENABLE_PROFILING`
     - Defined when profiling support is enabled. Activates performance monitoring and timing instrumentation throughout the codebase.

Usage in Code
~~~~~~~~~~~~~

These macros are used with preprocessor conditionals to enable or disable specific functionality. You should use these macros to ensure that your code can compile and run correctly depending on the available libraries and features. 
When using these macros, you should always consider fallbacks or alternative implementations when a feature is not available. This ensures that your code remains robust and can handle cases where optional dependencies are not present.
Here are common usage patterns:

**Conditional Compilation for Optional Dependencies**

.. code-block:: cpp

    #ifdef MEMILIO_HAS_JSONCPP
    #include "memilio/io/epi_data.h"
    #include "memilio/io/result_io.h"
    
    // JSON-based parameter I/O functions
    template <typename FP = double>
    IOResult<void> read_divi_data(const std::string& path, 
                                  const std::vector<int>& vregion, 
                                  Date date,
                                  std::vector<FP>& vnum_icu) {
    }
    #endif // MEMILIO_HAS_JSONCPP

**HDF5-Specific Functionality**

.. code-block:: cpp

    #ifdef MEMILIO_HAS_HDF5
    template <class Model>
    IOResult<void> export_input_data_county_timeseries(
        std::vector<Model> models, 
        const std::string& results_dir,
        /* ... other parameters ... */) {
        // HDF5-based time series export
    }
    #else
    template <class Model>
    IOResult<void> export_input_data_county_timeseries(
        std::vector<Model> models, 
        const std::string& results_dir,
        /* ... other parameters ... */) {
        return failure(StatusCode::UnknownError, 
                      "HDF5 not available");
    }
    #endif // MEMILIO_HAS_HDF5

**Parallel Computing**

.. code-block:: cpp

    #ifdef MEMILIO_ENABLE_OPENMP
    #include <omp.h>
    
    void parallel_simulation() {
        #pragma omp parallel for
        for (int i = 0; i < num_regions; ++i) {
            // Parallel execution of regional simulations
        }
    }
    #endif // MEMILIO_ENABLE_OPENMP

**Check Multiple Features Simultaneously**

.. code-block:: cpp

    // Check if both JSON and HDF5 are available
    #if defined(MEMILIO_HAS_JSONCPP) && defined(MEMILIO_HAS_HDF5)
    IOResult<void> read_and_export_data() {
        // Implementation using both JSON input and HDF5 output
    }
    #elif defined(MEMILIO_HAS_JSONCPP)
    IOResult<void> read_and_export_data() {
        // JSON-only implementation
    }
    #else
    IOResult<void> read_and_export_data() {
        return failure(StatusCode::UnknownError, 
                      "Neither JSON nor HDF5 available");
    }
    #endif
