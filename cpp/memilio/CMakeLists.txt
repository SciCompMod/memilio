# Configure internal header based on CMake settings
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config_internal.h.in" 
               "${CMAKE_CURRENT_BINARY_DIR}/memilio/config_internal.h"
               )

add_library(memilio
    # Config & Core
    config.h
    memilio/config_internal.h # Generated file

    # Data & Epidemiology
    data/analyze_result.h
    data/analyze_result.cpp
    epidemiology/age_group.h
    epidemiology/populations.h
    epidemiology/damping.cpp
    epidemiology/damping.h
    epidemiology/uncertain_matrix.h
    epidemiology/uncertain_matrix.cpp
    epidemiology/contact_matrix.h
    epidemiology/contact_matrix.cpp
    epidemiology/damping_sampling.h
    epidemiology/damping_sampling.cpp
    epidemiology/dynamic_npis.h
    epidemiology/dynamic_npis.cpp
    epidemiology/lct_infection_state.h
    epidemiology/lct_populations.h
    epidemiology/adoption_rate.h
    epidemiology/simulation_day.h

    # Geography
    geography/regions.h
    geography/regions.cpp
    geography/holiday_data.ipp # Note: .ipp files are usually included by .h/.cpp

    # Compartments & Simulation
    compartments/compartmentalmodel.h
    compartments/flow_model.h
    compartments/simulation.h
    compartments/flow_simulation.h
    compartments/parameter_studies.h

    # IO
    io/default_serialize.h
    io/default_serialize.cpp
    io/io.h
    io/io.cpp
    io/hdf5_cpp.h # Helper, might not need explicit listing if header-only/internal
    io/json_serializer.h
    io/json_serializer.cpp
    io/binary_serializer.h
    io/binary_serializer.cpp
    io/history.h
    io/mobility_io.h
    io/mobility_io.cpp
    io/result_io.h
    io/result_io.cpp
    io/epi_data.h
    io/epi_data.cpp
    io/cli.h

    # Math
    math/euler.cpp
    math/euler.h
    math/smoother.h
    math/adapt_rk.h
    math/adapt_rk.cpp
    math/stepper_wrapper.h
    math/stepper_wrapper.cpp
    math/integrator.h
    math/integrator.cpp
    math/eigen.h
    math/eigen_util.h
    math/matrix_shape.h
    math/matrix_shape.cpp
    math/interpolation.h
    math/interpolation.cpp
    math/time_series_functor.h
    math/time_series_functor.cpp

    # Mobility & Graphs
    mobility/metapopulation_mobility_instant.h
    mobility/metapopulation_mobility_instant.cpp
    mobility/metapopulation_mobility_stochastic.h
    mobility/metapopulation_mobility_stochastic.cpp
    mobility/graph_simulation.h
    mobility/graph_simulation.cpp
    mobility/graph.h
    mobility/graph.cpp

    # Utilities
    utils/visitor.h
    utils/uncertain_value.h
    utils/uncertain_value.cpp
    utils/index_range.h
    utils/index.h
    utils/flow.h
    utils/custom_index_array.h
    utils/memory.h
    utils/parameter_distributions.h
    utils/time_series.h
    utils/time_series.cpp
    utils/span.h
    utils/span.cpp
    utils/type_safe.h
    utils/type_safe.cpp
    utils/parameter_set.h
    utils/parameter_set.cpp
    utils/date.h
    utils/date.cpp
    utils/random_number_generator.h
    utils/random_number_generator.cpp
    utils/miompi.h
    utils/miompi.cpp
    utils/mioomp.h
    utils/type_list.h

    #ad library
    ad/ad.hpp
)

# Include directories for build and install interfaces
target_include_directories(memilio PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..> # Project root include dir (e.g., for memilio/...)
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>    # For configured headers (config_internal.h)
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>  # Standard install include dir
)

# Link core dependencies
target_link_libraries(memilio PUBLIC
    spdlog::spdlog       # Logging
    Eigen3::Eigen        # Linear algebra
    Boost::boost         # General Boost utilities (if needed directly)
    Boost::filesystem    # Filesystem operations
    Boost::disable_autolinking # MSVC specific, ensure provided by FindBoost
    Random123            # Random number generation (ensure target exists)
)

# Apply project-wide compiler settings
memilio_set_warning_flags(memilio)
memilio_set_sanitizer_flags(memilio)


# Optional dependencies
if(MEMILIO_HAS_HDF5)
    message(STATUS "Linking memilio against HDF5")
    target_link_libraries(memilio PUBLIC HDF5::HDF5)
endif()

if(MEMILIO_HAS_JSONCPP)
    message(STATUS "Linking memilio against JsonCpp")
    target_link_libraries(memilio PUBLIC JsonCpp::JsonCpp)

if(MEMILIO_ENABLE_MPI)
    message(STATUS "Linking memilio against MPI")
    target_link_libraries(memilio PUBLIC MPI::MPI_CXX) 
endif()

if(MEMILIO_ENABLE_OPENMP)
    message(STATUS "Enabling OpenMP for memilio (interface)")
    target_link_libraries(memilio INTERFACE OpenMP::OpenMP_CXX)
endif()
