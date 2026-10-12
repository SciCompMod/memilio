configure_file(config_internal.h.in memilio/config_internal.h)

add_library(memilio
    config.h
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
    geography/distance.h
    geography/regions.h
    geography/regions.cpp
    geography/geolocation.h
    geography/rtree.h
    geography/holiday_data.ipp
    compartments/compartmental_model.h
    compartments/flow_model.h
    compartments/flow_simulation_base.h
    compartments/flow_simulation.h
    compartments/simulation_base.h
    compartments/simulation.h
    compartments/stochastic_simulation.h
    compartments/stochastic_model.h
    compartments/parameter_studies.h
    io/binary_serializer.h
    io/binary_serializer.cpp
    io/cli.cpp
    io/cli.h
    io/default_serialize.h
    io/default_serialize.cpp
    io/directories.cpp
    io/directories.h
    io/epi_data.h
    io/epi_data.cpp
    io/hdf5_cpp.h
    io/io.h
    io/io.cpp
    io/json_serializer.h
    io/json_serializer.cpp
    io/history.h
    io/mobility_io.h
    io/mobility_io.cpp
    io/parameters_io.h
    io/parameters_io.cpp
    io/result_io.h
    io/result_io.cpp
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
    math/math_utils.h
    mobility/metapopulation_mobility_instant.h
    mobility/metapopulation_mobility_instant.cpp
    mobility/metapopulation_mobility_stochastic.h
    mobility/metapopulation_mobility_stochastic.cpp
    mobility/graph_simulation.h
    mobility/graph_simulation.cpp
    mobility/graph.h
    mobility/graph.cpp
    mobility/graph_builder.h
    mobility/graph_builder.cpp
    timer/auto_timer.h
    timer/basic_timer.cpp
    timer/basic_timer.h
    timer/definitions.cpp
    timer/definitions.h
    timer/list_printer.h
    timer/named_timer.h
    timer/registration.h
    timer/timer_registrar.h
    utils/visitor.h
    utils/uncertain_value.h
    utils/uncertain_value.cpp
    utils/index_range.h
    utils/index.h
    utils/flow.h
    utils/custom_index_array.h
    utils/memory.h
    utils/parameter_distributions.h
    utils/abstract_parameter_distribution.h
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
    utils/mioomp.cpp
    utils/string_literal.h
    utils/type_list.h
    utils/back_inserter_second_element.h
)

add_library(memilio::memilio ALIAS memilio)

# headers are included as "memilio/...", configured headers are found in the binary dir; AD has its own target
memilio_target_headers(memilio
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/.." "${CMAKE_CURRENT_BINARY_DIR}"
    EXCLUDE "/memilio/ad/include/"
    FILES "${CMAKE_CURRENT_BINARY_DIR}/memilio/config_internal.h"
)

target_link_libraries(memilio PUBLIC spdlog::spdlog Eigen3::Eigen Boost::boost Boost::disable_autolinking Random123 AD::AD)
# parameter_distributions.h uses the ABM's random number generator, which is part of the abm model ("abm/...")
# TODO: remove this dependency of the core library on a model
target_include_directories(memilio PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/models>)

# C++20 is required to use memilio, also for projects that link against it
target_compile_features(memilio PUBLIC cxx_std_20)
memilio_target_warnings(memilio)
target_compile_options(memilio PUBLIC $<$<CXX_COMPILER_ID:MSVC>:/bigobj>)

if(MEMILIO_HAS_HDF5)
    target_link_libraries(memilio PUBLIC HDF5::HDF5)
endif()

if(MEMILIO_HAS_JSONCPP)
    target_link_libraries(memilio PUBLIC JsonCpp::JsonCpp)
endif()

if(MEMILIO_ENABLE_MPI)
    target_link_libraries(memilio PUBLIC MPI::MPI_CXX)
endif()

if(MEMILIO_ENABLE_OPENMP)
    target_link_libraries(memilio PUBLIC OpenMP::OpenMP_CXX)
endif()

if(MEMILIO_ENABLE_LIKWID_MARKER)
    target_link_libraries(memilio PUBLIC likwid::likwid)
    target_compile_definitions(memilio PUBLIC LIKWID_PERFMON)
endif()
