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
    geography/regions.h
    geography/regions.cpp
    epidemiology/simulation_day.h
    geography/holiday_data.ipp
    compartments/compartmentalmodel.h
    compartments/flow_model.h
    compartments/simulation.h
    compartments/flow_simulation.h
    compartments/parameter_studies.h
    io/default_serialize.h
    io/default_serialize.cpp
    io/io.h
    io/io.cpp
    io/hdf5_cpp.h
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
    mobility/metapopulation_mobility_instant.h
    mobility/metapopulation_mobility_instant.cpp
    mobility/metapopulation_mobility_stochastic.h
    mobility/metapopulation_mobility_stochastic.cpp
    mobility/graph_simulation.h
    mobility/graph_simulation.cpp
    mobility/graph.h
    mobility/graph.cpp
    utils/visitor.h
    utils/uncertain_value.h
    utils/uncertain_value.cpp
    utils/index_range.h
    utils/index.h
    utils/flow.h
    utils/custom_index_array.h
    utils/memory.h
    utils/parameter_distributions.h
    utils/parameter_distribution_wrapper.h
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
)

target_include_directories(memilio PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # includes configured by cmake
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_features(memilio PUBLIC cxx_std_17)
target_link_libraries(memilio PUBLIC spdlog::spdlog Eigen3::Eigen Boost::boost Boost::filesystem Boost::disable_autolinking Random123 AD::AD)
target_compile_options(memilio
    PRIVATE
    ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS}
    PUBLIC
    $<$<CXX_COMPILER_ID:MSVC>: /bigobj>
)

if(MEMILIO_HAS_HDF5)
    target_link_libraries(memilio PUBLIC ${HDF5_C_LIBRARIES})
    target_include_directories(memilio PUBLIC ${HDF5_INCLUDE_DIRS})
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
