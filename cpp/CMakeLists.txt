cmake_minimum_required(VERSION 3.4)

project(epidemiology VERSION 0.1.0)

option(EPI_USE_BUNDLED_SPDLOG "Use spdlog bundled with epi" ON)
option(EPI_USE_BUNDLED_EIGEN "Use eigen bundled with epi" ON)
option(EPI_BUILD_EPI_IO "Build epidemiology_io library (requires HDF5)" OFF)
option(EPI_USE_BUNDLED_BOOST "Use boost bundled with epi (only for epi-io)" ON)
option(EPI_USE_BUNDLED_JSONCPP "Use jsoncpp bundled with epi (only for epi-io)" ON)

mark_as_advanced(EPI_USE_BUNDLED_SPDLOG)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# code coverage analysis
# Note: this only works under linux and with make
# Ninja creates different directory names which do not work together with this scrupt
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    option (EPI_TEST_COVERAGE "Enable GCov coverage analysis (adds a 'coverage' target)" OFF)
    if (EPI_TEST_COVERAGE)
        message(STATUS "Coverage enabled")
        include(CodeCoverage)
        append_coverage_compiler_flags()
        setup_target_for_coverage_lcov(
            NAME coverage
            EXECUTABLE runUnitTests
            EXCLUDE "${CMAKE_SOURCE_DIR}/tests*" "${CMAKE_SOURCE_DIR}/examples*" "${CMAKE_BINARY_DIR}/*" "/usr*"
        )
    endif()
endif()

include(thirdparty/CMakeLists.txt)

add_library(epidemiology
    epidemiology/math/euler.cpp
    epidemiology/math/euler.h
    epidemiology/math/smoother.h    
    epidemiology/utils/ScalarType.h
    epidemiology/model/compartmentalmodel.h
    epidemiology/model/populations.h
    epidemiology/model/simulation.h
    epidemiology/secir/age_group.h
    epidemiology/secir/infection_state.h
    epidemiology/secir/damping.cpp
    epidemiology/secir/damping.h
    epidemiology/secir/secir_params.h
    epidemiology/secir/secir.h
    epidemiology/secir/seir.h
    epidemiology/secir/seir.cpp
    epidemiology/math/adapt_rk.cpp
    epidemiology/math/adapt_rk.h
    epidemiology/math/integrator.h
    epidemiology/math/integrator.cpp
    epidemiology/migration/migration.h
    epidemiology/migration/migration.cpp
    epidemiology/migration/graph_simulation.h
    epidemiology/migration/graph_simulation.cpp
    epidemiology/utils/graph.h
    epidemiology/utils/graph.cpp
    epidemiology/utils/eigen.h
    epidemiology/utils/eigen_util.h
    epidemiology/utils/visitor.h
    epidemiology/secir/parameter_space.h
    epidemiology/secir/parameter_space.cpp
    epidemiology/utils/uncertain_value.h
    epidemiology/utils/uncertain_value.cpp
    epidemiology/secir/uncertain_matrix.h    
    epidemiology/secir/uncertain_matrix.cpp
    epidemiology/secir/contact_matrix.h    
    epidemiology/secir/contact_matrix.cpp
    epidemiology/secir/damping_sampling.h    
    epidemiology/secir/damping_sampling.cpp
    epidemiology/secir/dynamic_npis.h    
    epidemiology/secir/dynamic_npis.cpp
    epidemiology/secir/implicit_euler.h    
    epidemiology/secir/implicit_euler.cpp
    epidemiology/utils/index.h
    epidemiology/utils/custom_index_array.h
    epidemiology/utils/memory.h
    epidemiology/utils/parameter_distributions.h
    epidemiology/utils/time_series.h
    epidemiology/utils/time_series.cpp
    epidemiology/utils/span.h
    epidemiology/utils/span.cpp
    epidemiology/utils/type_safe.h
    epidemiology/utils/type_safe.cpp
    epidemiology/utils/parameter_set.h
    epidemiology/utils/parameter_set.cpp
    epidemiology/utils/date.h
    epidemiology/utils/date.cpp
    epidemiology/utils/regions.h
    epidemiology/utils/regions.cpp
    epidemiology/utils/holiday_data_de.ipp
    epidemiology/utils/matrix_shape.h
    epidemiology/utils/matrix_shape.cpp
    epidemiology/utils/random_number_generator.h
    epidemiology/utils/random_number_generator.cpp
    epidemiology/abm/location.cpp
    epidemiology/abm/location.h
    epidemiology/abm/simulation.cpp
    epidemiology/abm/simulation.h
    epidemiology/abm/person.cpp
    epidemiology/abm/person.h
    epidemiology/abm/world.cpp
    epidemiology/abm/world.h
    epidemiology/abm/state.h
    epidemiology/abm/location_type.h
    epidemiology/abm/parameters.h
    epidemiology/abm/parameters.cpp
    epidemiology/abm/migration_rules.cpp
    epidemiology/abm/migration_rules.h
    epidemiology/abm/lockdown_rules.cpp
    epidemiology/abm/lockdown_rules.h
    epidemiology/abm/testing_scheme.cpp
    epidemiology/abm/testing_scheme.h
    epidemiology/abm/world_builder.cpp
    epidemiology/abm/world_builder.h
    epidemiology/secir/parameter_studies.h
    epidemiology/secir/analyze_result.h
    epidemiology/secir/analyze_result.cpp
)

target_include_directories(epidemiology PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)


target_compile_features(epidemiology PUBLIC cxx_std_14)
target_link_libraries(epidemiology PUBLIC spdlog::spdlog Eigen3::Eigen Boost::boost)
target_compile_options(epidemiology PRIVATE
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>: -Wall -Wextra -Werror -Wshadow --pedantic-errors -Wno-deprecated-copy>
    $<$<CXX_COMPILER_ID:MSVC>: /W4 /WX>)

if (EPI_BUILD_EPI_IO)

    add_library(epidemiology_io
      epidemiology_io/secir_result_io.cpp
      epidemiology_io/secir_result_io.h
      epidemiology_io/secir_parameters_io.cpp
      epidemiology_io/secir_parameters_io.h
      epidemiology_io/io.cpp
      epidemiology_io/io.h
      epidemiology_io/mobility_io.cpp
      epidemiology_io/mobility_io.h
      epidemiology_io/json_serializer.h
      epidemiology_io/json_serializer.cpp
    )



    target_include_directories(epidemiology_io PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )

    target_link_libraries(epidemiology_io PUBLIC epidemiology jsoncpp Boost::disable_autolinking
                          PRIVATE ${HDF5_CXX_LIBRARIES} Boost::filesystem)
    target_include_directories(epidemiology_io PRIVATE ${HDF5_INCLUDE_DIRS})

    target_compile_definitions(epidemiology_io PUBLIC "-DHAVE_EPI_IO")
    target_compile_options(epidemiology_io PRIVATE
        $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>: -Wall -Wextra -Werror -Wshadow --pedantic-errors -Wno-deprecated-copy>
        $<$<CXX_COMPILER_ID:MSVC>: /W4 /WX>)

endif(EPI_BUILD_EPI_IO)

include(GNUInstallDirs)

add_subdirectory(examples)
add_subdirectory(tests)
# add_subdirectory(docs)

install(TARGETS epidemiology
        EXPORT epidemiology-targets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY epidemiology DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

include(CMakePackageConfigHelpers)


configure_package_config_file(
    ${CMAKE_CURRENT_LIST_DIR}/cmake/epidemiology-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/epidemiology-config.cmake
INSTALL_DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/epidemiology
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/epidemiology-config-version.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

#install(
#    EXPORT epidemiology-targets
#    FILE epidemiology-targets.cmake
#    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/epidemiology
#)

install (
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/epidemiology-config-version.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/epidemiology-config.cmake"
  DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/epidemiology
)
