cmake_minimum_required(VERSION 3.4)

project(epidemiology VERSION 0.1.0)

option(EPI_USE_BUNDLED_SPDLOG "Use spdlog bundled with epi" ON)
option(EPI_USE_BUNDLED_EIGEN "Use eigen bundled with epi" ON)
option(EPI_BUILD_EPI_IO "Build epidemiology_io library (requires TiXI and HDF5)" OFF)

mark_as_advanced(EPI_USE_BUNDLED_SPDLOG)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(thirdparty/CMakeLists.txt)

add_library(epidemiology
    epidemiology/tensor_helpers.h
    epidemiology/tensor_helpers.cpp
    epidemiology/euler.cpp
    epidemiology/euler.h
    epidemiology/damping.cpp
    epidemiology/damping.h
    epidemiology/secir.cpp
    epidemiology/secir.h
    epidemiology/seir.cpp
    epidemiology/seir.h
    epidemiology/adapt_rk.cpp
    epidemiology/adapt_rk.h
    epidemiology/integrator.h
    epidemiology/integrator.cpp
    epidemiology/migration.h
    epidemiology/migration.cpp
    epidemiology/graph.h
    epidemiology/graph.cpp
    epidemiology/eigen_util.h
    epidemiology/populations.h
    epidemiology/populations.cpp
    epidemiology/visitor.h
    epidemiology/parameter_studies/parameter_space.h
    epidemiology/parameter_studies/parameter_space.cpp
    epidemiology/uncertain_value.cpp
    epidemiology/uncertain_value.h
    epidemiology/uncertain_matrix.cpp
    epidemiology/uncertain_matrix.h    
    epidemiology/memory.h
    epidemiology/parameter_studies/parameter_distributions.h
)

target_include_directories(epidemiology PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)


target_compile_features(epidemiology PRIVATE cxx_std_14)
target_link_libraries(epidemiology PUBLIC spdlog::spdlog Eigen3::Eigen)

if (EPI_BUILD_EPI_IO)
    
    add_library(epidemiology_io
      epidemiology_io/secir_result_io.cpp
      epidemiology_io/secir_result_io.h
      epidemiology_io/secir_parameters_io.cpp
      epidemiology_io/secir_parameters_io.h
      epidemiology_io/io.cpp
      epidemiology_io/io.h
      epidemiology_io/twitter_migration_io.cpp
      epidemiology_io/twitter_migration_io.h
    )

    target_include_directories(epidemiology_io PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )

    target_link_libraries(epidemiology_io PUBLIC epidemiology tixi3
                          PRIVATE ${HDF5_CXX_LIBRARIES} )
    target_include_directories(epidemiology_io PRIVATE ${HDF5_INCLUDE_DIRS})

    target_compile_definitions(epidemiology_io PUBLIC "-DHAVE_EPI_IO")

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
