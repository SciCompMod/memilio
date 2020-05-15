cmake_minimum_required(VERSION 3.4)

project(epidemiology VERSION 0.1.0)

option(EPI_USE_BUNDLED_SPDLOG "Use spdlog bundled with epi" ON)
mark_as_advanced(EPI_USE_BUNDLED_SPDLOG)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(thirdparty)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)


add_library(epidemiology
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
)
target_include_directories(epidemiology PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_features(epidemiology PRIVATE cxx_std_14)
target_link_libraries(epidemiology PUBLIC spdlog::spdlog
                                          Eigen3::Eigen)

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
