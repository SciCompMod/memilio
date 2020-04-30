cmake_minimum_required(VERSION 3.4)

project(epidemiology VERSION 0.1.0)

option(EPI_USE_BUNDLED_SPDLOG "Use spdlog bundled with epi" ON)
mark_as_advanced(EPI_USE_BUNDLED_SPDLOG)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

add_subdirectory(thirdparty)

add_library(epidemiology INTERFACE
)
target_include_directories(epidemiology INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_features(epidemiology INTERFACE cxx_std_11)
target_link_libraries(epidemiology INTERFACE spdlog::spdlog_header_only)

include(GNUInstallDirs)

add_subdirectory(examples)
add_subdirectory(tests)
# add_subdirectory(docs)

install(TARGETS epidemiology
        EXPORT epidemiology-targets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
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
