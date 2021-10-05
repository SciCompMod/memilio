cmake_minimum_required(VERSION 3.10)

project(memilio VERSION 0.1.0)

option(EPI_BUILD_TESTS "Build memilio unit tests." ON)
option(EPI_BUILD_EXAMPLES "Build memilio examples." ON)
option(EPI_BUILD_MODELS "Build memilio models." ON)
option(EPI_BUILD_SIMULATIONS "Build memilio unit tests." ON)
option(EPI_USE_BUNDLED_SPDLOG "Use spdlog bundled with epi" ON)
option(EPI_USE_BUNDLED_EIGEN "Use eigen bundled with epi" ON)
option(EPI_USE_BUNDLED_BOOST "Use boost bundled with epi (only for epi-io)" ON)
option(EPI_USE_BUNDLED_JSONCPP "Use jsoncpp bundled with epi (only for epi-io)" ON)
option(EPI_SANITIZE_ADDRESS "Enable address sanitizer." OFF)
option(EPI_SANITIZE_UNDEFINED "Enable undefined behavior sanitizer." OFF)

mark_as_advanced(EPI_USE_BUNDLED_SPDLOG EPI_SANITIZE_ADDRESS EPI_SANITIZE_UNDEFINED)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# code coverage analysis
# Note: this only works under linux and with make
# Ninja creates different directory names which do not work together with this scrupt
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    option (EPI_TEST_COVERAGE "Enable GCov coverage analysis (adds a 'coverage' target)" OFF)
    mark_as_advanced(EPI_TEST_COVERAGE)
    if (EPI_TEST_COVERAGE)
        message(STATUS "Coverage enabled")
        include(CodeCoverage)
        append_coverage_compiler_flags()
        setup_target_for_coverage_lcov(
            NAME coverage
            EXECUTABLE memilio-test
            EXCLUDE "${CMAKE_SOURCE_DIR}/tests*" "${CMAKE_SOURCE_DIR}/simulations*" "${CMAKE_SOURCE_DIR}/examples*" "${CMAKE_BINARY_DIR}/*" "/usr*"
        )
    endif()
endif()

if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7))
    if(EPI_SANITIZE_ADDRESS)
        string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address")
        string(APPEND CMAKE_LINKER_FLAGS_DEBUG  " -fsanitize=address")
    endif(EPI_SANITIZE_ADDRESS)

    if(EPI_SANITIZE_UNDEFINED)
        string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=undefined")
        string(APPEND CMAKE_LINKER_FLAGS_DEBUG  " -fsanitize=undefined")
    endif(EPI_SANITIZE_UNDEFINED)
    
    if(EPI_SANITIZE_ADDRESS OR EPI_SANITIZE_UNDEFINED)
        string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fno-omit-frame-pointer -fno-sanitize-recover=all")
        string(APPEND CMAKE_LINKER_FLAGS_DEBUG  " -fno-omit-frame-pointer -fno-sanitize-recover=all")
    endif(EPI_SANITIZE_ADDRESS OR EPI_SANITIZE_UNDEFINED)
endif((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7))

include(thirdparty/CMakeLists.txt)

set(_EPI_CXX_FLAG_ENABLE_WARNING_ERRORS
    "$<$<CXX_COMPILER_ID:GNU>: -Wno-unknown-warning -Wall -Wextra -Werror -Wshadow --pedantic-errors -Wno-deprecated-copy>
    $<$<CXX_COMPILER_ID:CLANG>: -Wno-unknown-warning-option -Wall -Wextra -Werror -Wshadow --pedantic-errors -Wno-deprecated>
    $<$<CXX_COMPILER_ID:MSVC>: /W4 /WX>")

add_subdirectory(memilio)
if (EPI_BUILD_MODELS)
    add_subdirectory(models/abm)
    add_subdirectory(models/secir)
    add_subdirectory(models/seir)
endif()
if (EPI_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
if (EPI_BUILD_TESTS)
    add_subdirectory(tests)
endif()
if (EPI_BUILD_SIMULATIONS)
    add_subdirectory(simulations)
endif()

include(GNUInstallDirs)

install(TARGETS memilio
        EXPORT memilio-targets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY memilio DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN memilio/*/*.h)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/memilio DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN memilio/*/*.h)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_CURRENT_LIST_DIR}/cmake/memilio-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/memilio-config.cmake
INSTALL_DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/memilio
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/memilio-config-version.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

install (
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/memilio-config-version.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/memilio-config.cmake"
  DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/memilio
)
