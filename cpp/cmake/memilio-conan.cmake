include(conan)

set(MEMILIO_CONAN_BUILD "missing" CACHE STRING "Semicolon-separated list of packages that are built from source by conan instead of loaded as binary packages. Equivalent to 'conan install --build <package>', see conan documentation. Default 'missing', i.e. packages are only built if binary package is not available.")
mark_as_advanced(MEMILIO_CONAN_BUILD)

#memilio_conan_install(REQUIRES x y z OPTIONS a b c) installs the conan packages x, y, and z with options a, b, and c.
#packages are specified as the usual conan package reference, i.e. <name>/<version>.
#options are specified the usual conan way, i.e. <packagename>:<option>=<value>.
macro(memilio_conan_install)
    cmake_parse_arguments(_MCI "" "" "REQUIRES;OPTIONS" ${ARGN})

    conan_cmake_configure(
        REQUIRES ${_MCI_REQUIRES}
        GENERATORS cmake_find_package_multi
        OPTIONS ${_MCI_OPTIONS})

    # install conan packages
    get_property(_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
    if (_MULTI_CONFIG)
        message(STATUS "Memilio: Downloading packages for all configurations.")
        foreach(TYPE ${CMAKE_CONFIGURATION_TYPES})
            conan_cmake_autodetect(_MEMILIO_SETTINGS BUILD_TYPE ${TYPE})
            conan_cmake_install(
                PATH_OR_REFERENCE .
                BUILD ${MEMILIO_CONAN_BUILD}
                REMOTE conancenter
                SETTINGS ${_MEMILIO_SETTINGS})
        endforeach()
    else()
        message(STATUS "Memilio: Downloading packages.")
        if (NOT CMAKE_BUILD_TYPE)
            message(WARNING "Memilio: Conan requires CMAKE_BUILD_TYPE to be set. Setting default to Release. Overwrite with e.g. -DCMAKE_BUILD_TYPE=Debug")
            set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
        endif()
        conan_cmake_autodetect(_MEMILIO_SETTINGS)
        conan_cmake_install(
            PATH_OR_REFERENCE .
            BUILD ${MEMILIO_CONAN_BUILD}
            REMOTE conancenter
            SETTINGS ${_MEMILIO_SETTINGS})
    endif()

    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_BINARY_DIR})
    list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_BINARY_DIR})
endmacro()
