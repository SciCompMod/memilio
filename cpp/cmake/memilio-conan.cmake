include(conan)

macro(memilio_conan_install)
    cmake_parse_arguments(_MCI "" "" "REQUIRES" ${ARGN})

    conan_cmake_configure(
        REQUIRES ${_MCI_REQUIRES}
        GENERATORS cmake_find_package_multi)

    # install conan packages
    get_property(_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
    if (_MULTI_CONFIG)
        message(STATUS "Memilio: Downloading packages for all configurations.")
        foreach(TYPE ${CMAKE_CONFIGURATION_TYPES})
            conan_cmake_autodetect(_MEMILIO_SETTINGS BUILD_TYPE ${TYPE})
            conan_cmake_install(
                PATH_OR_REFERENCE .
                BUILD missing
                REMOTE conancenter
                SETTINGS ${_MEMILIO_SETTINGS})
        endforeach()
    else()
        message(STATUS "Memilio: Downloading packages.")
        if (NOT CMAKE_BUILD_TYPE)
            message(STATUS "MEMILIO WORKAROUND: Conan requires CMAKE_BUILD_TYPE to be set. Setting default to Release. Overwrite with e.g. -DCMAKE_BUILD_TYPE=<Debug>")
            set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
        endif()
        conan_cmake_autodetect(_MEMILIO_SETTINGS)
        conan_cmake_install(
            PATH_OR_REFERENCE .
            BUILD missing
            REMOTE conancenter
            SETTINGS ${_MEMILIO_SETTINGS})
    endif()

    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_BINARY_DIR})
    list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_BINARY_DIR})
endmacro()