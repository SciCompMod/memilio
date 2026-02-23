option(MEMILIO_USE_BUNDLED_PYBIND11 "Use pybind11 bundled with this library." ON)
mark_as_advanced(MEMILIO_USE_BUNDLED_PYBIND11)

include(GNUInstallDirs) # set to gnu folders. No cache variable so this is not global

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_PDB_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Windows has no RPATH (CMAKE_SKIP_RPATH defaults to ON there), so set this only on UNIX and APPLE.
if(APPLE)
    set(CMAKE_INSTALL_RPATH
        "@loader_path"
        "@loader_path/lib"
    )
elseif(UNIX)
    set(CMAKE_INSTALL_RPATH
        "$ORIGIN"
        "$ORIGIN/lib"
    )
endif()

# keep build RPATH separate and avoid putting build paths into installed wheels
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

# Use a static spdlog so the wheel stays self-contained.
set(SPDLOG_BUILD_SHARED OFF CACHE BOOL "" FORCE)

if(MEMILIO_USE_BUNDLED_PYBIND11)
    # Fetch pybind11
    include(FetchContent)
    FetchContent_Declare(pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11
        GIT_TAG v3.0.0
    )
    FetchContent_GetProperties(pybind11)

    if(NOT pybind11_POPULATED)
        FetchContent_Populate(pybind11)
        add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
    endif()
else()
    find_package(pybind11 REQUIRED)
endif()

# a list of all "LINKED_LIBRARIES" that are given to add_pymio_module. will contain duplicates
# used for wheel installation
set(PYMIO_MEMILIO_LIBS_LIST)

# Function to add a pybind11 module
function(add_pymio_module target_name)
    set(options)
    set(oneValueArgs)
    set(multiValueArgs LINKED_LIBRARIES SOURCES)
    cmake_parse_arguments(PYBIND11_MODULE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    pybind11_add_module(${target_name} MODULE ${PYBIND11_MODULE_SOURCES})
    target_link_libraries(${target_name} PRIVATE ${PYBIND11_MODULE_LINKED_LIBRARIES})
    target_include_directories(${target_name} PRIVATE memilio/simulation/bindings)

    # Collect linked libraries for later installation
    list(APPEND PYMIO_MEMILIO_LIBS_LIST "${PYBIND11_MODULE_LINKED_LIBRARIES}")
    set(PYMIO_MEMILIO_LIBS_LIST ${PYMIO_MEMILIO_LIBS_LIST} PARENT_SCOPE)

    # Set RPATH so extension can find the C++ libraries in lib subdirectory
    set_target_properties(${target_name} PROPERTIES
        INSTALL_RPATH "$ORIGIN/lib"
    )

    # Make sure to install shared and static libraries (LIBRARY and ARCHIVE) and the executables (RUNTIME).
    # Destination "." refers to the value of wheel.install-dir in the pyproject.toml.
    install(TARGETS ${target_name}
        LIBRARY DESTINATION .
        ARCHIVE DESTINATION .
        RUNTIME DESTINATION .
    )
endfunction()

# build python extensions
add_pymio_module(_simulation_abm
    LINKED_LIBRARIES memilio abm
    SOURCES memilio/simulation/bindings/models/abm.cpp
)

add_pymio_module(_simulation
    LINKED_LIBRARIES memilio
    SOURCES memilio/simulation/bindings/simulation.cpp
    memilio/simulation/bindings/epidemiology/damping_sampling.cpp
    memilio/simulation/bindings/epidemiology/uncertain_matrix.cpp
    memilio/simulation/bindings/mobility/metapopulation_mobility_instant.cpp
    memilio/simulation/bindings/utils/date.cpp
    memilio/simulation/bindings/utils/logging.cpp
    memilio/simulation/bindings/utils/time_series.cpp
    memilio/simulation/bindings/utils/parameter_distributions.cpp
    memilio/simulation/bindings/utils/uncertain_value.cpp
)

add_pymio_module(_simulation_osir
    LINKED_LIBRARIES memilio ode_sir
    SOURCES memilio/simulation/bindings/models/osir.cpp
)

add_pymio_module(_simulation_oseir
    LINKED_LIBRARIES memilio ode_seir
    SOURCES memilio/simulation/bindings/models/oseir.cpp
)

add_pymio_module(_simulation_oseirdb
    LINKED_LIBRARIES memilio ode_seirdb
    SOURCES memilio/simulation/bindings/models/oseirdb.cpp
)

add_pymio_module(_simulation_osecir
    LINKED_LIBRARIES memilio ode_secir
    SOURCES memilio/simulation/bindings/models/osecir.cpp
)

add_pymio_module(_simulation_osecirvvs
    LINKED_LIBRARIES memilio ode_secirvvs
    SOURCES memilio/simulation/bindings/models/osecirvvs.cpp
)

add_pymio_module(_simulation_ssir
    LINKED_LIBRARIES memilio sde_sir
    SOURCES memilio/simulation/bindings/models/ssir.cpp
)

add_pymio_module(_simulation_ssirs
    LINKED_LIBRARIES memilio sde_sirs
    SOURCES memilio/simulation/bindings/models/ssirs.cpp
)
add_pymio_module(_simulation_omseirs4
    LINKED_LIBRARIES memilio ode_mseirs4
    SOURCES memilio/simulation/bindings/models/omseirs4.cpp
)

# install all shared memilio libraries, which were given as "LINKED_LIBRARIES" to add_pymio_module
list(REMOVE_DUPLICATES PYMIO_MEMILIO_LIBS_LIST)

# Set RPATH for all C++ libraries so they can find each other
foreach(target_name IN LISTS PYMIO_MEMILIO_LIBS_LIST)
    if(TARGET ${target_name})
        set_target_properties(${target_name} PROPERTIES
            INSTALL_RPATH "$ORIGIN;$ORIGIN/lib"
        )
    endif()
endforeach()

install(TARGETS ${PYMIO_MEMILIO_LIBS_LIST}
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION lib
)

# If jsoncpp comes from the bundled third-party dependencies, copy its shared
# library into the wheel under every file name importers might look for. Wheels
# drop symlinks, so place real copies for each alias.
if(TARGET jsoncpp_lib)
    set_target_properties(jsoncpp_lib PROPERTIES
        INSTALL_RPATH "$ORIGIN"
    )
    install(FILES
        $<TARGET_FILE:jsoncpp_lib>
        $<TARGET_SONAME_FILE:jsoncpp_lib>
        $<TARGET_LINKER_FILE:jsoncpp_lib>
        DESTINATION lib
    )
endif()

if(TARGET spdlog)
    get_target_property(_SPDLOG_LIB_TYPE spdlog TYPE)

    if(_SPDLOG_LIB_TYPE STREQUAL "SHARED_LIBRARY")
        set_target_properties(spdlog PROPERTIES
            INSTALL_RPATH "$ORIGIN"
        )
        install(FILES
            $<TARGET_FILE:spdlog>
            $<TARGET_SONAME_FILE:spdlog>
            $<TARGET_LINKER_FILE:spdlog>
            DESTINATION lib
        )
    endif()
endif()
