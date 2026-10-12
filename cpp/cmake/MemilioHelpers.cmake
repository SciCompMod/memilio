# Helper functions used by the CMake files of memilio.

# Applies memilio's warning flags to a target, see MEMILIO_ENABLE_WARNINGS and MEMILIO_ENABLE_WARNINGS_AS_ERRORS.
# Warnings as errors can also be disabled for a single build using `cmake --compile-no-warning-as-error`.
function(memilio_target_warnings target)
    target_compile_options(${target} PRIVATE ${MEMILIO_CXX_WARNING_FLAGS})
    set_target_properties(${target} PROPERTIES COMPILE_WARNING_AS_ERROR ${MEMILIO_ENABLE_WARNINGS_AS_ERRORS})
endfunction()

# Adds all headers found in the current source directory to the public HEADERS file set of a target.
# The file set provides the include directories (BASE_DIRS) for building and installing, and installs the headers.
# Headers are collected by globbing, so that every header is installed, even if it is not listed in the sources.
#
# memilio_target_headers(<target> BASE_DIRS <dirs>... [EXCLUDE <regex>] [FILES <additional headers>...])
function(memilio_target_headers target)
    cmake_parse_arguments(PARSE_ARGV 1 ARG "" "EXCLUDE" "BASE_DIRS;FILES")
    file(GLOB_RECURSE headers CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp" "${CMAKE_CURRENT_SOURCE_DIR}/*.ipp")
    if(ARG_EXCLUDE)
        list(FILTER headers EXCLUDE REGEX "${ARG_EXCLUDE}")
    endif()
    get_target_property(type ${target} TYPE)
    if(type STREQUAL "INTERFACE_LIBRARY")
        set(scope INTERFACE)
    else()
        set(scope PUBLIC)
    endif()
    target_sources(${target} ${scope} FILE_SET HEADERS BASE_DIRS ${ARG_BASE_DIRS} FILES ${headers} ${ARG_FILES})
endfunction()

# Adds a model library, e.g. `memilio_add_model(ode_seir model.h model.cpp ...)`.
# Models are available as `<name>` and `memilio::<name>`, link against memilio and are installed with memilio.
function(memilio_add_model name)
    add_library(${name} ${ARGN})
    add_library(memilio::${name} ALIAS ${name})
    target_link_libraries(${name} PUBLIC memilio)
    # models are included as e.g. "ode_seir/model.h", so the base directory is the models directory
    memilio_target_headers(${name} BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/..")
    memilio_target_warnings(${name})
    set_property(GLOBAL APPEND PROPERTY MEMILIO_MODEL_TARGETS ${name})
endfunction()
