set (BOOST_DIR ${PROJECT_BINARY_DIR}/boost_1_75_0)
set (BOOST_ARCHIVE ${PROJECT_SOURCE_DIR}/thirdparty/boost_1_75_0.tar.gz)

if (NOT EXISTS ${BOOST_DIR})
    message(STATUS "Extracting boost")
    execute_process (
        COMMAND ${CMAKE_COMMAND} -E  tar xzf ${BOOST_ARCHIVE}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    )
endif()

# boost
add_library(boost INTERFACE)
add_library(Boost::boost ALIAS boost)
target_include_directories(boost INTERFACE
    $<BUILD_INTERFACE:${BOOST_DIR}>
)

add_library(boost_disable_autolink INTERFACE)
target_compile_definitions(boost_disable_autolink INTERFACE BOOST_ALL_NO_LIB)
add_library(Boost::disable_autolinking ALIAS boost_disable_autolink)

add_library(boost_filesystem 
    ${BOOST_DIR}/libs/filesystem/src/codecvt_error_category.cpp
    ${BOOST_DIR}/libs/filesystem/src/directory.cpp
    ${BOOST_DIR}/libs/filesystem/src/exception.cpp
    ${BOOST_DIR}/libs/filesystem/src/operations.cpp
    ${BOOST_DIR}/libs/filesystem/src/path.cpp
    ${BOOST_DIR}/libs/filesystem/src/path_traits.cpp
    ${BOOST_DIR}/libs/filesystem/src/portability.cpp
    ${BOOST_DIR}/libs/filesystem/src/unique_path.cpp
    ${BOOST_DIR}/libs/filesystem/src/utf8_codecvt_facet.cpp
    ${BOOST_DIR}/libs/filesystem/src/windows_file_codecvt.cpp
)

target_link_libraries(boost_filesystem PUBLIC boost_disable_autolink boost)
set_property(TARGET boost_filesystem PROPERTY POSITION_INDEPENDENT_CODE ON)
add_library(Boost::filesystem ALIAS boost_filesystem)

set(Boost_FOUND ON)
