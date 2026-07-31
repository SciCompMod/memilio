#
#
# Downloads GTest and provides a helper macro to add tests. Add make check, as well, which
# gives output on failed tests without having to set an environment variable.
#
#
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

message(STATUS "Downloading GoogleTest library")
include(FetchContent)
FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.16.0)
FetchContent_MakeAvailable(googletest)

if(CMAKE_CONFIGURATION_TYPES)
    add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
        --force-new-ctest-process --output-on-failure
        --build-config "$<CONFIGURATION>")
else()
    add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
        --force-new-ctest-process --output-on-failure)
endif()

set_target_properties(check PROPERTIES FOLDER "Scripts")

# Target must already exist
macro(add_gtest TESTNAME)
    target_link_libraries(${TESTNAME} PUBLIC gtest gmock gtest_main)

    if(GOOGLE_TEST_INDIVIDUAL)    
        include(GoogleTest)
        gtest_discover_tests(${TESTNAME}
            TEST_PREFIX "${TESTNAME}."
            PROPERTIES FOLDER "Tests"
            DISCOVERY_TIMEOUT 30)
    else()
        add_test(${TESTNAME} ${TESTNAME})
        set_target_properties(${TESTNAME} PROPERTIES FOLDER "Tests")
    endif()
endmacro()

mark_as_advanced(
    gmock_build_tests
    gtest_build_samples
    gtest_build_tests
    gtest_disable_pthreads
    gtest_force_shared_crt
    gtest_hide_internal_symbols
    BUILD_GMOCK
    BUILD_GTEST
)

set_target_properties(gtest gtest_main gmock gmock_main
    PROPERTIES FOLDER "Extern")
