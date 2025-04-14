#-------------------------------------------------------
# Prerequisite Checks & Setup
#-------------------------------------------------------
include(GoogleTest)
set(MEMILIO_GTEST_VERSION "v1.12.1" CACHE STRING "GoogleTest version to use")

#-------------------------------------------------------
# GoogleTest Dependency Setup (using FetchContent)
#-------------------------------------------------------
message(STATUS "Configuring GoogleTest dependency for tests...")

# Force shared CRT for GoogleTest on Windows - set BEFORE MakeAvailable.
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        ${MEMILIO_GTEST_VERSION}
    GIT_SHALLOW    TRUE
)

# Make googletest targets available (downloads, configures, adds subdirectory)
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS TRUE)
FetchContent_MakeAvailable(googletest)
unset(CMAKE_SUPPRESS_DEVELOPER_WARNINGS)

#-------------------------------------------------------
# Test Source File Organization
#-------------------------------------------------------
set(CORE_TESTS
    testmain.cpp
    test_eigen_util.cpp
    test_stl_util.cpp
    test_uncertain.cpp
    test_math_floating_point.cpp
    test_math_time_series_functor.cpp
    test_type_safe.cpp
    test_custom_index_array.cpp
    test_matrix_shape.cpp
    test_metaprogramming.cpp
    test_utils.cpp
    test_date.cpp
    test_time_series.cpp
    test_transform_iterator.cpp
    test_history.cpp
    test_random_number_generator.cpp
)

set(EBM_TESTS
    test_populations.cpp
    test_odeseir.cpp
    test_odesir.cpp
    test_numericalIntegration.cpp
    test_smoother.cpp
    test_damping.cpp
    test_odesecir.cpp
    test_odesecirvvs.cpp
    test_odesecirts.cpp
    test_sde_sir.cpp
    test_sde_sirs.cpp
    test_sde_seirvv.cpp
    test_mobility.cpp
    test_odesecir_ageres.cpp
    test_odeseair.cpp
    test_parameter_studies.cpp
    test_dynamic_npis.cpp
    test_flows.cpp
    test_feedback.cpp
    test_parameter_set.cpp
    test_damping_sampling.cpp
    test_compartmentsimulation.cpp
    test_ide_seir.cpp
    test_ide_secir.cpp
    test_ide_secir_ageres.cpp
    test_state_age_function.cpp
    test_lct_secir.cpp
    test_lct_initializer_flows.cpp
    test_glct_secir.cpp
    test_smm_model.cpp
)

# set(NETWORK_TESTS
#     test_graph.cpp
#     test_graph_simulation.cpp
#     test_graph_abm.cpp
# )

# set(ABM_TESTS
#     test_abm_household.cpp
#     test_abm_infection.cpp
#     test_abm_location.cpp
#     test_abm_lockdown_rules.cpp
#     test_abm_masks.cpp
#     test_abm_mobility_rules.cpp
#     test_abm_person.cpp
#     test_abm_simulation.cpp
#     test_abm_serialization.cpp
#     test_abm_testing_strategy.cpp
#     test_abm_model.cpp
#     test_d_abm_model.cpp
#     abm_helpers.h
#     abm_helpers.cpp
# )

# set(IO_TESTS
#     test_regions.cpp
#     test_io_cli.cpp
#     test_io_framework.cpp
#     test_binary_serializer.cpp
#     test_mobility_io.cpp
#     test_contact_matrix.cpp
#     test_analyze_result.cpp
# )

set(TEST_UTILS
    actions.h
    distributions_helpers.h
    distributions_helpers.cpp
    matchers.cpp
    matchers.h
    random_number_test.h
    sanitizers.cpp
    temp_file_register.h
)

# Conditionally add tests dependent on optional libraries
set(JSONCPP_TESTS "")
if(MEMILIO_HAS_JSONCPP)
    list(APPEND JSONCPP_TESTS
        test_json_serializer.cpp
        test_epi_data_io.cpp
        test_lct_parameters_io.cpp
        test_ide_parameters_io.cpp
    )
endif()

set(JSONCPP_HDF5_TESTS "")
if(MEMILIO_HAS_JSONCPP AND MEMILIO_HAS_HDF5)
    list(APPEND JSONCPP_HDF5_TESTS
        test_save_parameters.cpp
        test_save_results.cpp
    )
endif()

# Combine all sources
set(ALL_TEST_SOURCES
    ${CORE_TESTS}
    ${EBM_TESTS}
    ${NETWORK_TESTS}
    ${ABM_TESTS}
    ${IO_TESTS}
    ${TEST_UTILS}
    ${JSONCPP_TESTS}
    ${JSONCPP_HDF5_TESTS}
)

#-------------------------------------------------------
# Test Executable Definition
#-------------------------------------------------------
add_executable(memilio-test ${ALL_TEST_SOURCES})

#-------------------------------------------------------
# Target Configuration (Includes, Links, Flags)
#-------------------------------------------------------
target_include_directories(memilio-test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(memilio-test PRIVATE
    # Memilio Models (they lready link to the required libraries publicly)
    memilio
    ode_secir ode_seir ode_secirvvs ode_secirts ode_seair
    ide_seir ide_secir lct_secir glct_secir
    abm graph_abm d_abm smm
    # Testing Framework (gmock already includes gtest_main)
    GTest::gtest GTest::gmock 
)

# Apply project-wide flags via helper functions (defined in root)
memilio_set_warning_flags(memilio-test)
memilio_set_coverage_flags(memilio-test)


#-------------------------------------------------------
# Test Discovery and CTest Integration
#-------------------------------------------------------
gtest_discover_tests(memilio-test
    TEST_PREFIX "memilio."
    PROPERTIES
        FOLDER "Tests"
        ENVIRONMENT "CTEST_OUTPUT_ON_FAILURE=1"
    DISCOVERY_TIMEOUT 120
)

#-------------------------------------------------------
# Custom Test Runner Target (`memilio-check`)
#-------------------------------------------------------
if(CMAKE_CONFIGURATION_TYPES) # Multi-config generators
    add_custom_target(memilio-check ALL
        COMMAND ${CMAKE_CTEST_COMMAND}
                --force-new-ctest-process --output-on-failure
                --build-config "$<CONFIGURATION>"
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running memilio tests ($<CONFIG>) with output on failure..."
        VERBATIM
    )
else() # Single-config generators
    add_custom_target(memilio-check ALL
        COMMAND ${CMAKE_CTEST_COMMAND}
                --force-new-ctest-process --output-on-failure
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running memilio tests with output on failure..."
        VERBATIM
    )
endif()

#-------------------------------------------------------
# IDE Cleanup (Hide Internal/External Targets)
#-------------------------------------------------------
mark_as_advanced(
    gtest_force_shared_crt
    gtest_hide_internal_symbols
    INSTALL_GTEST
    BUILD_GMOCK
    BUILD_GTEST
)