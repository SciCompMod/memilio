enable_testing()
# Downloads GTest and provides a helper macro to add tests. Add make check, as well, which
# gives output on failed tests without having to set an environment variable.


# Download and unpack googletest at configure time
include(FetchContent)
FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.13.0 )
FetchContent_MakeAvailable(googletest)

if(MSVC)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif()

if(BUILD_SHARED_LIBS)
    target_compile_options(gtest PUBLIC -DGTEST_CREATE_SHARED_LIBRARY=1)
endif()

include(GoogleTest)
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

macro(add_gtest TESTNAME)
    target_link_libraries(${TESTNAME} PUBLIC gmock gtest)

    gtest_discover_tests(${TESTNAME}
        TEST_PREFIX "${TESTNAME}."
        PROPERTIES FOLDER "Tests"
        DISCOVERY_TIMEOUT 30)
endmacro()

set(TESTSOURCES
    testmain.cpp
    test_populations.cpp
    test_odeseir.cpp
    test_odesir.cpp
    test_numericalIntegration.cpp
    test_smoother.cpp
    test_damping.cpp
    test_odesecir.cpp
    test_odesecirvvs.cpp
    test_mobility.cpp
    test_date.cpp
    test_eigen_util.cpp
    test_odesecir_ageres.cpp
    test_parameter_studies.cpp
    test_graph.cpp
    test_graph_simulation.cpp
    test_stl_util.cpp
    test_uncertain.cpp
    test_random_number_generator.cpp
    test_time_series.cpp
    test_abm_household.cpp
    test_abm_infection.cpp
    test_abm_location.cpp
    test_abm_lockdown_rules.cpp
    test_abm_masks.cpp
    test_abm_migration_rules.cpp
    test_abm_person.cpp
    test_abm_simulation.cpp
    test_abm_testing_strategy.cpp
    test_abm_world.cpp
    test_analyze_result.cpp
    test_contact_matrix.cpp
    test_type_safe.cpp
    test_custom_index_array.cpp
    test_flows.cpp
    test_parameter_set.cpp
    test_matrix_shape.cpp
    test_damping_sampling.cpp
    test_dynamic_npis.cpp
    test_regions.cpp
    test_io_cli.cpp
    test_io_framework.cpp
    test_binary_serializer.cpp
    test_compartmentsimulation.cpp
    test_mobility_io.cpp
    test_transform_iterator.cpp
    test_metaprogramming.cpp
    test_ide_seir.cpp
    test_history.cpp
    test_utils.cpp
    abm_helpers.h
    abm_helpers.cpp
    test_ide_secir.cpp
    test_state_age_function.cpp
    distributions_helpers.h
    distributions_helpers.cpp
    actions.h
    matchers.h
    temp_file_register.h
    sanitizers.cpp
)

if(MEMILIO_HAS_JSONCPP)
set(TESTSOURCES ${TESTSOURCES}
test_json_serializer.cpp
test_epi_data_io.cpp
)
endif()

if(MEMILIO_HAS_JSONCPP AND MEMILIO_HAS_HDF5)
set(TESTSOURCES ${TESTSOURCES}
test_save_parameters.cpp
test_save_results.cpp
)
endif()

add_executable(memilio-test ${TESTSOURCES})
target_include_directories(memilio-test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(memilio-test PRIVATE memilio ode_secir ode_seir ode_secirvvs ide_seir ide_secir abm)
target_compile_options(memilio-test PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

# make unit tests find the test data files
file(TO_CMAKE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data" MEMILIO_TEST_DATA_DIR)
configure_file(test_data_dir.h.in test_data_dir.h)
target_include_directories(memilio-test PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

#compile for shared library
if(BUILD_SHARED_LIBS)   
    target_compile_definitions(memilio-test PRIVATE -DGTEST_LINKED_AS_SHARED_LIBRARY=1)
endif()

#add Executable to googletest
add_gtest(memilio-test)
