enable_testing()
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(GOOGLE_TEST_INDIVIDUAL ON)
include(AddGoogleTest)

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
    test_odesecirts.cpp
    test_sde_sir.cpp
    test_sde_sirs.cpp
    test_sde_seirvv.cpp
    test_mobility.cpp
    test_date.cpp
    test_eigen_util.cpp
    test_odesecir_ageres.cpp
    test_odeseair.cpp
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
    test_abm_mobility_rules.cpp
    test_abm_person.cpp
    test_abm_simulation.cpp
    test_abm_serialization.cpp
    test_abm_testing_strategy.cpp
    test_abm_model.cpp
    test_math_floating_point.cpp
    test_math_time_series_functor.cpp
    test_analyze_result.cpp
    test_contact_matrix.cpp
    test_type_safe.cpp
    test_custom_index_array.cpp
    test_d_abm_model.cpp
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
    test_history.cpp
    test_utils.cpp
    test_ide_seir.cpp
    test_ide_secir.cpp
    test_ide_secir_ageres.cpp
    test_state_age_function.cpp
    test_lct_secir.cpp
    test_lct_initializer_flows.cpp
    test_glct_secir.cpp
    test_ad.cpp
    test_smm_model.cpp
    abm_helpers.h
    abm_helpers.cpp
    actions.h
    distributions_helpers.h
    distributions_helpers.cpp
    matchers.cpp
    matchers.h
    random_number_test.h
    sanitizers.cpp
    temp_file_register.h
)

if(MEMILIO_HAS_JSONCPP)
    set(TESTSOURCES ${TESTSOURCES}
        test_json_serializer.cpp
        test_epi_data_io.cpp
        test_lct_parameters_io.cpp
        test_ide_parameters_io.cpp
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
target_link_libraries(memilio-test PRIVATE memilio ode_secir ode_seir ode_secirvvs ode_secirts ode_seair ide_seir ide_secir lct_secir glct_secir abm gtest_main AD::AD dabm smm)
target_compile_options(memilio-test PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

# make unit tests find the test data files
file(TO_CMAKE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data" MEMILIO_TEST_DATA_DIR)
configure_file(test_data_dir.h.in test_data_dir.h)
target_include_directories(memilio-test PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

add_gtest(memilio-test)
