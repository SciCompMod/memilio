enable_testing()
include(GoogleTest)

include(memilio-conan)
find_package(GTest 1.10.0 CONFIG REQUIRED)
memilio_target_alternatives(TARGET GTest::gtest ALTERNATIVES gtest::gtest)
memilio_target_alternatives(TARGET GTest::gmock ALTERNATIVES gtest::gmock)

set(TESTSOURCES
  test_household.cpp
  testmain.cpp
  test_populations.cpp
  test_seir_js_compare.cpp
  test_numericalIntegration.cpp
  test_smoother.cpp
  test_damping.cpp
  test_secir.cpp
  test_implicit_euler.cpp
  test_mobility.cpp
  test_date.cpp
  test_eigen_util.cpp
  test_secir_ageres.cpp
  test_parameter_studies.cpp
  test_graph.cpp
  test_graph_simulation.cpp
  test_stl_util.cpp
  test_uncertain.cpp
  test_time_series.cpp
  test_abm.cpp
  test_analyze_result.cpp
  test_contact_matrix.cpp
  test_type_safe.cpp
  test_custom_index_array.cpp
  test_parameter_set.cpp
  test_matrix_shape.cpp
  test_damping_sampling.cpp
  test_dynamic_npis.cpp
  test_regions.cpp
  test_io_framework.cpp
  test_compartmentsimulation.cpp
  test_mobility_io.cpp
  test_transform_iterator.cpp
  distributions_helpers.h
  distributions_helpers.cpp
  actions.h
  matchers.h
  temp_file_register.h
)

if (MEMILIO_HAS_JSONCPP)
    set(TESTSOURCES ${TESTSOURCES}
        test_json_serializer.cpp
    )
endif()
if (MEMILIO_HAS_JSONCPP AND MEMILIO_HAS_HDF5)        
    set(TESTSOURCES ${TESTSOURCES}
        test_save_parameters.cpp
        test_save_results.cpp
    )
endif()

add_executable(memilio-test ${TESTSOURCES})
target_include_directories(memilio-test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(memilio-test PRIVATE memilio secir seir abm GTest::gtest GTest::gmock) 
target_compile_options(memilio-test PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})            
gtest_discover_tests(memilio-test)

#make unit tests find the test data files
file(TO_CMAKE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data" MEMILIO_TEST_DATA_DIR)
configure_file(test_data_dir.h.in test_data_dir.h)
target_include_directories(memilio-test PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
