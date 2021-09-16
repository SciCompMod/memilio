cmake_minimum_required(VERSION 3.4)
project(epidemiology-tests)

if (NOT TARGET "epidemiology")
  find_package(epidemiology REQUIRED)
endif()

enable_testing()
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(GOOGLE_TEST_INDIVIDUAL ON)
include(AddGoogleTest)

set(TESTSOURCES
  testmain.cpp
  test_populations.cpp
  test_seir_js_compare.cpp
  test_numericalIntegration.cpp
  test_smoother.cpp
  test_damping.cpp
  test_secir.cpp
  test_implicit_euler.cpp
  test_migration.cpp
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
  distributions_helpers.h
  distributions_helpers.cpp
  actions.h
  matchers.h
)

if (TARGET epidemiology_io)
  set(TESTSOURCES ${TESTSOURCES}
    temp_file_register.h
    test_save_parameters.cpp
    test_save_results.cpp
    test_mobility_io.cpp
    test_json_serializer.cpp
  )
endif()

add_executable(runUnitTests ${TESTSOURCES})
target_include_directories(runUnitTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(runUnitTests PRIVATE epidemiology gtest_main)
target_compile_options(runUnitTests PRIVATE
    $<$<CXX_COMPILER_ID:GNU>: -Wall -Wextra -Werror -Wshadow --pedantic-errors -Wno-deprecated-copy>
    $<$<CXX_COMPILER_ID:CLANG>: -Wall -Wextra -Werror -Wshadow --pedantic-errors -Wno-deprecated>
    $<$<CXX_COMPILER_ID:MSVC>: /W4 /WX>)
target_compile_features(runUnitTests PUBLIC cxx_std_14)

if (TARGET epidemiology_io)
  target_link_libraries(runUnitTests PRIVATE epidemiology_io)
endif()

#make unit tests find the test data files
file(TO_CMAKE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data" EPI_TEST_DATA_DIR)
configure_file(test_data_dir.h.in test_data_dir.h)
target_include_directories(runUnitTests PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

add_gtest(runUnitTests)
