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
  test_tensor_helpers.cpp
  test_populations.cpp
  test_seir_js_compare.cpp
  test_numericalIntegration.cpp
  test_damping.cpp
  test_secir.cpp
  test_migration.cpp
  test_eigen_util.cpp
  test_secir_ageres.cpp
  test_parameter_studies.cpp
  test_graph.cpp
  test_stl_util.cpp
  test_uncertain.cpp
  distributions_helpers.h
  distributions_helpers.cpp
  gtest_helpers.h
)

if (TARGET epidemiology_io)
  set(TESTSOURCES ${TESTSOURCES}
    test_save_parameters.cpp
    test_save_results.cpp
    test_twitter.cpp
  )
endif()

add_executable(runUnitTests ${TESTSOURCES})
target_include_directories(runUnitTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(runUnitTests PRIVATE epidemiology gtest_main)

if (TARGET epidemiology_io)
  target_link_libraries(runUnitTests PRIVATE epidemiology_io)
endif()

#make unit tests find the test data files
file(TO_CMAKE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data" EPI_TEST_DATA_DIR)
configure_file(test_data_dir.h.in test_data_dir.h)
target_include_directories(runUnitTests PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

add_gtest(runUnitTests)
