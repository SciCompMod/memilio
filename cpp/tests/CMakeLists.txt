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
  test_seir_js_compare.cpp
  test_numericalIntegration.cpp
  test_damping.cpp
  test_secir.cpp
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(runUnitTests ${TESTSOURCES})
target_link_libraries(runUnitTests PRIVATE epidemiology gtest_main)
add_gtest(runUnitTests)
