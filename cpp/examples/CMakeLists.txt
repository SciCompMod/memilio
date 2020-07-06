cmake_minimum_required(VERSION 3.4)

project(epidemiology-examples)

if (NOT TARGET "epidemiology")
  find_package(epidemiology REQUIRED)
endif()
if (NOT TARGET "epidemiology_io")
  find_package(epidemiology_io REQUIRED)
endif()

add_executable(euler_test euler_test.cpp)
target_link_libraries(euler_test PRIVATE epidemiology)

add_executable(rungekutta_test rungekutta_test.cpp)
target_link_libraries(rungekutta_test PRIVATE epidemiology)

add_executable(seir seir.cpp)
target_link_libraries(seir PRIVATE epidemiology)

add_executable(secir secir.cpp)
target_link_libraries(secir PRIVATE epidemiology)

add_executable(secir_ageres secir_ageres.cpp)
target_link_libraries(secir_ageres PRIVATE epidemiology)

add_executable(secir_ensemble secir_ensemble.cpp)
target_link_libraries(secir_ensemble PRIVATE epidemiology epidemiology_io)

add_executable(group_seir group_seir.cpp)
target_link_libraries(group_seir PRIVATE epidemiology)
