cmake_minimum_required(VERSION 3.4)

project(epidemiology-examples)

if (NOT TARGET "epidemiology")
  find_package(epidemiology REQUIRED)
endif()

add_executable(euler_test euler_test.cpp)
target_link_libraries(euler_test PRIVATE epidemiology)

add_executable(rungekutta_test rungekutta_test.cpp)
target_link_libraries(rungekutta_test PRIVATE epidemiology)

add_executable(seir seir.cpp)
target_link_libraries(seir PRIVATE epidemiology)

add_executable(secir secir.cpp)
target_link_libraries(secir PRIVATE epidemiology)

add_executable(group_seir group_seir.cpp)
target_link_libraries(group_seir PRIVATE epidemiology)