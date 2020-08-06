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

add_executable(secir_ageres secir_ageres.cpp)
target_link_libraries(secir_ageres PRIVATE epidemiology)

add_executable(graph graph.cpp)
target_link_libraries(graph PRIVATE epidemiology)

add_executable(abm abm.cpp)
target_link_libraries(abm PRIVATE epidemiology)

if (TARGET "epidemiology_io")
  add_executable(read_graph read_graph.cpp)
  target_link_libraries(read_graph PRIVATE epidemiology epidemiology_io)

  target_link_libraries(secir_ageres PRIVATE epidemiology_io)

  add_executable(parameter_study_secir parameter_study_secir.cpp)
  target_link_libraries(parameter_study_secir PRIVATE epidemiology epidemiology_io)
endif()

#add_executable(group_seir group_seir.cpp)
#target_link_libraries(group_seir PRIVATE epidemiology)
