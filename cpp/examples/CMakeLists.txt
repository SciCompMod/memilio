cmake_minimum_required(VERSION 3.4)

project(epidemiology-examples)

if (NOT TARGET "epidemiology")
  find_package(epidemiology REQUIRED)
endif()

add_executable(euler_test euler_test.cpp)
target_link_libraries(euler_test PRIVATE epidemiology)

add_executable(parameter_sampling parameter_sampling.cpp)
target_link_libraries(parameter_sampling PRIVATE epidemiology)

add_executable(adapt_rk_test adapt_rk_test.cpp)
target_link_libraries(adapt_rk_test PRIVATE epidemiology)

add_executable(seir seir.cpp)
target_link_libraries(seir PRIVATE epidemiology)

add_executable(secir secir.cpp)
target_link_libraries(secir PRIVATE epidemiology)

add_executable(abm_from_statistical_data abm_from_statistical_data.cpp)
target_link_libraries(abm_from_statistical_data PRIVATE epidemiology) 

add_executable(secir_ageres secir_ageres.cpp)
target_link_libraries(secir_ageres PRIVATE epidemiology)

add_executable(graph graph.cpp)
target_link_libraries(graph PRIVATE epidemiology)

add_executable(abm abm.cpp)
target_link_libraries(abm PRIVATE epidemiology)

if (TARGET "epidemiology_io")

  file(TO_CMAKE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../data" EPI_DATA_DIR)
  configure_file(data_dir.h.in data_dir.h)

  add_executable(read_graph read_graph.cpp)
  target_link_libraries(read_graph PRIVATE epidemiology epidemiology_io)
  target_include_directories(read_graph PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

  target_link_libraries(secir_ageres PRIVATE epidemiology_io)

  add_executable(twitter_migration twitter_migration.cpp)
  target_link_libraries(twitter_migration PRIVATE epidemiology epidemiology_io)

  add_executable(parameter_study_secir parameter_study_secir.cpp)
  target_link_libraries(parameter_study_secir PRIVATE epidemiology epidemiology_io)

  add_executable(serialize serialize.cpp)
  target_link_libraries(serialize PRIVATE epidemiology epidemiology_io)
endif()

#add_executable(group_seir group_seir.cpp)
#target_link_libraries(group_seir PRIVATE epidemiology)
