#configure directory that contains the data files used by examples
file(TO_CMAKE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../data" EPI_DATA_DIR)
configure_file(data_dir.h.in data_dir.h)

add_executable(euler_example euler_test.cpp)
target_link_libraries(euler_example PRIVATE memilio)

add_executable(parameter_sampling_example parameter_sampling.cpp)
target_link_libraries(parameter_sampling_example PRIVATE memilio secir)

add_executable(adapt_rk_example adapt_rk_test.cpp)
target_link_libraries(adapt_rk_example PRIVATE memilio)

add_executable(seir_example seir.cpp)
target_link_libraries(seir_example PRIVATE memilio seir)

add_executable(secir_example secir.cpp)
target_link_libraries(secir_example PRIVATE memilio secir)

add_executable(secir_ageres_example secir_ageres.cpp)
target_link_libraries(secir_ageres_example PRIVATE memilio secir)

add_executable(graph_example graph.cpp)
target_link_libraries(graph_example PRIVATE memilio seir)

add_executable(abm_example abm.cpp)
target_link_libraries(abm_example PRIVATE memilio abm)

if (MEMILIO_HAS_JSONCPP)
  add_executable(read_graph_example read_graph.cpp)
  target_link_libraries(read_graph_example PRIVATE memilio secir)
  target_include_directories(read_graph_example PRIVATE ${CMAKE_CURRENT_BINARY_DIR}) #configured headers
endif()

add_executable(twitter_migration_example twitter_migration.cpp)
target_link_libraries(twitter_migration_example PRIVATE memilio secir)

if (MEMILIO_HAS_HDF5 AND MEMILIO_HAS_JSONCPP)
  add_executable(parameter_study_secir_example parameter_study_secir.cpp)
  target_link_libraries(parameter_study_secir_example PRIVATE memilio secir)
endif()

if (MEMILIO_HAS_JSONCPP)
  add_executable(serialize_example serialize.cpp)
  target_link_libraries(serialize_example PRIVATE memilio)
endif()
