# configure directory that contains the data files used by examples
file(TO_CMAKE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../data" MEMILIO_DATA_DIR)
configure_file(data_dir.h.in data_dir.h)

add_executable(euler_example euler_test.cpp)
target_link_libraries(euler_example PRIVATE memilio)

add_executable(ode_secir_parameter_sampling_example ode_secir_parameter_sampling.cpp)
target_link_libraries(ode_secir_parameter_sampling_example PRIVATE memilio ode_secir)

add_executable(adapt_rk_example adapt_rk_test.cpp)
target_link_libraries(adapt_rk_example PRIVATE memilio)

add_executable(seir_example seir.cpp)
target_link_libraries(seir_example PRIVATE memilio ode_seir)

add_executable(seair_example seair.cpp)
target_link_libraries(seair_example PRIVATE memilio ode_seair)


add_executable(ode_secir_example ode_secir.cpp)
target_link_libraries(ode_secir_example PRIVATE memilio ode_secir)

add_executable(secirvvs_example secirvvs.cpp)
target_link_libraries(secirvvs_example PRIVATE memilio ode_secirvvs)

add_executable(ode_secir_ageres_example ode_secir_ageres.cpp)
target_link_libraries(ode_secir_ageres_example PRIVATE memilio ode_secir)

add_executable(graph_example graph.cpp)
target_link_libraries(graph_example PRIVATE memilio ode_seir)

add_executable(abm_minimal_example abm_minimal.cpp)
target_link_libraries(abm_minimal_example PRIVATE memilio abm)

if(MEMILIO_HAS_JSONCPP)
  add_executable(ode_secir_read_graph_example ode_secir_read_graph.cpp)
  target_link_libraries(ode_secir_read_graph_example PRIVATE memilio ode_secir)
  target_include_directories(ode_secir_read_graph_example PRIVATE ${CMAKE_CURRENT_BINARY_DIR}) # configured headers
endif()

add_executable(twitter_migration_example twitter_migration.cpp)
target_link_libraries(twitter_migration_example PRIVATE memilio ode_secir)

if(MEMILIO_HAS_HDF5 AND MEMILIO_HAS_JSONCPP)
  add_executable(ode_secir_parameter_study_example ode_secir_parameter_study.cpp)
  target_link_libraries(ode_secir_parameter_study_example PRIVATE memilio ode_secir)
endif()

if(MEMILIO_HAS_JSONCPP)
  add_executable(serialize_example serialize.cpp)
  target_link_libraries(serialize_example PRIVATE memilio)
endif()

add_executable(ide_seir_example ide_seir.cpp)
target_link_libraries(ide_seir_example PRIVATE memilio ide_seir)

if(MEMILIO_HAS_HDF5)
  add_executable(ode_secir_save_results ode_secir_save_results.cpp)
  target_link_libraries(ode_secir_save_results PRIVATE memilio ode_secir)
endif()
