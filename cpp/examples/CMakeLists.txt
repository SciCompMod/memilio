# configure directory that contains the data files used by examples
file(TO_CMAKE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../data" MEMILIO_DATA_DIR)
configure_file(data_dir.h.in data_dir.h)

add_executable(euler_example euler_test.cpp)
target_link_libraries(euler_example PRIVATE memilio)
target_compile_options(euler_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ode_secir_parameter_sampling_example ode_secir_parameter_sampling.cpp)
target_link_libraries(ode_secir_parameter_sampling_example PRIVATE memilio ode_secir)
target_compile_options(ode_secir_parameter_sampling_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(adapt_rk_example adapt_rk_test.cpp)
target_link_libraries(adapt_rk_example PRIVATE memilio)
target_compile_options(adapt_rk_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ode_seir_example ode_seir.cpp)
target_link_libraries(ode_seir_example PRIVATE memilio ode_seir)
target_compile_options(ode_seir_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ode_secir_example ode_secir.cpp)
target_link_libraries(ode_secir_example PRIVATE memilio ode_secir)
target_compile_options(ode_secir_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ode_secirvvs_example ode_secirvvs.cpp)
target_link_libraries(ode_secirvvs_example PRIVATE memilio ode_secirvvs)
target_compile_options(ode_secirvvs_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ode_secir_ageres_example ode_secir_ageres.cpp)
target_link_libraries(ode_secir_ageres_example PRIVATE memilio ode_secir)
target_compile_options(ode_secir_ageres_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(graph_example graph.cpp)
target_link_libraries(graph_example PRIVATE memilio ode_seir)
target_compile_options(graph_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(graph_stochastic_mobility_example graph_stochastic_mobility.cpp)
target_link_libraries(graph_stochastic_mobility_example PRIVATE memilio ode_secir)
target_compile_options(graph_stochastic_mobility_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(abm_minimal_example abm_minimal.cpp)
target_link_libraries(abm_minimal_example PRIVATE memilio abm)
target_compile_options(abm_minimal_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(abm_history_example abm_history_object.cpp)
target_link_libraries(abm_history_example PRIVATE memilio abm)
target_compile_options(abm_history_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(twitter_migration_example twitter_migration.cpp)
target_link_libraries(twitter_migration_example PRIVATE memilio ode_secir)

add_executable(ide_seir_example ide_seir.cpp)
target_link_libraries(ide_seir_example PRIVATE memilio ide_seir)
target_compile_options(twitter_migration_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ide_secir_example ide_secir.cpp)
target_link_libraries(ide_secir_example PRIVATE memilio ide_secir)
target_compile_options(ide_secir_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

if(MEMILIO_HAS_JSONCPP)
  add_executable(ode_secir_read_graph_example ode_secir_read_graph.cpp)
  target_link_libraries(ode_secir_read_graph_example PRIVATE memilio ode_secir)
  target_include_directories(ode_secir_read_graph_example PRIVATE ${CMAKE_CURRENT_BINARY_DIR}) # configured headers
  target_compile_options(ode_secir_read_graph_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
endif()

if(MEMILIO_HAS_HDF5 AND MEMILIO_HAS_JSONCPP)
  add_executable(ode_secir_parameter_study_example ode_secir_parameter_study.cpp)
  target_link_libraries(ode_secir_parameter_study_example PRIVATE memilio ode_secir)
  target_compile_options(ode_secir_parameter_study_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
endif()

if(MEMILIO_HAS_JSONCPP)
  add_executable(cli_example cli.cpp)
  target_link_libraries(cli_example PRIVATE memilio)
  target_compile_options(cli_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
endif()

if(MEMILIO_HAS_JSONCPP)
  add_executable(serialize_example serialize.cpp)
  target_link_libraries(serialize_example PRIVATE memilio)
  target_compile_options(serialize_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
endif()

if(MEMILIO_HAS_HDF5)
  add_executable(ode_secir_save_results_example ode_secir_save_results.cpp)
  target_link_libraries(ode_secir_save_results_example PRIVATE memilio ode_secir)
  target_compile_options(ode_secir_save_results_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
endif()
