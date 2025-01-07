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

add_executable(ode_seir_ageres_example ode_seir_ageres.cpp)
target_link_libraries(ode_seir_ageres_example PRIVATE memilio ode_seir)
target_compile_options(ode_seir_ageres_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ode_sir_example ode_sir.cpp)
target_link_libraries(ode_sir_example PRIVATE memilio ode_sir)
target_compile_options(ode_sir_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(sde_sir_example sde_sir.cpp)
target_link_libraries(sde_sir_example PRIVATE memilio sde_sir)
target_compile_options(sde_sir_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(sde_sirs_example sde_sirs.cpp)
target_link_libraries(sde_sirs_example PRIVATE memilio sde_sirs)
target_compile_options(sde_sirs_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(sde_seirvv_example sde_seirvv.cpp)
target_link_libraries(sde_seirvv_example PRIVATE memilio sde_seirvv)
target_compile_options(sde_seirvv_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ode_seair_example ode_seair.cpp)
target_link_libraries(ode_seair_example PRIVATE memilio ode_seair AD::AD)
target_compile_options(ode_seair_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

if(MEMILIO_ENABLE_IPOPT)
    add_executable(ode_seair_optimization ode_seair_optimization.cpp)
    target_link_libraries(ode_seair_optimization PRIVATE memilio ode_seair AD::AD Ipopt::Ipopt)
    target_compile_options(ode_seair_optimization PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
endif()

add_executable(ad_square_example ad_square_example.cpp)
target_link_libraries(ad_square_example PRIVATE memilio AD::AD)
target_compile_options(ad_square_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ad_odeint_example ad_odeint_example.cpp)
target_link_libraries(ad_odeint_example PRIVATE memilio AD::AD Boost::boost)
target_compile_options(ad_odeint_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ode_sir_ageres_example ode_sir_ageres.cpp)
target_link_libraries(ode_sir_ageres_example PRIVATE memilio ode_sir)
target_compile_options(ode_sir_ageres_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(seir_flows_example ode_seir_flows.cpp)
target_link_libraries(seir_flows_example PRIVATE memilio ode_seir)
target_compile_options(seir_flows_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ode_secir_example ode_secir.cpp)
target_link_libraries(ode_secir_example PRIVATE memilio ode_secir)
target_compile_options(ode_secir_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ode_secir_contact_changes ode_secir_contact_changes.cpp)
target_link_libraries(ode_secir_contact_changes PRIVATE memilio ode_secir)
target_compile_options(ode_secir_contact_changes PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ode_secirvvs_example ode_secirvvs.cpp)
target_link_libraries(ode_secirvvs_example PRIVATE memilio ode_secirvvs)
target_compile_options(ode_secirvvs_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ode_secirts_example ode_secirts.cpp)
target_link_libraries(ode_secirts_example PRIVATE memilio ode_secirts)
target_compile_options(ode_secirts_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ode_secir_ageres_example ode_secir_ageres.cpp)
target_link_libraries(ode_secir_ageres_example PRIVATE memilio ode_secir)
target_compile_options(ode_secir_ageres_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ode_secir_graph_example ode_secir_graph.cpp)
target_link_libraries(ode_secir_graph_example PRIVATE memilio ode_secir)
target_compile_options(ode_secir_graph_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(graph_stochastic_mobility_example graph_stochastic_mobility.cpp)
target_link_libraries(graph_stochastic_mobility_example PRIVATE memilio ode_secir)
target_compile_options(graph_stochastic_mobility_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(abm_minimal_example abm_minimal.cpp)
target_link_libraries(abm_minimal_example PRIVATE memilio abm)
target_compile_options(abm_minimal_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(abm_history_example abm_history_object.cpp)
target_link_libraries(abm_history_example PRIVATE memilio abm)
target_compile_options(abm_history_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(twitter_mobility_example twitter_mobility.cpp)
target_link_libraries(twitter_mobility_example PRIVATE memilio ode_secir)
target_compile_options(twitter_mobility_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ide_seir_example ide_seir.cpp)
target_link_libraries(ide_seir_example PRIVATE memilio ide_seir)
target_compile_options(ide_seir_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ide_secir_example ide_secir.cpp)
target_link_libraries(ide_secir_example PRIVATE memilio ide_secir)
target_compile_options(ide_secir_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(ide_secir_ageres_example ide_secir_ageres.cpp)
target_link_libraries(ide_secir_ageres_example PRIVATE memilio ide_secir)
target_compile_options(ide_secir_ageres_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(lct_secir_example lct_secir.cpp)
target_link_libraries(lct_secir_example PRIVATE memilio lct_secir)
target_compile_options(lct_secir_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(glct_secir_example glct_secir.cpp)
target_link_libraries(glct_secir_example PRIVATE memilio glct_secir)
target_compile_options(glct_secir_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

add_executable(history_example history.cpp)
target_link_libraries(history_example PRIVATE memilio)
target_compile_options(history_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

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

  add_executable(ode_secir_parameter_study_graph ode_secir_parameter_study_graph.cpp)
  target_link_libraries(ode_secir_parameter_study_graph PRIVATE memilio ode_secir)
  target_compile_options(ode_secir_parameter_study_graph PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
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

if(MEMILIO_HAS_JSONCPP)
    add_executable(ide_initialization_example ide_initialization.cpp)
    target_link_libraries(ide_initialization_example PRIVATE memilio ide_secir)
    target_compile_options(ide_initialization_example PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
endif()
