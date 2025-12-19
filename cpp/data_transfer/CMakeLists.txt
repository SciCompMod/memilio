if(MEMILIO_HAS_JSONCPP AND MEMILIO_HAS_HDF5)
    add_executable(generate_graph_from_data generate_graph_from_data.cpp)
    target_link_libraries(generate_graph_from_data PRIVATE memilio ode_secirvvs Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(generate_graph_from_data PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    # disable warnings in order to compile with gcc 14.1.1
    add_executable(generate_extrapolated generate_extrapolated.cpp)
    target_link_libraries(generate_extrapolated PRIVATE memilio ode_secirvvs Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(generate_extrapolated PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    # disable warnings in order to compile with gcc 14.1.1
    add_executable(simulate_scenario simulate_scenario.cpp)
    target_link_libraries(simulate_scenario PRIVATE memilio ode_secirvvs Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(simulate_scenario PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    # disable warnings in order to compile with gcc 14.1.1
    add_executable(simulate_optimal_control_scenario simulate_optimal_control_scenario.cpp)
    target_link_libraries(simulate_optimal_control_scenario PRIVATE memilio ode_secirvvs Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(simulate_optimal_control_scenario PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(simulate_lha_scenario simulate_lha_scenario.cpp)
    target_link_libraries(simulate_lha_scenario PRIVATE memilio ode_secirvvs ) # Boost::filesystem ${HDF5_C_LIBRARIES}
    target_compile_options(simulate_lha_scenario PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
endif()
