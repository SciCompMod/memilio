if(MEMILIO_HAS_JSONCPP AND MEMILIO_HAS_HDF5)
    add_executable(2020_npis_wildtype 2020_npis_sarscov2_wildtype_germany.cpp)
    target_link_libraries(2020_npis_wildtype PRIVATE memilio ode_secir Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(2020_npis_wildtype PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(2021_vaccination_delta 2021_vaccination_sarscov2_delta_germany.cpp)
    target_link_libraries(2021_vaccination_delta PRIVATE memilio ode_secirvvs Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(2021_vaccination_delta PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(abm_simulation abm.cpp)
    target_link_libraries(abm_simulation PRIVATE memilio abm Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(abm_simulation PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(munich_graph_sim munich_graph_sim)
    target_link_libraries(munich_graph_sim PRIVATE memilio ode_secir abm Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(munich_graph_sim PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(abm_braunschweig abm_braunschweig.cpp)
    target_link_libraries(abm_braunschweig PRIVATE memilio abm Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(abm_braunschweig PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(ide_convergence_rate ide_convergence_rate.cpp)
    target_link_libraries(ide_convergence_rate PRIVATE memilio ode_secir ide_secir Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(ide_convergence_rate PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(ide_changepoints ide_changepoints.cpp)
    target_link_libraries(ide_changepoints PRIVATE memilio ode_secir ide_secir Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(ide_changepoints PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(ide_covid_scenario ide_covid_scenario.cpp)
    target_link_libraries(ide_covid_scenario PRIVATE memilio ode_secir ide_secir Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(ide_covid_scenario PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
endif()
