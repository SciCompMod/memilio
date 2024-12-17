if(MEMILIO_HAS_JSONCPP AND MEMILIO_HAS_HDF5)
    add_executable(ide_convergence_rate ide_convergence_rate.cpp)
    target_link_libraries(ide_convergence_rate PRIVATE memilio ode_secir ide_secir Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(ide_convergence_rate PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(ide_changepoints ide_changepoints.cpp)
    target_link_libraries(ide_changepoints PRIVATE memilio ode_secir ide_secir Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(ide_changepoints PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(compute_parameters compute_parameters.cpp)
    target_compile_options(compute_parameters PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
    
    add_executable(ide_covid_scenario ide_covid_scenario.cpp)
    target_link_libraries(ide_covid_scenario PRIVATE memilio ode_secir ide_secir Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(ide_covid_scenario PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
endif()
