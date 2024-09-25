if(MEMILIO_HAS_JSONCPP AND MEMILIO_HAS_HDF5)
    add_executable(2020_npis_wildtype 2020_npis_sarscov2_wildtype_germany.cpp)
    target_link_libraries(2020_npis_wildtype PRIVATE memilio ode_secir Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(2020_npis_wildtype PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(2021_vaccination_delta 2021_vaccination_sarscov2_delta_germany.cpp)
    target_link_libraries(2021_vaccination_delta PRIVATE memilio ode_secirvvs Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(2021_vaccination_delta PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(2022_omicron_late_phase_mobility_vvs 2022_omicron_late_phase_mobility_vvs.cpp)
    target_link_libraries(2022_omicron_late_phase_mobility_vvs PRIVATE memilio ode_secirvvs Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(2022_omicron_late_phase_mobility_vvs PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(2022_omicron_late_phase_mobility_ts 2022_omicron_late_phase_mobility_ts.cpp)
    target_link_libraries(2022_omicron_late_phase_mobility_ts PRIVATE memilio ode_secirts Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(2022_omicron_late_phase_mobility_ts PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(amelag_sim 2022_amelag.cpp)
    target_link_libraries(amelag_sim PRIVATE memilio ode_secirts Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(amelag_sim PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(amelag_single Amelag_single.cpp)
    target_link_libraries(amelag_single PRIVATE memilio ode_secirts Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(amelag_single PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(amelag_single_vvs amelag_single_vvs.cpp)
    target_link_libraries(amelag_single_vvs PRIVATE memilio ode_secirvvs Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(amelag_single_vvs PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(2022_oktober_festival 2022_oktober_festival.cpp)
    target_link_libraries(2022_oktober_festival PRIVATE memilio ode_secirts Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(2022_oktober_festival PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(abm_simulation abm.cpp)
    target_link_libraries(abm_simulation PRIVATE memilio abm Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(abm_simulation PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(munich_graph_sim munich_graph_sim)
    target_link_libraries(munich_graph_sim PRIVATE memilio ode_secir abm Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(munich_graph_sim PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(abm_braunschweig abm_braunschweig.cpp)
    target_link_libraries(abm_braunschweig PRIVATE memilio abm Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(abm_braunschweig PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
endif()
