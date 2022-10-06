if(MEMILIO_HAS_JSONCPP AND MEMILIO_HAS_HDF5)
    add_executable(2020_npis_wildtype 2020_npis_wildtype.cpp)
    target_link_libraries(2020_npis_wildtype PRIVATE memilio secir Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(2020_npis_wildtype PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})

    add_executable(2021_vaccination_delta 2021_vaccination_delta.cpp)
    target_link_libraries(2021_vaccination_delta PRIVATE memilio ode_secirvvs Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(2021_vaccination_delta PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
endif()
