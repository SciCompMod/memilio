if(MEMILIO_HAS_JSONCPP AND MEMILIO_HAS_HDF5)

    find_package(LIKWID REQUIRED)
    add_executable(2020_npis_wildtype 2020_npis_sarscov2_wildtype_germany.cpp)
    target_include_directories(2020_npis_wildtype PRIVATE ${LIKWID_INCLUDE_DIRS})
    target_link_libraries(2020_npis_wildtype PRIVATE memilio ode_secir Boost::filesystem ${HDF5_C_LIBRARIES}${LIKWID_LIBRARIES})
    target_compile_options(2020_npis_wildtype PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
    target_compile_definitions(2020_npis_wildtype PRIVATE "-DLIKWID_PERFMON")         

    #find_package(LIKWID REQUIRED)
    add_executable(2021_vaccination_delta 2021_vaccination_sarscov2_delta_germany.cpp)
    target_include_directories(2021_vaccination_delta PRIVATE ${LIKWID_INCLUDE_DIRS})
    target_link_libraries(2021_vaccination_delta PRIVATE memilio ode_secirvvs Boost::filesystem ${HDF5_C_LIBRARIES}${LIKWID_LIBRARIES})
    target_compile_options(2021_vaccination_delta PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
    target_compile_definitions(2021_vaccination_delta PRIVATE "-DLIKWID_PERFMON")         

    #find_package(LIKWID REQUIRED)
    add_executable(abm_simulation abm.cpp)
    target_include_directories(abm_simulation PRIVATE ${LIKWID_INCLUDE_DIRS})
    target_link_libraries(abm_simulation PRIVATE memilio abm Boost::filesystem ${HDF5_C_LIBRARIES} ${LIKWID_LIBRARIES})
    target_compile_options(abm_simulation PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
    target_compile_definitions(abm_simulation PRIVATE "-DLIKWID_PERFMON")         

endif()

