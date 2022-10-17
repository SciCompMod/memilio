if(MEMILIO_HAS_JSONCPP AND MEMILIO_HAS_HDF5)
    add_executable(paper_202011 paper_202011.cpp)
    if(MEMILIO_HAS_MPI)
        target_include_directories(paper_202011 PRIVATE ${MPI_CXX_INCLUDE_DIRS})
        target_link_libraries(paper_202011 PRIVATE memilio secir Boost::filesystem ${HDF5_C_LIBRARIES}
                              ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
        target_compile_options(paper_202011 PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS}
                               ${MPI_CXX_COMPILE_FLAGS})
    else()
        target_link_libraries(paper_202011 PRIVATE memilio secir Boost::filesystem ${HDF5_C_LIBRARIES})
        target_compile_options(paper_202011 PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
    endif()

    add_executable(paper_202107 paper_202107.cpp)
    target_link_libraries(paper_202107 PRIVATE memilio ode_secirvvs Boost::filesystem ${HDF5_C_LIBRARIES})
    target_compile_options(paper_202107 PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
endif()
