if (TARGET epidemiology_io)
    add_executable(paper_202011 paper_202011.cpp)
    target_link_libraries(paper_202011 PRIVATE epidemiology epidemiology_io Boost::filesystem)

    add_executable(paper_202107 paper_202107.cpp)
    target_link_libraries(paper_202107 PRIVATE epidemiology epidemiology_io Boost::filesystem)

    add_executable(paper_202107_one_region paper_202107_one_region.cpp)
    target_link_libraries(paper_202107_one_region PRIVATE epidemiology epidemiology_io Boost::filesystem)

endif (TARGET epidemiology_io)
