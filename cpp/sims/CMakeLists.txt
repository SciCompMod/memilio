if (TARGET epidemiology_io)
    add_executable(paper_202011 paper_202011.cpp)
    target_link_libraries(paper_202011 PRIVATE epidemiology epidemiology_io Boost::filesystem)
endif (TARGET epidemiology_io)