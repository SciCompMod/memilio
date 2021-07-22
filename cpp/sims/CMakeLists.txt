if (TARGET epidemiology_io)
    add_executable(paper1 paper1.cpp)
    target_link_libraries(paper1 PRIVATE epidemiology epidemiology_io Boost::filesystem)
endif (TARGET epidemiology_io)