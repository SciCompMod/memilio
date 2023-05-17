add_library(graph_abm
    graph_world.cpp
    graph_world.h
)
target_link_libraries(graph_abm PUBLIC memilio)
target_include_directories(graph_abm PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(graph_abm PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
