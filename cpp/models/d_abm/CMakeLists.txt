add_library(d_abm
    model.h
    model.cpp
    simulation.h
    simulation.cpp
    parameters.h
)
target_link_libraries(d_abm PUBLIC memilio)
target_include_directories(d_abm PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(d_abm PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
