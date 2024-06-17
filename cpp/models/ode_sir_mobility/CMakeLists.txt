add_library(ode_sir_mobility
    infection_state.h
    model.h
    model.cpp
    parameters.h
    regions.h
)
target_link_libraries(ode_sir_mobility PUBLIC memilio)
target_include_directories(ode_sir_mobility PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ode_sir_mobility PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
