add_library(ide_sir_analytical_S_deriv
    infection_state.h
    model.h
    model.cpp
    simulation.h
    simulation.cpp
    parameters.h
    gregory_weights.h 
    gregory_weights.cpp
)
target_link_libraries(ide_sir_analytical_S_deriv PUBLIC memilio)
target_include_directories(ide_sir_analytical_S_deriv PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ide_sir_analytical_S_deriv PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
