add_library(ide_sir
    infection_state.h
    model.h
    model.cpp
    model_simplified.h
    model_simplified.cpp
    simulation.h
    simulation.cpp
    simulation_simplified.h
    simulation_simplified.cpp
    parameters.h
    gregory_weights.h 
    gregory_weights.cpp
)
target_link_libraries(ide_sir PUBLIC memilio)
target_include_directories(ide_sir PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ide_sir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
