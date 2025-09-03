add_library(ide_endemic_sird
    infection_state.h
    model.h
    model.cpp
    normalized_model.h
    normalized_model.cpp
    simulation.h
    simulation.cpp
    parameters.h
    computed_parameters.h
)
target_link_libraries(ide_endemic_sird PUBLIC memilio)
target_include_directories(ide_endemic_sird PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ide_endemic_secir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
