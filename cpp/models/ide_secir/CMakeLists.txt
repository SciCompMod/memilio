add_library(ide_secir
    infection_state.h
    model_ide.h
    model_ide.cpp
    simulation.h
    simulation.cpp
    parameters.h
    initialize_from_ode.h
    initialize_from_ode.cpp
    parameters_io.h
    parameters_io.cpp
)
target_link_libraries(ide_secir PUBLIC memilio)
target_include_directories(ide_secir PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ide_secir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
