add_library(ide_secir
    infection_state.h
    model.h
    model.cpp
    simulation.h
    simulation.cpp
    parameters.h
    parameters_io.h
)
target_link_libraries(ide_secir PUBLIC memilio)
target_include_directories(ide_secir PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ide_secir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
