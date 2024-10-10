add_library(ode_seir_mobility_improved
    infection_state.h
    model.h
    model.cpp
    parameters.h
    regions.h
)
target_link_libraries(ode_seir_mobility_improved PUBLIC memilio)
target_include_directories(ode_seir_mobility_improved PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ode_seir_mobility_improved PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
