add_library(ode_seirs
    infection_state.h
    model.h
    model.cpp
    parameters.h
)
target_link_libraries(ode_seirs PUBLIC memilio)
target_include_directories(ode_seirs PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ode_seirs PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
