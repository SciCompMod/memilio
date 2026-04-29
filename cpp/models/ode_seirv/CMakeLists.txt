add_library(ode_seirv
    infection_state.h
    model.h
    model.cpp
    parameters.h
)
target_link_libraries(ode_seirv PUBLIC memilio)
target_include_directories(ode_seirv PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ode_seirv PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
