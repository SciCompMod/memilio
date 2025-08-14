add_library(ode_mseirs4
    infection_state.h
    parameters.h
    model.h
    model.cpp
    README.md
)
target_link_libraries(ode_mseirs4 PUBLIC memilio)
target_include_directories(ode_mseirs4 PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ode_mseirs4 PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
