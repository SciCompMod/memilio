add_library(ode_secirvvs
    analyze_result.h
    analyze_result.cpp
    infection_state.h
    parameter_space.h
    parameter_space.cpp
    parameters.h
    parameters_io.h
    parameters_io.cpp
    model.h
    model.cpp
)
target_link_libraries(ode_secirvvs PUBLIC memilio)
target_include_directories(ode_secirvvs PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ode_secirvvs PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
