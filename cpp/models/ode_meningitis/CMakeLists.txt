add_library(ode_meningitis
    analyze_result.h
    analyze_result.cpp
    infection_state.h
    parameter_space.h
    parameter_space.cpp
    parameters.h
    model.h
    model.cpp
)
target_link_libraries(ode_meningitis PUBLIC memilio)
target_include_directories(ode_meningitis PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ode_meningitis PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
