add_library(ode_secir
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
    reproduction_nb_secir.h
    reproduction_nb_secir.cpp
)
target_link_libraries(ode_secir PUBLIC memilio)
target_include_directories(ode_secir PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ode_secir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
