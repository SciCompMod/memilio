add_library(ode_secirts
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
    state_estimators.h
)
target_link_libraries(ode_secirts PUBLIC memilio)
target_include_directories(ode_secirts PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ode_secirts PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
