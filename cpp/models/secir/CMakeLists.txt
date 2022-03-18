add_library(secir
    infection_state.h
    analyze_result.h
    analyze_result.cpp
    implicit_euler.h
    implicit_euler.cpp
    parameter_space.h
    parameter_space.cpp
    parameter_studies.h
    secir_params.h
    secir_parameters_io.h
    secir_parameters_io.cpp
    secir.h
    secir.cpp
)
target_link_libraries(secir PUBLIC memilio)
target_include_directories(secir PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(secir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})