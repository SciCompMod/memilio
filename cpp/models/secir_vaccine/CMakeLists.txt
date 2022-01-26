add_library(secir_vaccine
    age_group.h
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
    secir_result_io.h
    secir_result_io.cpp
    secir.h
    secir.cpp
)
target_link_libraries(secir_vaccine PUBLIC memilio)
target_include_directories(secir_vaccine PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(secir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
