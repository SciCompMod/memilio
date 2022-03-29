add_library(secir_vaccine
    infection_state.h
    analyze_result.h
    analyze_result.cpp
    parameter_space.h
    parameter_space.cpp
    secir_params.h
    secir_parameters_io.h
    secir_parameters_io.cpp
    secir.h
    secir.cpp
)
target_link_libraries(secir_vaccine PUBLIC memilio)
target_include_directories(secir_vaccine PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(secir_vaccine PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
