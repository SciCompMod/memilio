add_library(lct_secir_2_disease
    infection_state.h
    model.h
    model.cpp
    parameters.h
)
target_link_libraries(lct_secir_2_disease PUBLIC memilio)
target_include_directories(lct_secir_2_disease PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(lct_secir_2_disease PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
