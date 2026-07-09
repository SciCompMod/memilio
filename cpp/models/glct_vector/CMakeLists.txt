add_library(glct_vector
    infection_state.h
    model.h
    model.cpp
    parameters.h
)
target_link_libraries(glct_vector PUBLIC memilio)
target_include_directories(glct_vector PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(glct_vector PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})