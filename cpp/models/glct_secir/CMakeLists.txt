add_library(glct_secir
    infection_state.h
    model.h
    model.cpp
    parameters.h
)
target_link_libraries(glct_secir PUBLIC memilio)
target_include_directories(glct_secir PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(glct_secir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})