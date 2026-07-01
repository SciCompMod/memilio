add_library(ode_seirdb
    infection_state.h
    model.h
    model.cpp
    parameters.h
)
target_link_libraries(ode_seirdb PUBLIC memilio)
target_include_directories(ode_seirdb PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ode_seirdb PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
