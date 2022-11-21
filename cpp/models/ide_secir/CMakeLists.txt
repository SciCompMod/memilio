add_library(ode_secir
    infection_state.h
    model.h
    model.cpp
)
target_link_libraries(ode_secir PUBLIC memilio)
target_include_directories(ode_secir PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ode_secir PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
