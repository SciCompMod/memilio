add_library(ode_metapop_wang
    infection_state.h
    model.h
    model.cpp
    parameters.h
    regions.h
)
target_link_libraries(ode_metapop_wang PUBLIC memilio)
target_include_directories(ode_metapop_wang PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(ode_metapop_wang PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
