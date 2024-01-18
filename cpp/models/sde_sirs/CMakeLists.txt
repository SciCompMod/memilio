add_library(sde_sirs
    infection_state.h
    model.h
    model.cpp
    parameters.h
)
target_link_libraries(sde_sirs PUBLIC memilio)
target_include_directories(sde_sirs PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_options(sde_sirs PRIVATE ${MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS})
